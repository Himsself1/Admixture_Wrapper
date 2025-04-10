# * Libraries

list_of_packages <- c(
  "yaml", "foreach", "doMC"
)

for (i in list_of_packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = T)
  }
}

# * Command line arguments

args<-commandArgs(TRUE)
yaml_input <- args[1]

# ** Verify arguments

if( file.exists( yaml_input ) ){
  input_params <- read_yaml( yaml_input )
}else{
  cat(sprintf( "YAML file: %s not found. \n", yaml_input ))
  stop("Stop 1")
}

input_files <- list.files(path = input_params$path_to_data,
           pattern = input_params$prefix)
if( length(input_files) < 3 ){
  cat(sprintf("Input files missing or incomplete. \n"))
  stop("Stop 2")
}else{
  initial_file <- grep( "bed|ped|geno", input_files, value = TRUE )
  ## This may be wrong when prefix is matched by more than one triplet.
  if( length(input_for_admixture) > 1 ){
    cat(sprintf("More than one inputs match the prefix. Please rename. \n"))
    stop("Stop 3")
  }
}

# ** Create output folders

out_dir_full_name <- paste0(c(input_params$output_folder, input_params$run_name), collapse = '/')
out_dir_for_data <- paste0(c(out_dir_full_name, "/data/"), collapse = '')
out_dir_for_plots <- paste0(c(out_dir_full_name, "/plots/"), collapse = '')
out_dir_for_stats <- paste0(c(out_dir_full_name, "/stats/"), collapse = '')
dir.create( out_dir_for_plots, recursive = TRUE )
dir.create( out_dir_for_stats, recursive = TRUE )

# * Call Plink for trimming

## Need to modify this in order to remove relatives! ##

# ** Remove Relatives

no_family_prefix <- paste0(
  c(outdir_for_data, "/",
    input_params$prefix,
    "_no_family"), collapse = ''
)

command_for_plink_relatives <- paste0(
  c(
    "plink2.0 --bfile", initial_file,
    "--remove", input_params$family_file,
    "--out", no_family_prefix,
    "--allow-no-sex --keep-allele-order --no-pheno"
  ),
  collapse = ' '
)
print(command_for_plink_relatives)
## system(command_for_plink_relatives)

input_for_filtering <- grep(
  "bed|ped",
  list.files(output_for_data, pattern = no_family_prefix, full_names = T),
  value = T
)

# ** Filter SNPs in LD

filter_prefix <- paste0(
  c(outdir_for_data, "/",
    input_params$prefix,
    "_filtered"), collapse = ''
)

## This command outputs a file containing the SNPs that pass the filtering.
command_for_plink_filtering <- paste0(
  c( 
    "plink2.0 --bfile", input_for_filtering,
    "--indep-pairwise 200 25 08", "--maf 0.05",
    "--out", filter_prefix,
    "--allow-no-sex --keep-allele-order --no-pheno"
  ), collapse = ' '
)
print(command_for_plink_filtering)
## system(command_for_plink_filtering)

# ** Trim SNPs

## Get the name of the pruned_in snp file
pruned_in <- grep(
  "prune.in",
  list.files(output_for_data, pattern = filter_prefix, full_names = T),
  value = T
)

trimmed_prefix <- paste0(
  c(outdir_for_data, "/",
    input_params$prefix,
    "_trimmed"), collapse = ''
)

## This command generates a new .bed file that contains only the SNPs that pass the filtering.
command_for_plink_trimming <- paste0(
  c(
    "plink2.0 --bfile", input_for_trimming,
    "--extract", pruned_in, "--make-bed",
    "--out", trimmed_prefix,
    "--allow-no-sex --keep-allele-order --no-pheno"
  ),
  collapse = " "
)

print(command_for_plink_trimming)
## system(command_for_plink_trimming)

# * Running ADMIXTURE

registerDoMC(input_params$threads)

# ** Create input and output names

cv_error_file <- paste0( c(out_dir_for_stats, "/", input_params$name, "cv_errors.csv"), collapse = "" )

admixture_output_names <- c()
for(i in 2:input_params$cv_error){
  temp_output <- paste0(c(
    out_dir_for_stats, "/",
    input_params$name,
    "K_", i, ".haploid.out"
  ), collapse = '')
  admixture_output_names <- c(admixture_output_names, temp_output)
}

# ** System calls

setwd( out_dir_for_stats )

foreach(i = 2:input_params$cv_error) %dopar% {
  cmd <- noquote(paste0(c(
    "admixture32 --cv ",
    input_for_admixture,
    ' --haploid="*" --seed time | tee ',
    admixture_output_names[i-1]
  ), collapse = ""))
  ## system(cmd)
  print(cmd)
}

## CV Errors
for (i in 1:length(admixture_output_names)) {
  cmd <- noquote(paste0(c(
    "grep -h CV ",
    admixture_output_names[i],
    " | awk '{print $3 $4}' | sed 's/(K=\\([0-9]*\\)):\\([0-9]\\.[0-9]*\\)/\\1,\\2/g'",
    ## This regex builds a simple csv with 2 columns: K,CV_error
    ## R needs to escape the escape character ('\') in order to include it in a system command.
    " >> ", cv_error_file
  ), collapse = ""))
  ## system(cmd)
  print(cmd)
}

# * Plotting

