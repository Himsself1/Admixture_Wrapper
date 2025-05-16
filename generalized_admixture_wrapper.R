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
           pattern = input_params$prefix, full.names = T)
if( length(input_files) < 3 ){
  cat(sprintf("Input files missing or incomplete. \n"))
  stop("Stop 2")
}else{
  initial_file <- grep("bed$|ped$", input_files, value = TRUE)
  ## print(initial_file)
  initial_prefix <- gsub("\\.bed|\\.ped", "", initial_file)
  ## This may be wrong when prefix is matched by more than one triplet.
  if( length(initial_file) > 1 ){
    cat(sprintf("More than one inputs match the prefix. Please rename. \n"))
    stop("Stop 3")
  }
}

if( input_params$project_excluded == TRUE ){
  if( input_params$family_file < 1 ){
    cat(sprintf("Projection of excluded samples was requested but no file was provided"))
    cat(sprintf("Please include a 'family_file'"))
    stop("Stop 4")
  }
}
  
# ** Create output folders
initial_dir <- getwd()
out_dir_full_name <- file.path(input_params$output_folder, input_params$run_name)
out_dir_for_data <- file.path(out_dir_full_name, "data")
out_dir_for_plots <- file.path(out_dir_full_name, "plots") 
out_dir_for_stats <- file.path(out_dir_full_name, "stats")
dir.create(out_dir_for_plots, recursive = TRUE)
dir.create(out_dir_for_stats, recursive = TRUE)
dir.create(out_dir_for_data, recursive = TRUE)
if(input_params$project_excluded == TRUE ){
  out_dir_for_excluded <- file.path(out_dir_for_data, "excluded")
  dir.create(out_dir_for_excluded)
}

# * Call Plink for trimming

# ** Remove Relatives if file is provided

if( input_params$family_file > 1 ){
  no_family_prefix <- file.path(
    out_dir_for_data,
    paste0(c(
      input_params$prefix,
      "_no_family"), collapse = ''
      )
  )
  
  command_for_plink_relatives <- paste0(
    c(
      "plink2 --bfile", initial_prefix,
      "--remove", input_params$family_file,
      "--out", no_family_prefix,
      "--chr 1-22 --make-bed --allow-no-sex --keep-allele-order --no-pheno --set-all-var-ids '@_#_$r_$a' "
    ),
    collapse = ' '
  )
  print(command_for_plink_relatives)
  system(command_for_plink_relatives)

# *** Build dataset of excluded individuals for projection if requested

  if( input_params$project_excluded == TRUE ){
    excluded_prefix <- file.path(
      out_dir_for_excluded,
      paste0(c(
        input_params$prefix,
        "_excluded"), collapse = ''
        )
    )

    command_for_plink_excluded <- paste0(
      c(
        "plink2 --bfile", initial_prefix,
        "--keep", input_params$family_file,
        "--out", excluded_prefix,
        "--chr 1-22 --make-bed --allow-no-sex --keep-allele-order --no-pheno --set-all-var-ids '@_#_$r_$a' "
      ),
      collapse = ' '
    )
    print(command_for_plink_excluded)
    system(command_for_plink_excluded)
  }
  
} else {
  no_family_prefix <- file.path(out_dir_for_data, input_params$prefix)
}

# ** Filter SNPs in LD

filter_prefix <- file.path(
  out_dir_for_data,
  paste0(c(
    input_params$prefix,
    "_filtered"), collapse = ''
    )
)

if( input_params$ld_prune == TRUE ){
## This command outputs a file containing the SNPs that pass the filtering.
  command_for_plink_filtering <- paste0(
    c( 
      "plink2 --bfile", no_family_prefix,
      "--indep-pairwise", input_params$ld_window, input_params$ld_step, input_params$ld_r,
      "--maf 0.05",
      "--geno", input_params$geno,
      "--out", filter_prefix,
      "--allow-no-sex --keep-allele-order --no-pheno"
    ),
    collapse = ' '
  )
  print(command_for_plink_filtering)
  system(command_for_plink_filtering)
  pruned_in <- grep(
    "prune.in$",
    list.files(out_dir_for_data, pattern = basename(filter_prefix), full.names = T),
    value = T
  )
  ## Get the name of the pruned_in snp file
  trimmed_prefix <- file.path(
    out_dir_for_data,
    paste0(c(input_params$prefix,
             "_trimmed"), collapse = '')
  )
  ## This command generates a new .bed file that contains only the SNPs that pass the filtering.
  command_for_plink_trimming <- paste0(
    c(
      "plink2 --bfile", no_family_prefix,
      "--extract", pruned_in, "--make-bed",
      "--out", trimmed_prefix,
      "--allow-no-sex --keep-allele-order --no-pheno"
    ),
    collapse = " "
  )
  print(command_for_plink_trimming)
  system(command_for_plink_trimming)
  # These are the names of input files for ADMIXTURE and plots.
  fam_file <- grep(
    ".fam$",
    list.files(out_dir_for_data, pattern = basename(trimmed_prefix), full.names = T),
    value = T
  )
  bim_file <- grep(
    ".bim$",
    list.files(out_dir_for_data, pattern = basename(trimmed_prefix), full.names = T),
    value = T
  )
  bed_file <- grep(
    ".bed$",
    list.files(out_dir_for_data, pattern = basename(trimmed_prefix), full.names = T),
    value = T
  )

  ## Filter the 'exluded' dataset.
  if( input_params$project_excluded = TRUE ){

    excluded_trimmed_prefix <- file.path(
      out_dir_for_excluded,
      paste0(c(input_params$prefix,
               "_excluded_trimmed"), collapse = '')
    )
    ## Filter the 'excluded' dataset using the SNPs from the original.
    command_for_plink_trimming_for_excluded_individuals <- paste0(
    c(
      "plink2 --bfile", excluded_prefix,
      "--extract", pruned_in, "--make-bed",
      "--out", excluded_trimmed_prefix,
      "--allow-no-sex --keep-allele-order --no-pheno"
    ),
    collapse = " "
  )
    print(command_for_plink_trimming_for_excluded_individuals)
    system(command_for_plink_trimming_for_excluded_individuals)

    fam_file <- grep(
      ".fam$",
      list.files(out_dir_for_excluded, pattern = basename(excluded_trimmed_prefix), full.names = T),
      value = T
    )
    bim_file <- grep(
      ".bim$",
      list.files(out_dir_for_excluded, pattern = basename(excluded_trimmed_prefix), full.names = T),
      value = T
    )
    bed_file <- grep(
      ".bed$",
      list.files(out_dir_for_excluded, pattern = basename(excluded_trimmed_prefix), full.names = T),
      value = T
    )
  }
} else {
  # These are the names of input files for ADMIXTURE and plots if no prunning was performed
  trimmed_prefix <- no_family_prefix
  fam_file <- grep(
    ".fam$",
    list.files(out_dir_for_data, pattern = basename(no_family_prefix), full.names = T),
    value = T
  )
  bim_file <- grep(
    ".bim$",
    list.files(out_dir_for_data, pattern = basename(no_family_prefix), full.names = T),
    value = T
  )
  bed_file <- grep(
    ".bed$",
    list.files(out_dir_for_data, pattern = basename(no_family_prefix), full.names = T),
    value = T
  )
}


# * Running ADMIXTURE

registerDoMC(input_params$threads)

# ** Create input and output names

cv_error_file <- file.path(
  out_dir_for_stats,
  paste0(c(input_params$name, "_cv_errors.csv"), collapse = "" )
)

admixture_output_names <- c()
for(i in 2:input_params$max_K){
  temp_output <- paste0(c(
    out_dir_for_stats,
    input_params$name,
    "K_", i, ".out"
  ), collapse = '')
  admixture_output_names <- c(admixture_output_names, temp_output)
}

# ** System calls

setwd( out_dir_for_stats )
sporos <- runif(input_params$max_K)

foreach(i = 2:input_params$max_K) %dopar% {
  cmd <- noquote(paste0(c(
    "admixture32 --cv",
    bed_file, i,
    '--seed time | tee',
    admixture_output_names[i-1]
  ), collapse = " "))
  print(cmd)
  system(cmd)
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
  print(cmd)
  system(cmd)
}

# * Plotting 

setwd(initial_dir)

plot_cv_error_command <- paste0(c(
  "Rscript plot_CV_error.R",
  "-plot_folder", out_dir_for_plots,
  "-input_file", cv_error_file,
  "-name", input_params$run_name
), collapse = ' ')

print(plot_cv_error_command)
system(plot_cv_error_command)

## Make command for plotting script

plotting_command <- paste0(c(
  "Rscript generalized_plotting.R -input_folder", out_dir_for_stats,
  "-plot_folder", out_dir_for_plots,
  "-prefix", trimmed_prefix,
  "-label_file", fam_file,
  "-meta", input_params$metadata_file,
  "-name", input_params$run_name
), collapse = ' ')

print(plotting_command)
system(plotting_command)
