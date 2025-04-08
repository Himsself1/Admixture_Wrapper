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
if( length(eig_files) < 3 ){
  cat(sprintf("Input files missing or incomplete. \n"))
  stop("Stop 2")
}else{
  input_for_admixture <- grep( "bed|ped|geno", input_files, value = TRUE )
  ## This may be wrong when prefix is matched by more than one triplet.
  if( length(input_for_admixture) > 1 ){
    cat(sprintf("More than one inputs match the prefix. Please rename. \n"))
    stop("Stop 3")
  }
  
}

# ** Create output folders

out_dir_full_name <- paste0( c(input_params$output_folder, input_params$run_name), collapse = '/')
out_dir_for_plots <- paste0(out_dir_full_name, "/plots/", collapse = '')
out_dir_for_stats <- paste0(out_dir_full_name, "/stats/", collapse = '')
dir.create( out_dir_for_plots, recursive = TRUE )
dir.create( out_dir_for_stats, recursive = TRUE )

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
  system(cmd)
  ## print(cmd)
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
  system(cmd)
}

# * Plotting
