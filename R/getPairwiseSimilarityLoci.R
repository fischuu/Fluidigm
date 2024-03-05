#' @title Get Pairwise Similarity Loci
#'
#' @description
#' This function is a wrapper to a perl script that determines the loci of pairwise similarities. It is designed to work with data from the 'PLINK' software, which is commonly used in bioinformatics for whole genome association analysis.
#'
#' @param file A string. This is the path to the previously created PLINK ped file (with .ped extension). The ped file contains genotype information in a format that can be used for further analysis.
#' @param verbose A logical. If TRUE, the function will print detailed messages during its execution to help you understand what it's doing at each step. Default is TRUE.
#' @param verbosity An integer. This parameter controls the level of verbosity. The higher the number, the more detailed the messages. Default is 1.
#'
#' @details
#' This function first checks the input parameters. It then constructs the output file names and the command to run the perl script.
#' The command is executed using the system function. If the 'verbose' parameter is set to TRUE, the function will print a message
#' when it has finished running.
#'
#' The perl script was written by Doug Scofield, see references.
#'
#' @references
#' The original code this function is based on can be found at: GitHub https://github.com/douglasgscofield/bioinfo/blob/main/scripts/plink-pairwise-loci.pl
#'
#' @examples
#' \dontrun{
#'   outdir <- tempdir()
#'   getPairwiseSimilarityLoci(file = file.path(outdir, "example_data.csv.GOOD"))
#' }
#'
#' @return This function does not return a value in the R environment. Instead, it creates an output file with the
#' '.pairs' extension in the same directory as the input file. This output file contains the results of the pairwise
#' similarity loci analysis.
#'
#' @export

getPairwiseSimilarityLoci <- function(file, verbose=TRUE, verbosity=1){

  ## Verbose output of input parameters
  #################################################

  if (verbosity >= 2) {
    cat("Input parameters:\n")
    cat("file: ", file, "\n")
    cat("verbose: ", verbose, "\n")
    cat("verbosity: ", verbosity, "\n")
  }

  ### Input check
  ##############################################################################
  # Check if 'file' is a character string
  if(!is.character(file) || length(file) != 1){
    stop("'file' must be a single character string representing the path to the .ped file.")
  }

  # Check if 'verbose' is logical
  if(!is.logical(verbose) || length(verbose) != 1){
    stop("'verbose' must be a single logical value.")
  }

  # Check if 'verbosity' is numeric
  if(!is.numeric(verbosity) || length(verbosity) != 1){
    stop("'verbosity' must be a single numeric value.")
  }

  # Convert 'verbose' to numeric if it's logical
  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)
  verbose <- verbosity

  # Set 'intern.param' based on 'verbose'
  ifelse(verbose>1, intern.param <- FALSE, intern.param <- TRUE)

  # Construct output file names
  out <- paste0(file,".pairs")
  pedfile <- paste0(file, ".ped")
  mibsfile <- paste0(file,".mibs")

  # Find the package and script directories
#  package.dir <- find.package('Fluidigm')
  package.dir <- system.file(package = "Fluidigm", lib.loc = .libPaths()[1])
  perl.dir <- file.path(package.dir ,'perl')
  script <- file.path(perl.dir,'pairwise-loci.pl')

  # Construct the command to run the perl script
  perlCommand <- paste0("perl ",script," --ped ",pedfile," --mibs ",mibsfile," > ", out)

  # Print the command if 'verbose' is greater than 1
  if(verbose>1){
    message("Run the following PERL command:\n", perlCommand)
  }

  # Run the command
    system(perlCommand, ignore.stdout=FALSE, ignore.stderr=intern.param)

  # Print a completion message if 'verbose' is greater than 0
  if(verbose>0)message("### Get pairwise similarity loci: DONE! ",date(),"\n","##############################################################\n")
}
