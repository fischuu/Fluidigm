#' @title Get Pairwise Similarity Loci
#'
#' @description
#' This function is a wrapper to a perl script to determine the loci of pairwise similarities
#'
#' @param file Path to the previous plink created ped file (with ped extension)
#' @param verbose Should the output be verbose, logical
#' @param verbosity Level of verbosity, set to higher number for more details
#'
#' #' @details
#' Additional details...
#'
#' @examples
#' \dontrun{
#'   getPairwiseSimilarityLoci()
#' }
#'
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

getPairwiseSimilarityLoci <- function(file, verbose=TRUE, verbosity=1){

  ### Input check
  ##############################################################################
  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)
  verbose <- verbosity
  ifelse(verbose>1, intern.param <- FALSE, intern.param <- TRUE)

  out <- paste0(file,".pairs")
  pedfile <- paste0(file, ".ped")
  mibsfile <- paste0(file,".mibs")

  package.dir <- find.package('Fluidigm')
  perl.dir <- file.path(package.dir,'perl')
  script <- file.path(perl.dir,'pairwise-loci.pl')

  perlCommand <- paste0("perl ",script," --ped ",pedfile," --mibs ",mibsfile," > ", out)

  if(verbose>1){
    cat("Run the following PERL command:\n", perlCommand)
  }
  # Run Plink
    system(perlCommand, intern=intern.param)

  if(verbose>0)cat("\n ### Get pairwise similarity loci: DONE! ",date(),"\n","##############################################################\n")
}
