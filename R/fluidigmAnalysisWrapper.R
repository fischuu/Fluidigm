#' Run the fluidigm analysis script together
#'
#' This function is a wrapper for the whole analysis
#'
#' @param file Path to the fluidigm input file
#' @param out Out file name, keep empty to keep the original basename
#' @param db Filepath to database file
#' @param ymap Filepath to PlateD_withY.map file
#' @param ymap Filepath to PlateD_withoutY.map file
#' @param keep.rep numeric, keep only this n-fold replicates, default n=2
#' @param y.marker logical, Remove markers located on y-chromosome
#' @param plots logical, shall plots be created?
#' @param verbose Should the output be verbose, logical or numerical
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

fluidigmAnalysisWrapper <- function(file, out=NA, db, ymap=NA, woymap,  keep.rep=2, neg_controls=NA, y.marker=NA, x.marker=NA, plots=TRUE, rearrange=FALSE, group=NA, verbose=TRUE){

  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)

  if(is.na(y.marker)) stop("You do not provide a vector with marker names to y.marker!")

  filename <- basename(file)
  dirname <- dirname(file)
  pedfile <- gsub("csv$","ped",filename)
  plinkfile <- gsub("csv$","GOOD",filename)

  fluidigm2PLINK(file=file, out=out, ymap=ymap, verbose=verbose, plots=plots, rearrange=rearrange)

  out <- estimateErrors(file=file.path(dirname,pedfile), keep.rep=keep.rep, y.marker=y.marker, x.marker=x.marker, neg_controls=neg_controls, plots=plots, verbose=verbose)

  calculatePairwiseSimilarities(file=file.path(dirname,plinkfile), db=db, map=woymap, verbose=verbose)

  getPairwiseSimilarityLoci(file=file.path(dirname,plinkfile), verbose=verbose)

  similarityMatrix(file=paste0(file.path(dirname,plinkfile),"_oDB"),
                   group=group,
                   verbose=verbose)

  if(verbose>0)cat("\n ### All DONE!",date(),"\n","##############################################################\n")

  out
}
