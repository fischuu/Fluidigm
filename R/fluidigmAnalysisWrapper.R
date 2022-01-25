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
#' @param remove.y logical, Remove markers located on y-chromosome
#' @param plots logical, shall plots be created?
#' @param verbose Should the output be verbose, logical or numerical
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

fluidigmAnalysisWrapper <- function(file, out=NA, db, ymap="PlateD_withY.map", woymap,  keep.rep=2, remove.y=TRUE, plots=TRUE, rearrange=FALSE, group=NA, verbose=TRUE){

  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)

  filename <- basename(file)
  dirname <- dirname(file)
  pedfile <- gsub("csv$","ped",filename)
  plinkfile <- gsub("csv$","GOOD",filename)

  fluidigm2PLINK(file=file, out=out, ymap=ymap, verbose=verbose, plots=plots, rearrange=rearrange)

  out <- estimateErrors(file.path(dirname,pedfile), keep.rep=keep.rep, remove.y=remove.y, plots=plots, verbose=verbose)

  calculatePairwiseSimilarities(file=file.path(dirname,plinkfile), db=db, map=woymap)

  getPairwiseSimilarityLoci(file=file.path(dirname,plinkfile))

  similarityMatrix(file=paste0(file.path(dirname,plinkfile),"_oDB"),
                   group=group)

  if(verbose>0)cat("\n ### All DONE!",date(),"\n","##############################################################\n")

  out
}
