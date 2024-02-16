#' @title Run the fluidigm analysis script together
#'
#' @description
#' This function is a wrapper for the whole analysis
#'
#' @param file Path to the fluidigm input file
#' @param out Out file name, keep empty to keep the original basename
#' @param db Filepath to database file
#' @param appendSamplesToDB Logical, shoud new samples be added to database
#' @param map Filepath to PlateDnoY.map file
#' @param keep.rep numeric, keep only this n-fold replicates, default n=2
#' @param neg_controls Names of negative controls
#' @param y.marker Y markers for sexing
#' @param x.marker X markers for sexing
#' @param sp.marker Markers used for species-identification
#' @param plots logical, shall plots be created?
#' @param allele_error Threshold for RERUN on Allele errors
#' @param marker_dropout Threshold for RERUN on Marker dropout
#' @param no_marker Number of markers
#' @param male.y Threshold for sexing, male y-chromosome markers
#' @param male.hetX Threshold for sexing, heterozygote x-chr markers
#' @param female.y Threshold for sexing, female y-chromosome markers
#' @param female.Xtot Threshold for sexing, total female x-chr markers
#' @param female.hetXtot Threshold for sexing, heterozygote x-chr markers
#' @param warning.noYtot  Threshold for sexing, when should warning be triggered
#' @param warning.noHetXtot  Threshold for sexing, when should warning be triggered
#' @param rearrange Logical, rearrange the ped/map output in order of ymap
#' @param group Sample identified for sample-wise statistics
#' @param fixNames logical, remove whitespaces from sample names automatically
#' @param sexing Logical, if sexing should be performed
#' @param verbose Should the output be verbose, logical or numerical
#' @param verbosity Level of verbosity, set to higher number for more details
#'
#' #' @details
#' Additional details...
#'
#' @examples
#' \dontrun{
#'   fluidigmAnalysisWrapper()
#' }
#'
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

fluidigmAnalysisWrapper <- function(file, out=NA, db=NA, appendSamplesToDB=FALSE, map=NA, keep.rep=1,
                                    neg_controls=NA, y.marker=NA, x.marker=NA, sp.marker=NA, plots=TRUE,
                                    allele_error=5, marker_dropout=15, no_marker=50,
                                    male.y=3, male.hetX=0, female.y=0, female.Xtot=8, female.hetXtot=3,
                                    warning.noYtot=2, warning.noHetXtot=3,
                                    rearrange=TRUE, group=NA, fixNames=TRUE, sexing=TRUE, verbose=TRUE, verbosity=1){

  if(!verbose & verbosity > 0) verbosity <- 0
  verbose <- verbosity
  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)

  if(sexing) if(is.na(y.marker)) stop("You do not provide a vector with marker names to y.marker!")

  filename <- basename(file)
  dirname <- dirname(file)
  pedfile <- gsub("csv$","ped",filename)
  plinkfile <- gsub("csv$","GOOD",filename)
  path_plinkfile <- file.path(dirname,plinkfile)
  if(dirname==".") path_plinkfile <- gsub("./", "", path_plinkfile)

  if(file.exists(file.path(dirname, db))){
    newDB <- FALSE
  } else {
    newDB <- TRUE
  }


  fluidigm2PLINK(file=file,
                 out=out,
                 map=map,
                 plots=plots,
                 rearrange=rearrange,
                 fixNames=fixNames,
                 verbose=verbose,
                 verbosity=verbosity)

  out <- estimateErrors(file=file.path(dirname,pedfile),
                        db=db,
                        appendSamplesToDB=appendSamplesToDB,
                        keep.rep=keep.rep,
                        y.marker=y.marker,
                        x.marker=x.marker,
                        sp.marker=sp.marker,
                        neg_controls=neg_controls,
                        allele_error=allele_error,
                        marker_dropout=marker_dropout,
                        no_marker=no_marker,
                        male.y=male.y,
                        male.hetX=male.hetX,
                        female.y=female.y,
                        female.Xtot=female.Xtot,
                        female.hetXtot=female.hetXtot,
                        warning.noYtot=warning.noYtot,
                        warning.noHetXtot=warning.noHetXtot,
                        plots=plots,
                        sexing=sexing,
                        verbose=verbose,
                        verbosity=verbosity)

  if(newDB) db <- NA

  calculatePairwiseSimilarities(file=path_plinkfile,
                                db=db,
                                sexing=sexing,
                                verbose=verbose,
                                verbosity=verbosity)

  if(!newDB) path_plinkfile <- paste0(path_plinkfile,"_oDB")

  getPairwiseSimilarityLoci(file=path_plinkfile,
                            verbose=verbose,
                            verbosity=verbosity)

  similarityMatrix(file=path_plinkfile,
                   group=group,
                   verbose=verbose,
                   verbosity=verbosity)

  if(verbose>0)cat("\n ### All DONE!",date(),"\n","##############################################################\n")

  out
}
