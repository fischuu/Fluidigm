#' @title Run the Fluidigm Analysis Script Together
#'
#' @description
#' This function serves as a wrapper for the entire analysis pipeline. It takes a Fluidigm input file and performs several operations including conversion to PLINK format, error estimation, calculation of pairwise similarities, determination of pairwise similarity loci, and calculation of the similarity matrix.
#'
#' @param file A string specifying the path to the Fluidigm input file.
#' @param out A string specifying the output file name. If left empty, the original basename of the input file will be used.
#' @param db A string specifying the filepath to the database file. If not provided, the function will proceed with the existing data.
#' @param appendSamplesToDB A logical indicating whether new samples should be added to the database. Default is FALSE.
#' @param map A string specifying the filepath to the PlateDnoY.map file. If not provided, the function will use the map file with the same name as the ped file.
#' @param keep.rep A numeric value indicating the number of replicates to keep. Default is 1.
#' @param neg_controls A vector specifying the names of negative controls. Default is NA.
#' @param y.marker A vector specifying the Y markers for sexing. Default is NA.
#' @param x.marker A vector specifying the X markers for sexing. Default is NA.
#' @param sp.marker A vector specifying the markers used for species identification. Default is NA.
#' @param plots A logical indicating whether plots should be created. Default is TRUE.
#' @param allele_error A numeric value specifying the threshold for RERUN on Allele errors. Default is 5.
#' @param marker_dropout A numeric value specifying the threshold for RERUN on Marker dropout. Default is 15.
#' @param no_marker A numeric value specifying the number of markers. Default is 50.
#' @param male.y A numeric value specifying the threshold for sexing, male y-chromosome markers. Default is 3.
#' @param male.hetX A numeric value specifying the threshold for sexing, heterozygote x-chr markers. Default is 0.
#' @param female.y A numeric value specifying the threshold for sexing, female y-chromosome markers. Default is 0.
#' @param female.Xtot A numeric value specifying the threshold for sexing, total female x-chr markers. Default is 8.
#' @param female.hetXtot A numeric value specifying the threshold for sexing, heterozygote x-chr markers. Default is 3.
#' @param warning.noYtot A numeric value specifying the threshold for sexing, when should warning be triggered. Default is 2.
#' @param warning.noHetXtot A numeric value specifying the threshold for sexing, when should warning be triggered. Default is 3.
#' @param rearrange A logical indicating whether the ped/map output should be rearranged in order of provided map file. Default is TRUE.
#' @param group A string specifying the sample identifier for statistics. Default is NA.
#' @param fixNames A logical indicating whether whitespaces from sample names should be automatically removed. Default is TRUE.
#' @param sexing A logical indicating whether sexing should be performed. Default is FALSE.
#' @param verbose A logical or numerical value indicating whether the output should be verbose. Default is TRUE.
#' @param verbosity A numerical value indicating the level of verbosity. Set to a higher number for more details. Default is 1.
#'
#' @details
#' The function first checks the input parameters and sets default values if necessary. It then runs the following functions in order:
#' - fluidigm2PLINK: Converts the Fluidigm data to PLINK format.
#' - estimateErrors: Estimates errors in the PLINK ped files.
#' - calculatePairwiseSimilarities: Calculates pairwise similarities between genotypes.
#' - getPairwiseSimilarityLoci: Determines the loci of pairwise similarities.
#' - similarityMatrix: Calculates the similarity matrix.
#' The function prints a completion message when all operations are done.
#'
#' @examples
#' \dontrun{
#'   fluidigmAnalysisWrapper(file="path/to/your/file.csv", map="path/to/your/mapfile.map")
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

#  if(!verbose & verbosity > 0) verbosity <- 0
#  verbose <- verbosity
#  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)

  if(sexing) if(is.na(y.marker)) stop("You do not provide a vector with marker names to y.marker!")

  filename <- basename(file)
  dirname <- dirname(file)
  if(!is.na(out)){
    fluidigm_file <- paste0(out, ".csv")
  } else {
    fluidigm_file <- filename
  }
  pedfile <- gsub("csv$","ped",fluidigm_file)
  plinkfile <- gsub("csv$","GOOD",fluidigm_file)
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

  out <- estimateErrors(file=pedfile,
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

  if(verbose>0)message("\n ### All DONE!",date(),"\n","##############################################################\n")

  out
}
