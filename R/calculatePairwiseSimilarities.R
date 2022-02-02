#' Run plink to Calculate Pairwise Similarities
#'
#' This function is a wrapper to plink to calculate pairwise similarities
#'
#' @param file Path to the filtered ped/map file pair (without ped/map file extension)
#' @param db Path to an existing genotype database
#' @param map Filepath to PlateDnoY.map file
#' @param out Path to the output, can be emptu
#' @param verbose Should the output be verbose, logical or numerical
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

calculatePairwiseSimilarities <- function(file, db=NA, map=NA, out=NA, verbose=TRUE){

   ### Input check
   ##############################################################################
    ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)
    ifelse(verbose>1, intern.param <- FALSE, intern.param <- TRUE)
    if(is.na(out)) out <- paste0(file,"_oDB")
    if(is.na(map)) map <- paste0(file,".map")

    if(is.na(db)) out <- file

    file <- gsub("^./", "", file)
    out <- gsub("^./", "", out)

   # Running plink
   # Merge the genotype output with the existing genotype database (example database included here "DataBase.ped"):
     if(is.na(db)){
         if(verbose>0) cat("No database.ped provided, we just continue with the existing data!\n")
     } else {
         if(verbose>0) cat("A database.ped is provided, we combine it with the existing data!\n")
         plinkCommand <- paste0("plink --noweb --file ",file," --merge ",db," ",map," --recode --out ",out)

         if(verbose>1){
             cat("Run the following PLINK command:\n", plinkCommand, "\n")
         }
         # Run Plink
         system(plinkCommand, intern=intern.param)
     }



   # Calculate pairwise similarities for all samples (and database individuals):
     plinkCommand <- paste0("plink --noweb --file ",out," --cluster --matrix --out ", out)
     if(verbose>1){
       cat("Run the following PLINK command:\n", plinkCommand)
     }
     # Run Plink
     system(plinkCommand, intern.param)

    if(verbose>0)cat("\n ### Calculating pairwise similarities: DONE! ",date(),"\n","##############################################################\n")
}
