#' @title Run plink to Calculate Pairwise Similarities
#'
#' @description
#' This function serves as a wrapper to the PLINK software, which is a free, open-source whole genome association analysis toolset.
#' It specifically uses PLINK to calculate pairwise similarities between genotypes.
#'
#' @param file A string representing the path to the filtered ped/map file pair (without ped/map file extension).
#' @param db A string representing the path to an existing genotype database. If not provided, the function will proceed with the existing data.
#' @param map A string representing the filepath to PlateDnoY.map file. If not provided, the function will use the map file with the same name as the ped file.
#' @param out A string representing the path to the output. If not provided, the output will be written to a file with the same name as the input file, appended with "_oDB".
#' @param sexing A logical value indicating whether the function should try to perform sexing. Default is FALSE.
#' @param verbose A logical value indicating whether the output should be verbose. Default is TRUE.
#' @param verbosity An integer representing the level of verbosity. Set to a higher number for more detailed output. Default is 1.
#'
#' @details
#' The function first checks the input parameters and sets default values if necessary. It then constructs and executes a PLINK command to merge
#' the genotype output with the existing genotype database if one is provided. Finally, it calculates pairwise similarities for all samples (and
#' database individuals) using another PLINK command. If the 'sexing' parameter is set to TRUE, the function will also attempt to determine the
#' sex of the individuals.
#'
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#'
#' @examples
#' \dontrun{
#'   calculatePairwiseSimilarities(file, db=NA, map=NA, out=NA, sexing=FALSE,
#'                                 verbose=TRUE, verbosity=1)
#' }
#'
#' @export


calculatePairwiseSimilarities <- function(file, db=NA, map=NA, out=NA, sexing=TRUE, verbose=TRUE, verbosity=1){

   ### Input checks
   ##############################################################################

   # Check if file has .ped or .map extension
   file_extension <- tools::file_ext(file)
   if (file_extension == "ped" || file_extension == "map") {
      # Remove the extension
      file <- tools::file_path_sans_ext(file)

      # Issue a warning if verbose is TRUE
      if (verbose) {
         warning("The file extension (.ped or .map) was removed from the input file name. The adjusted file name is: ", file)
      }
   }

   # Check if .ped file exists
   if (!file.exists(paste0(file, ".ped"))) {
      stop("The specified .ped file does not exist.")
   }

   # Check if .map file exists
   if (!file.exists(paste0(file, ".map"))) {
      stop("The specified .map file does not exist.")
   }

   # Check if db exists
   if (!is.na(db) && !file.exists(db)) {
     stop("The specified database does not exist.")
   }

   # Check if map exists
   if (!is.na(map) && !file.exists(map)) {
     stop("The specified map file does not exist.")
   }

 # Check if verbosity is a positive integer
    if (!is.numeric(verbosity) || verbosity < 0) {
      stop("The verbosity level should be a positive integer.")
    }

    if(!verbose & verbosity > 0) verbosity <- 0
    verbose <- verbosity
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
         if(verbose>0) message("No database.ped provided, we just continue with the existing data!\n")
     } else {
         if(verbose>0) message("A database.ped is provided, we combine it with the existing data!\n")
         plinkCommand <- paste0("plink --noweb --file ",file," --merge ",db," ",map," --recode --out ",out)

         if(verbose>1){
             message("Run the following PLINK command:\n", plinkCommand, "\n")
         }
         # Run Plink
         system(plinkCommand, intern=intern.param)
     }



   # Calculate pairwise similarities for all samples (and database individuals):
     if(sexing){
       plinkCommand <- paste0("plink --noweb --file ",out," --cluster --matrix --out ", out)
     } else {
       plinkCommand <- paste0("plink --noweb --file ",out," --cluster --matrix --allow-no-sex --out ", out)
     }
     if(verbose>1){
       message("Run the following PLINK command:\n", plinkCommand)
     }
     # Run Plink
     system(plinkCommand, intern.param)

    if(verbose>0) message("### Calculating pairwise similarities: DONE! ",date(),"\n","##############################################################\n")
}
