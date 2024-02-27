#' @title Convert Fluidigm Output to PLINK Format
#'
#' @description
#' This function takes a Fluidigm output file and converts it into a PLINK format. The PLINK format is a widely used genetic
#' variation data format. This function is useful for researchers who want to analyze Fluidigm data using tools that accept
#' PLINK format.
#'
#' @param file A string specifying the path to the input file in CSV format.
#' @param map A string specifying the filepath to the map-file that should be used.
#' @param out A string specifying the output file name. If left empty, the original basename of the input file will be used.
#' @param outdir A string specifying the output folder. If left empty the original folder path of the input file will be used.
#' @param plots A logical indicating whether additional figures for conversion should be plotted. Default is TRUE.
#' @param rearrange A logical indicating whether the ped/map output should be rearranged in order of provided map file. Default is TRUE.
#' @param missing.geno A character string specifying how missing values should be coded. Default is "0 0".
#' @param fixNames A logical indicating whether whitespaces from sample names should be automatically removed. Default is TRUE.
#' @param overwrite A logical indicating wheter the original map file should be overwritten or not. Default FALSE
#' @param verbose A logical or numerical value indicating whether the output should be verbose. Default is TRUE.
#' @param verbosity A numerical value indicating the level of verbosity. Set to a higher number for more details. Default is 1.
#'
#' @details
#' The function first checks the input parameters and then imports the Fluidigm data from the CSV file. It creates a new MAP file
#' based on the information provided in the given map file. The function then creates a PED file and exports both files.
#' If requested, the function also generates plots for genotyping success and additional summary statistics.
#'
#' This function uses the PLINK software. For more information about PLINK, please refer to the official documentation.
#'
#' It might be so that marker names do not match between the csv and map file. This might happen, if special characters like ':' and '_'
#' are used. An easy way to harmonize the files is to apply sed with the '-i' option to change the file directly like this:
#'
#' sed -i 's/:/_/g' \<filename\>
#'
#' In above command all ':' will be substituted by '_' in the file <filename>.
#'
#' @references
#' PLINK: Whole genome data analysis toolset - Harvard University
#'
#' @examples
#' \dontrun{
#'   file_path_csv <- system.file("extdata", "example_data.csv", package = "Fluidigm")
#'   file_path_map <- system.file("extdata", "PlateD_withY.map", package = "Fluidigm")
#'   outdir <- tempdir()
#'
#'   fluidigm2PLINK(file=file_path, map=file_path_map, outdir=outdir)
#' }
#'
#' @return A list containing the ped/map file pair and optional diagnostic plots.
#' @export


fluidigm2PLINK <- function(file=NA, map=NA, out=NA, outdir=NA, plots=TRUE, rearrange=TRUE, missing.geno="0 0",
                           fixNames=TRUE, overwrite=FALSE, verbose=TRUE, verbosity=1){

  ### Input checks
  ##############################################################################
  if(!verbose & verbosity > 0) verbosity <- 0
  verbose <- verbosity
  if(is.na(file)) stop("Please provide a csv file!")
  if(is.na(map)) stop("Please provide a map file!")
  if(!file.exists(file)) stop("The file 'file' does not exist, please check the path!")
  if(!file.exists(map)) stop("The file 'map' does not exist, please check the path!")
  if(!overwrite & is.na(out)){
    if(is.na(outdir)){
      stop("Neither a new name for output nor an alternative output folder provided and old output is set to 'no overwrite'. Please change either!")
    } else {
      if(verbose>1) message("New output folder provided, hence the same name can be used with the overwrite option set to TRUE.")
    }
  }


  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)


  if(verbose){
    if(length(grep(" ", file))>0) warning("There are whitespaces in your input file, this will most likely crash your plink system call!\n
                                          My suggestion, change the whitespaces with underscores ('_')  and rerun.")
  }

  ### Import the fluidigm data from csv file
  ##############################################################################

  # Extract the file and directory name
    dirname <- dirname(file)
    if(!is.na(outdir)) dirname <- outdir
    filename.in <- basename(file)
    if(!is.na(out)) {
      filename.out <- out
    } else {
      filename.out <- filename.in
    }
  # Read the fluidigm data into a data table
    tmp <- readLines(file)
    skip <- which(tmp=="SNP Converted Calls")
    dd1 <- data.table::fread(file, skip=skip +1, header = TRUE)
  # Turn the data table into a data frame
    dd <- as.data.frame(unclass(dd1))
    dd$V2 <- as.factor(dd$V2)
  # Add new level to V2
    levels(dd$V2) <- c(levels(dd$V2),"Chipblank")
  # Adjust the names for NTC
    dd$V2[dd$V2=="NTC"] <- "Chipblank"

  ### Create the MAP file
  ##############################################################################

  # snpids are in the first row, but first two columns are reserved for other input
    snpIDs <- colnames(dd)[-c(1:2)]
  # Create the map file
    ddmap <- data.frame(matrix(0, nrow = length(snpIDs), ncol = 4))
  # Add the SNP-IDs to the map
    ddmap$X2 <- snpIDs
    newOrder <- 1:nrow(ddmap)

  if(!is.na(map)){
    # Check that markers are in correct order
      truemap <- read.table(map, sep = "\t", col.names= c("X1", "MAP", "X3", "X4"), colClasses = "character")
      #comp <- as.factor(c(truemap$MAP)) == as.factor(c(ddmap$X2))
      comp <- truemap$MAP == ddmap$X2

      if(sum(comp) < nrow(truemap)){
        if(sum(is.element(truemap$MAP,ddmap$X2))<nrow(truemap)){
          if(verbose){
            message("Mismatching entries between MAP entries:\n-----------------------------------------------------------------\n")
            printThis <- cbind(ddmap$X2[!comp], truemap$MAP[!comp])
            colnames(printThis) <- c("file-input", "map-input")
            print(printThis)
            stop("ERROR: You need to fix first the marker names before you can proceed!")
          }
        }
        if(rearrange){
          if(verbose>1) message("Markers are not in the same order as in the provided map file, rearrange the output!")
          for(i in 1:nrow(truemap)){
            newOrderEntry <- which(truemap$MAP[i] == ddmap[,2])
            if(length(newOrderEntry)==0){
              if(verbose) stop("Entry ", truemap$MAP[i]," from map-file cannot be found in ", file, "\n")
            }
            newOrder[i] <- newOrderEntry
          }
          ddmap <- ddmap[newOrder,]
        } else {
          stop("ERROR!!! Markers are not in correct order. Either change the order or set rearrange=TRUE to adjust the order of fluidigm file.")
        }
      } else {
        if(verbose>1) message("Markers are in the same order as map file")
      }
  } else {
     if(verbose>1) message("No map file provided, create one based on marker IDs from csv file\n")
  }

  # Populate the new map file with information from the provided map file
    for(i in 1:nrow(ddmap)){
      if(truemap[i,1]!=0) ddmap[i,1] <- truemap[i,1]
      if(truemap[i,3]!=0) ddmap[i,3] <- truemap[i,3]
      if(truemap[i,4]!=0) ddmap[i,4] <- truemap[i,4]
    }
  # Export the .map file
    #map.filename <- gsub(".csv", ".map", filename)
    map.filename <- paste0(filename.out, ".map")
    if(verbose>1){
      message("Export", nrow(ddmap), "markers into the new created map file:",map.filename,"\n")
    }
    write.table(ddmap, file = file.path(dirname, map.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ### Create the PED file
  ##############################################################################
  # Input checks
  #  Maybe rearrange the file
    if(rearrange){
      dd <- dd[,c(1,2,newOrder+2)]
    }
    ddped1 <- dd
    number_samples <- nrow(ddped1)
    if(verbose>1) message("Number of samples in data:",number_samples, "\n")
    number_markers <- ncol(ddped1)-2
    if(verbose>1) message("Number of markers in data:",number_markers, "\n")

  # Add columns for sex-information
    ddped <- ddped1
    ddped$sex <- "NA"
    ddped$sexnum <- 0

  # Add columns for parental id and phenotype
    ddped$patID <- 0
    ddped$matID <- 0
    ddped$pheno <- 0

  # Rearrange to get right order -- use indexing
    p <- which(colnames(ddped)=="patID")
    m <- which(colnames(ddped)=="matID")
    s <- which(colnames(ddped)=="sexnum")
    ph <- which(colnames(ddped)=="pheno")
    l <- ncol(ddped)-5
    ddped2<-ddped[,c(1,2,p,m,s,ph,3:l)]
  # If order needed to be rearranned

  # TODO: This is not effective and for now rather naive
  # reformat genotype (and missing genotype) data
  # to replace ":" with " " in genotypes using double loop (over both columns and rows)
    for(i in 7:ncol(ddped2)){
      for(j in 1:nrow(ddped2)){
        levels(ddped2[,i]) <- unique(c(levels(ddped2[,i]), gsub(":"," ",ddped2[j,i])))
        ddped2[j,i] <- gsub(":"," ",ddped2[j,i])
      }
    }

  # to replace missing genotypes add "missing.geno" instead of Invalid, NTC and No Call
    for(k in 7:ncol(ddped2)){
      levels(ddped2[,k]) <- c(levels(ddped2[,k]),missing.geno)
    }
    ddped2[ddped2 == "Invalid"] <- missing.geno
    ddped2[ddped2 == "NTC"] <- missing.geno
    ddped2[ddped2 == "No Call"] <- missing.geno

  # Fix the sample names
    if(fixNames){
      ddped2$V2 <- sub(" ", "_", ddped2$V2)
    } else {
      if(verbose){
        if(length(grep(" ", ddped2$V2))>0) warning("There are whitespaces in the sample names, most likely the downstream analysis will fail!\n
                                                   Fix it manually and rerun this function or set 'fixNames=TRUE'")
      }
    }

  # export ped file
#    ped.filename <- gsub(".csv", ".ped", filename)
    ped.filename <- paste0(filename.out, ".ped")
    write.table(ddped2, file = file.path(dirname, ped.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ### Create summary statistics
  # call rate markers
    genomark <- table(ddped2[,7]!=missing.geno)["TRUE"]
    for(l in 8:ncol(ddped2)){
      b <- table(ddped2[,l]!=missing.geno)["TRUE"]
      genomark <- rbind(genomark,b)
    }
    genomark[is.na(genomark)]<- 0
    if(verbose>1){
      message("Summary of call rate markers\n--------------------------------------\n")
      print(summary(genomark))
    }
    if(plots){
      #fig1.filename <- sub(".csv", ".call_rate_markers.png", filename)
      fig1.filename <- paste0(filename.out, ".call_rate_markers.png")
      png(file.path(dirname, fig1.filename), width=1200, height=1200)
         hist(genomark, breaks=number_markers, xlim=c(0,ceiling(number_markers/10)*10 ), main='', xlab="Call rate, markers")
      dev.off()
    }

  # genotyping success samples
    genosamp <- table(ddped2[1,c(7:ncol(ddped2))]!=missing.geno)["TRUE"]
    for(m in 2:nrow(ddped2)){
      c <- table(ddped2[m,c(7:ncol(ddped2))]!=missing.geno)["TRUE"]
      genosamp <- rbind(genosamp,c)
    }
    genosamp[is.na(genosamp)]<- 0
    if(verbose>1){
      message("Summary of genotyping success\n--------------------------------------\n")
      print(summary(genosamp))
    }
    if(plots){
      #fig2.filename <- sub(".csv", ".genotyping_success_samples.png", filename)
      fig2.filename <- paste0(filename.out, ".genotyping_success_samples.png")
      png(file.path(dirname, fig2.filename), width=1200, height=1200)
         hist(genosamp, breaks=number_samples, xlim=c(0,ceiling(number_samples/10)*10), main='', xlab="Genotyping success, samples")
      dev.off()
    }

  # Some additional summary stats
    if(plots){
       #fig3.filename <- sub(".csv", ".additional_summary_stats.png", filename)
       fig3.filename <- paste0(filename.out, ".additional_summary_stats.png")
       png(file.path(dirname, fig3.filename), width=1200, height=1200)
       par(mfrow=c(2,2))
       plot(1,1, type="n", axes=FALSE, ylab=" ", xlab=" ")
       text(1,1, paste("Input file name: ", filename.in),font=2, pos=1)
       text(1,1, paste("Number of markers in MAP file: ", nrow(ddmap)), pos=1, offset=2)
       text(1,1, paste("Number of samples in PED file: ", nrow(ddped2)), pos=1, offset=3)
       text(1,1, paste("Output map file name: ", map.filename), pos=1, offset=4)
       text(1,1, paste("Output ped file name: ", ped.filename), pos=1, offset=5)
       hist(genomark, breaks=number_markers, xlim=c(0,ceiling(number_markers/10)*10 ), main='', xlab="Call rate, markers")
       hist(genosamp, breaks=number_samples, xlim=c(0,ceiling(number_samples/10)*10 ), main='', xlab="Genotyping success, samples")
       dev.off()
    }

  if(verbose>0)message("\n### Conversion: DONE! ",date(),"\n","##############################################################\n")
}
