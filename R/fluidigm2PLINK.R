#' Turn Fluidigm Output Into PLINK Format
#'
#' This function loads a fluidigm output and turns it into a PLINK format
#'
#' @param file Path to the input file
#' @param out Out file name, keep empty to keep the original basename
#' @param map Filepath to PlateD_withY.map file
#' @param plot Logical, plot additional figures for conversion
#' @param rearrange Logical, rearrange the ped/map output in order of ymap
#' @param verbose Should the output be verbose, logical or numerical
#' @return A ped/map file pair and optional diagnostic plots
#' @export

fluidigm2PLINK <- function(file, out=NA, map=NA, plots=TRUE, rearrange=TRUE, verbose=TRUE){
  ### Input checks
  ##############################################################################
  if(!file.exists(map)) stop("The file 'map' does not exist, please check the path!")
  ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)

  ### Import the fluidigm data from csv file
  ##############################################################################

  # Extract the file and directory name
    dirname <- dirname(file)
    filename <- basename(file)
  # Read the fluidigm data into a data table
    tmp <- readLines(file)
    skip <- which(tmp=="SNP Converted Calls")
    dd1 <- data.table::fread(file, skip=skip +1)
  # Turn the data table into a data fram
    dd <- as.data.frame(unclass(dd1))
  # Add new level to V2
    levels(dd$V2)<-c(levels(dd$V2),"Chipblank")
  # Adjust the names for NTC
    dd$V2[dd$V2=="NTC"] <- "Chipblank"

  ### Create the MAP file
  ##############################################################################

  # snpids are in the first row, but first two columns are reserved for other input
    snpIDs <- dd[1,c(-1,-2)]
  # Create the map file
    ddmap <- data.frame(matrix(0, nrow = length(snpIDs), ncol = 4))
  # Add the SNP-IDs to the map
    ddmap$X2 <- t(snpIDs)
    newOrder <- 1:nrow(ddmap)

  if(!is.na(map)){
    # Check that markers are in correct order
      truemap <- read.table(map, sep = "\t", col.names= c("X1", "MAP", "X3", "X4"), colClasses = "character")
      comp <- as.factor(c(truemap$MAP)) == as.factor(c(ddmap$X2))

      if(sum(cumprod(comp)) < nrow(ddmap)){
        if(rearrange){
          if(verbose>1) print("Markers are not in the same order as in the provided map file, rearrange the output!")
          for(i in 1:nrow(truemap)){
            newOrder[i] <- which(truemap$MAP[i] == ddmap[,2])
          }
          ddmap <- ddmap[newOrder,]
        } else {
          stop("ERROR!!! Markers are not in correct order. Either change the order or set rearrange=TRUE to adjust the order of fluidigm file.")
        }
      } else {
        if(verbose>1) print("Markers are in the same order as map file")
      }
  } else {
     if(verbose>1) cat("No map file provided, create one based on marker IDs from csv file\n")
  }
  # Populate the new map file with information from the provided map file
    for(i in 1:nrow(ddmap)){
      if(truemap[i,1]!=0) ddmap[i,1] <- truemap[i,1]
      if(truemap[i,3]!=0) ddmap[i,3] <- truemap[i,3]
      if(truemap[i,4]!=0) ddmap[i,4] <- truemap[i,4]
    }
  # Export the .map file
    map.filename <- gsub(".csv", ".map", filename)
    write.table(ddmap, file = file.path(dirname, map.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ### Create the PED file
  ##############################################################################
  # Input checks
  #  Maybe rearrange the file
    if(rearrange){
      dd <- dd[,c(1,2,newOrder+2)]
    }
    ddped1 <- dd[-1,]
    number_samples <- nrow(ddped1)
    if(verbose>1) cat("Number of samples in data:",number_samples, "\n")
    number_markers <- ncol(ddped1)-2
    if(verbose>1) cat("Number of markers in data:",number_markers, "\n")

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

  # to replace missing genotypes add "0 0" instead of Invalid, NTC and No Call
    for(k in 7:ncol(ddped2)){
      levels(ddped2[,k]) <- c(levels(ddped2[,k]),"0 0")
    }
    ddped2[ddped2 == "Invalid"] <- "0 0"
    ddped2[ddped2 == "NTC"] <- "0 0"
    ddped2[ddped2 == "No Call"] <- "0 0"

  # export ped file
    ped.filename <- gsub(".csv", ".ped", filename)
    write.table(ddped2, file = file.path(dirname, ped.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ### Create summary statistics
  # call rate markers
    genomark <- table(ddped2[,7]!="0 0")["TRUE"]
    for(l in 8:ncol(ddped2)){
      b <- table(ddped2[,l]!="0 0")["TRUE"]
      genomark <- rbind(genomark,b)
    }
    genomark[is.na(genomark)]<- 0
    if(verbose>1){
      cat("Summary of call rate markers\n
             --------------------------------------\n")
      print(summary(genomark))
    }
    if(plots){
      fig1.filename <- sub(".csv", ".call_rate_markers.png", filename)
      png(file.path(dirname, fig1.filename), width=1200, height=1200)
         hist(genomark, breaks=number_markers, xlim=c(0,ceiling(number_markers/10)*10 ), main='', xlab="Call rate, markers")
      dev.off()
    }

  # genotyping success samples
    genosamp <- table(ddped2[1,c(7:ncol(ddped2))]!="0 0")["TRUE"]
    for(m in 2:nrow(ddped2)){
      c <- table(ddped2[m,c(7:ncol(ddped2))]!="0 0")["TRUE"]
      genosamp <- rbind(genosamp,c)
    }
    genosamp[is.na(genosamp)]<- 0
    if(verbose>1){
      cat("Summary of genotyping success\n
             --------------------------------------\n")
      print(summary(genosamp))
    }
    if(plots){
      fig2.filename <- sub(".csv", ".genotyping_success_samples.png", filename)
      png(file.path(dirname, fig2.filename), width=1200, height=1200)
         hist(genosamp, breaks=number_samples, xlim=c(0,ceiling(number_samples/10)*10), main='', xlab="Genotyping success, samples")
      dev.off()
    }

  # Some additional summary stats
    if(plots){
       fig3.filename <- sub(".csv", ".additional_summary_stats.png", filename)
       png(file.path(dirname, fig3.filename), width=1200, height=1200)
       par(mfrow=c(2,2))
       plot(1,1, type="n", axes=FALSE, ylab=" ", xlab=" ")
       text(1,1, paste("Input file name: ", filename),font=2, pos=1)
       text(1,1, paste("Number of markers in MAP file: ", nrow(ddmap)), pos=1, offset=2)
       text(1,1, paste("Number of samples in PED file: ", nrow(ddped2)), pos=1, offset=3)
       text(1,1, paste("Output map file name: ", map.filename), pos=1, offset=4)
       text(1,1, paste("Output ped file name: ", ped.filename), pos=1, offset=5)
       hist(genomark, breaks=number_markers, xlim=c(0,ceiling(number_markers/10)*10 ), main='', xlab="Call rate, markers")
       hist(genosamp, breaks=number_samples, xlim=c(0,ceiling(number_samples/10)*10 ), main='', xlab="Genotyping success, samples")
       dev.off()
    }

  if(verbose>0)cat("\n ### Conversion: DONE! ",date(),"\n","##############################################################\n")
}
