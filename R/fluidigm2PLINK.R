#' Turn fluidigm into PLINK format
#'
#' This function loads a fluidigm output and turns it into a PLINK format
#'
#' @param file Path to the input file
#' @param out Out file name, keep empty to keep the original basename
#' @param ymap Filepath to PlateD_withY.map file
#' @param verbose Should the output be verbose, logical
#' @return Something
#' @export

fluidigm2PLINK <- function(file, out=NA, ymap="PlateD_withY.map", verbose=TRUE, plots=TRUE){
  ### Input checks
  ##############################################################################
  if(!file.exists(ymap)) stop("The file 'ymap' does not exist, please check the path!")

  ### Import the fluidigm data
  ##############################################################################

  # Extract the file and directory name
    dirname <- dirname(file)
    filename <- basename(file)
  # Read the fluidigm into a data table
    dd1 <- data.table::fread(file, skip="SNP Converted Calls")
  # Turn the data table into a data fram
    dd <- as.data.frame(unclass(dd1[-1,]))
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
  # Check that markers are in correct order
    truemap <- read.table(ymap, sep = "\t", col.names= c("X1", "MAP", "X3", "X4"), colClasses = "character")
    comp <- as.factor(c(truemap$MAP)) == as.factor(c(ddmap$X2))
    if(sum(cumprod(comp)) < nrow(ddmap)){
      stop("WARNING! Markers are not in correct order")
    } else {
      if(verbose) print("OK markers are in correct order")
    }
  # Export the .map file
    map.filename <- gsub(".csv", ".map", filename)
    write.table(ddmap, file = file.path(dirname, map.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  ### Create the PED file
  ##############################################################################
  # Input checks
    ddped1 <- dd[-1,]
    if(verbose) cat("Number of samples in data:",nrow(ddped1), "TODO: CHECK THIS!!!\n")
    if(verbose) cat("Number of markers in data:",ncol(ddped1), "TODO: CHECK THIS!!!\n")

  # Add columns for sex-information
    ddped <- ddped1
    ddped$sex <- "NA"
    ddped$sexnum <- 0

  # Add columns for parental id and phenotype
    ddped$patID <- 0
    ddped$matID <- 0
    ddped$pheno <- 0

  # Rearange to get right order -- use indexing
    p <- ncol(ddped)-2
    m <- ncol(ddped)-1
    s <- ncol(ddped)-3
    ph <- ncol(ddped)
    l <- ncol(ddped)-5
    ddped2<-ddped[,c(1,2,p,m,s,ph,3:l)]

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
    genomark<-table(ddped2[,7]!="0 0")[2]
    for(l in 8:ncol(ddped2)){
      b<-table(ddped2[,l]!="0 0")[2]
      genomark<-rbind(genomark,b)
    }
    genomark[is.na(genomark)]<- 0
    if(verbose){
      cat("Summary of call rate markers\n
             --------------------------------------\n")
      print(summary(genomark))
    }
    if(plots){
      fig1.filename <- sub(".csv", ".call_rate_markers.png", filename)
      png(file.path(dirname, fig1.filename), width=1200, height=1200)
         hist(genomark, breaks=96, xlim=c(0,96), main='', xlab="Call rate, markers")
      dev.off()
    }

  # genotyping success samples
    genosamp <- table(ddped2[1,c(7:ncol(ddped2))]!="0 0")[2]
    for(m in 2:nrow(ddped2)){
      c<-table(ddped2[m,c(7:ncol(ddped2))]!="0 0")[2]
      genosamp <- rbind(genosamp,c)
    }
    genosamp[is.na(genosamp)]<- 0
    if(verbose){
      cat("Summary of genotyping success\n
             --------------------------------------\n")
      print(summary(genosamp))
    }
    if(plots){
      fig2.filename <- sub(".csv", ".genotyping_success_samples.png", filename)
      png(file.path(dirname, fig2.filename), width=1200, height=1200)
         hist(genosamp, breaks=96, xlim=c(0,96), main='', xlab="Genotyping success, samples")
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
       hist(genomark, breaks=96, xlim=c(0,96), main='', xlab="Call rate, markers")
       hist(genosamp, breaks=96, xlim=c(0,96), main='', xlab="Genotyping success, samples")
       dev.off()
    }
  if(verbose)print("Conversion: DONE!\n")
}
