#' @title Process the PLINK ped files
#'
#' @description
#' This function processes the PLINK ped files
#'
#' @param file Path to the ped input file
#' @param db Name of the used database
#' @param appendSamplesToDB Logical, shoud new samples be added to database
#' @param keep.rep numeric, keep only this n-fold replicates, default n=2
#' @param y.marker Y markers for sexing
#' @param x.marker X markers for sexing
#' @param sp.marker Markers used for species-identification
#' @param plots logical, shall plots be created?
#' @param neg_controls Names of negative controls
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
#' @param sexing Logical, if sexing should be performed
#' @param verbose Should the output be verbose, logical or numerical
#' @param verbosity Level of verbosity, set to higher number for more details
#'
#' #' @details
#' Additional details...
#'
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#'
#' @examples
#' \dontrun{
#'   estimateErrors()
#' }
#' @export

estimateErrors <- function(file, db=NA, appendSamplesToDB=FALSE, keep.rep=1,
                           y.marker=NA, x.marker=NA, sp.marker=NA, plots=TRUE, neg_controls=NA,
                           allele_error=5, marker_dropout=15, no_marker=50,
                           male.y=3, male.hetX=0, female.y=0, female.Xtot=8, female.hetXtot=3,
                           warning.noYtot=2, warning.noHetXtot=3, sexing=TRUE, verbose=TRUE, verbosity=1){

  # Input checks
    if(!verbose & verbosity > 0) verbosity <- 0
    verbose <- verbosity
    ifelse(as.numeric(verbose)>0, verbose <- as.numeric(verbose) , verbose <- 0)
    if(sexing) if(is.na(y.marker)) stop("You do not provide a vector with marker names to y.marker!")

  # Welcome screen
  if(sexing){
    if(verbose>0){
      cat("Sex determination settings:\n")
      cat("-------------------------------\n")
      cat("MALE      : Call as MALE, if number of Y total (noYtot) >=", male.y,"(male.y - option) AND noHetXtot <=", male.hetX,"(male.hetX - option)\n")
      cat("FEMALE    : Call as FEMALE, if number of Y total (noYtot) equals == ",female.y,"(female.y - option) AND number of total X (noXtot)is >= ",female.Xtot,"(female.Xtot - option)\n")
      cat("FEMALE    : Call as FEMALE, if number of Y total (noYtot) equals == ",female.y,"(female.y - option) AND noHetXtot >= ",female.hetXtot, "(female.hetXtot - option)\n")
      cat("WARNING   : Call as WARNING, if if number of Y total (noYtot) >= ",warning.noYtot,"(warning.noYtot - option) and noHetXtot >= ",warning.noHetXtot,"(warning.noHetXtot - option\n")
      cat("UNCERTAIN : Call as Uncertain, if marker dropout > 25% \n")
      cat("Unhandled : All other cases will be marked as unhandled (this case should not happen...)\n")
    }
  }


  # Import sample genotypes:
    dirname <- dirname(file)
    filename <- basename(file)
    mapfile <- gsub("ped$","map",filename)
    basename <- sub(".CMR.ped","",filename)
    ped1 <- read.table(file, stringsAsFactors=TRUE) # make sure all genotype columns are treated as factors
    map1 <- read.table(file.path(dirname, mapfile), stringsAsFactors=TRUE) # make sure all genotype columns are treated as factors

  # Check to remove replicates (removes neg controls):
    n <- table(ped1$V2)
    replicates <- FALSE
    if(sum(n>1)>0) replicates <- TRUE
    remove_those <- names(n)[n!=keep.rep]
    if(!is.na(neg_controls)) remove_those <- c(remove_those, neg_controls)
    if(length(which(is.element(ped1$V2, remove_those)))>0){
      ped <- ped1[-which(is.element(ped1$V2, remove_those)),]
      ped$V2 <- factor(ped$V2)
      if(verbose>1){
        cat("Remove samples (based on too low/high repitition \n------------------------------------------------\n")
        for(i in 1:length(remove_those)){
          cat(remove_those[i],"\n")
        }
      }
    } else {
      if(verbose>0) warning("No samples were removed, please check if this is correct or if e.g. wrong sample names are provided for negative controls?!\n")
      ped <- ped1
    }


 # Remove y-chromosome based markers
  if(sum(!is.na(y.marker))>0 ){
    marker_pos <- NULL
    if(sum(is.numeric(y.marker))>0 | sum(y.marker=="Y")>0 ){
      marker_pos <- which(is.element(map1$V1, y.marker))
    } else {
      marker_pos <- which(is.element(map1$V2, y.marker))
    }
    if(length(marker_pos)==0) warning("No y-markers found. Please check your option y.marker to be either numeric and marker names. The map file needs to have then either information to be present.\n")
    if(verbose>1){
      cat("Remove markers\n---------------------------\n")
      for(i in 1:length(marker_pos)){
        cat(as.character(map1$V2[marker_pos[i]]),"\n")
      }
    }
    remove_y <- c(marker_pos*2+5, marker_pos*2+6)
    map1_wo_y <- map1[-marker_pos,]
    ped_wo_y <- ped[,-remove_y]
  } else {
  # If nothing should be removed, just copy the original files to process them further
    map1_wo_y <- map1
    ped_wo_y <- ped
  }

  ### Start main analyses for summary output
  ##############################################################################

  # Create new data frame
    nam <- as.character(ped_wo_y[,2])
    gensim <- data.frame(matrix(NA, nrow=ncol(ped_wo_y[,-c(1:6)]), ncol=length(unique(nam))))
    names(gensim) <- as.character(unique(nam))

  # Check if first and second run of sample have the same genotype
  # loop over all individuals (SAMPLE NAME ends up in header row)
    for(i in 1:length(unique(nam))){
      # get all genotypes for individual i
        ind <- unique(nam)[i]
        indgen <- ped_wo_y[ped_wo_y$V2==ind,]
        gen <- indgen[,-c(1:6)]

      # replace missing data with "NA"
        gen[gen==0] <- NA

      # check if both replicates are identical:
      if(nrow(gen)>1){
        replicates <- TRUE
        if(verbose>1) if(i==1) warning("Replicates found, object gensim will contain the agreement between replicate 1 and 2. ALL OTHERS WILL HERE BE IGNORED AT THE MOMENT!!!\n")
        for(j in 1:ncol(gen)){
          gensim[j,i]<- gen[1,j]==gen[2,j]
        }
      } else {
        if(verbose>1) if(i==1) warning("No replicates found, object gensim will contain only TRUE and NA\n")
        for(j in 1:ncol(gen)){
          gensim[j,i] <- gen[1,j]==gen[1,j]
        }
      }
    }

  # Create new data frame (summs) with summary statistics
    summs <- data.frame(matrix(NA, nrow=ncol(gensim), ncol=6))
    names(summs) <- c("Ind","Ntrue","Nfalse","Nna","Allele_error","Marker_dropout")
    summs$Ind <- unique(nam)

  # Now fill in error rates in data.frame (using a loop this time)
  # 5) Allele error: no mismatch/number alleles (alleles) *100 %
  # 6) Marker dropout: Number NA's/2 / number markers (alleles/2) *100 %
    for(k in 1:nrow(summs)){
      summs[k,2] <- sum(gensim[,k], na.rm=TRUE)
      summs[k,3] <- sum(gensim[,k]=="FALSE", na.rm=TRUE)
      summs[k,4] <- sum(is.na(gensim[,k]))
      summs[k,5] <- summs[k,3]/(summs[k,2]+summs[k,3])*100
      summs[k,6] <- summs[k,4]/(summs[k,2]+summs[k,3]+summs[k,4])*100
    }

    if(!replicates) if(verbose>1) warning("There are no replicates, not all values in summs are meaningful, also the marker classification is not based on the allele_error statistic!!!\n")

  # Add columns: "No_markers_repl1" and "No_markers_repl2"
    summs$No_markers_repl1 <- NA
    summs$No_markers_repl2 <- NA

  # Check if first and second run of sample have no genotype
  # loop over all individuals
    for(i in 1:length(unique(nam))){
      # get all genotypes for individual i
        ind <- unique(nam)[i]
        indgen <- ped_wo_y[ped_wo_y$V2==ind,]
        gen <- indgen[,-c(1:6)]

      # count number of cases with "0" as genotype (devide by 2 for no loci)
        nogeno1 <- sum(gen[1,]==0, na.rm=TRUE)/2
        nogeno2 <- sum(gen[2,]==0, na.rm=TRUE)/2

      # put these in the right place in table "summs"
        summs$No_markers_repl1[summs$Ind==ind] <- nogeno1
        summs$No_markers_repl2[summs$Ind==ind] <- nogeno2
      }

  ### Add category
  # "BAD" if marker dropout > 85% (needs to be tweaked to rerun also samples where one was good and one bad...)
  # "BAD" if both replicates have No_markers_repl1 AND No_markers_repl2 > 79 (85% of 93)
  # "RERUN" if marker dropout 85% - 15%, or error rate > 5%
  # "GOOD" else
    summs$categ <- as.factor("GOOD")
    levels(summs$categ) <- c("GOOD","RERUN","BAD")

    if(verbose > 1){
      cat("Marker classification is based on those thresholds:\n",
          "Allele errors between replicates (if present):", allele_error,"\n",
          "Marker dropout :", marker_dropout,"\n",
          "Number of missing marker per replicate:",no_marker , "\n")
    }
    summs$categ[summs$Allele_error > allele_error] <- "RERUN"
    summs$categ[summs$Marker_dropout > marker_dropout] <- "RERUN"
    if(replicates){
      summs$categ[summs$No_markers_repl1 > no_marker & summs$No_markers_repl2 > no_marker] <- "BAD"
    } else {
      summs$categ[summs$No_markers_repl1 > no_marker] <- "BAD"
    }


  ### Summarise error rates in figure
    culr <- ifelse(summs$categ=="GOOD", "blue", "red")
    culr[summs$categ=="RERUN"] <- "purple"

    if(plots){
      fig1.filename <- sub(".ped", ".allele_error_marker_dropout.png", filename)
      png(file.path(dirname, fig1.filename), width=1200, height=1200)
        par(mfrow=c(2,2), cex=1.2)
        hist(summs$Allele_error, breaks=50, xlim=c(0,100), main="", xlab="Allele Error (%)")
        lines(c(5,5), c(0,100), lty=2)
        hist(summs$Marker_dropout, breaks=50, xlim=c(0,100), main="", xlab="Marker Dropout (%)")
        lines(c(15,15), c(0,100), lty=2)
        lines(c(85,85), c(0,100), lty=5)
        plot(summs$Allele_error ~ summs$Marker_dropout,
             xlab="Marker dropout rate (%)", ylab="Allele error rate (%)", pch=21, bg=culr)
        lines(c(-5,15), c(5,5), lty=2)
        lines(c(15,15), c(-5,5), lty=2)
        lines(c(85,85), c(-5,100), lty=5)
        if(replicates){
          corA <- cor.test(summs$Marker_dropout, summs$Allele_error)
        } else {
          corA <- list()
          corA$estimate <- NA
          corA$p.value <- NA
        }

        plot(2,2, type="n", axes=FALSE, ylab=" ", xlab=" ")
        text(2,2, paste("Input file name: ", filename),font=1, pos=1)
        text(2,2, paste("R-value (for replicates): ", corA$estimate), pos=1, offset=2)
        text(2,2, paste("p-value (for replicates): ", corA$p.value), pos=1, offset=3)
      dev.off()
    }

  ### Sex assignment
  if(sexing){
    summs$noHetX1 <- NA
    summs$noHetX2 <- NA
    summs$noHetXtot <- NA
    summs$noX1 <- NA
    summs$noX2 <- NA
    summs$noXtot <- NA
    summs$noY1 <- NA
    summs$noY2 <- NA
    summs$noYtot <- NA

  # loop over all individuals (i)
    for(i in 1:length(unique(nam)))    {
    # get all genotypes for individual i
      ind <- unique(nam)[i]
      indgen <- ped[ped$V2==ind,]
      if(verbose>1 && i==1){
        cat("\nConsider the following columns in ped-file:\n")
        for(indY in 1:(length(remove_y)/2)){
          writeThis <- indY*2
          cat(remove_y[c(writeThis-1, writeThis)], "(", as.character(map1$V2[floor((remove_y[writeThis-1])/2)[1]-2]),")\n")
        }

      }
      genY <- indgen[,remove_y[1:(length(remove_y)/2)]]

     # count number of cases with something else than "0" as genotype
      if(!is.null(nrow(genY))) {
        Ygeno1 <- sum(genY[1,]!=0, na.rm=TRUE)
        Ygeno2 <- sum(genY[2,]!=0, na.rm=TRUE)
      } else {
        Ygeno1 <- sum(genY!=0, na.rm=TRUE)
        Ygeno2 <- NA
      }

    # put these in the right place in table "summs"
      summs$noY1[summs$Ind==ind] <- Ygeno1
      summs$noY2[summs$Ind==ind] <- Ygeno2
      summs$noYtot[summs$Ind==ind] <- sum(Ygeno1, Ygeno2, na.rm=TRUE)
    }

  # Check genotypes of X markers
  # loop over all individuals (i)
    if(sum(is.na(x.marker))>0) stop("No x.marker provided, please do this.")
    x_marker_pos <- NULL
    if(sum(is.numeric(x.marker))>0 | sum(x.marker=="X")>0){
       x_marker_pos <- which(is.element(map1$V1, x.marker))
    } else {
       x_marker_pos <- which(is.element(map1$V2, x.marker))
    }
    remove_x <- sort(c(x_marker_pos*2+5, x_marker_pos*2+6))

    if(verbose>1){
      cat("Base x-marker analysis on the following markers:\n")
      for(i in 1:length(x_marker_pos)){
        cat(as.character(map1$V2)[x_marker_pos[i]],"\n")
      }
      cat("\nConsider the following columns in ped-file:\n")
      cat(remove_x,"\n")
    }

    for(i in 1:length(unique(nam))){
    # get all genotypes for individual i
      ind <- unique(nam)[i]
      indgen <- ped[ped$V2==ind,]
      genX <- indgen[,remove_x]
      genX[] <- lapply(genX, factor) 				# 3 new lines - change all variables in genX to factor
      levs <- unique(unlist(lapply(genX,levels)))		# Get vector of all levels that appear in the data.frame
      genX <- data.frame(lapply(genX,factor,levels=levs))	# Set these as the levels for each column

    # count number of of heterozygote genotypes
      Xgeno1 <- 0
      Xgeno1[genX[1,1] != genX[1,2]]  <- Xgeno1+1
      Xgeno1[genX[1,3] != genX[1,4]]  <- Xgeno1+1
      Xgeno1[genX[1,5] != genX[1,6]]  <- Xgeno1+1
      Xgeno1[genX[1,7] != genX[1,8]]  <- Xgeno1+1
      Xgeno1[genX[1,9] != genX[1,10]]  <- Xgeno1+1
      Xgeno1[genX[1,11] != genX[1,12]]  <- Xgeno1+1

      if(replicates){
        Xgeno2 <- 0
      } else {
        Xgeno2 <- NA
      }
      Xgeno2[genX[2,1] != genX[2,2]]  <- Xgeno2+1
      Xgeno2[genX[2,3] != genX[2,4]]  <- Xgeno2+1
      Xgeno2[genX[2,5] != genX[2,6]]  <- Xgeno2+1
      Xgeno2[genX[2,7] != genX[2,8]]  <- Xgeno2+1
      Xgeno2[genX[2,9] != genX[2,10]]  <- Xgeno2+1
      Xgeno2[genX[2,11] != genX[2,12]]  <- Xgeno2+1

    # put these in the right place in table "summs"
      summs$noHetX1[summs$Ind==ind] <- Xgeno1
      summs$noHetX2[summs$Ind==ind] <- Xgeno2
      summs$noHetXtot[summs$Ind==ind] <- sum(Xgeno1,Xgeno2, na.rm=TRUE)

    # count total number of working X-markers
      Xtgeno1 <- sum(genX[1,]!=0, na.rm=TRUE)/2
      Xtgeno2 <- sum(genX[2,]!=0, na.rm=TRUE)/2
      if(!replicates) Xtgeno2 <- NA

    # put these in the right place in table "summs"
      summs$noX1[summs$Ind==ind] <- Xtgeno1
      summs$noX2[summs$Ind==ind] <- Xtgeno2
      summs$noXtot[summs$Ind==ind] <- sum(Xtgeno1,Xtgeno2, na.rm=TRUE)
    }

    # Add column for assigned sex, and fill in using criteria below
      summs$sex <- as.factor("Unhandled")
      levels(summs$sex) <- c("Uncertain", "Unhandled","Female","Male","WARNING")
    # Call as Male if number Ytot > 2 AND number HetXtot < 1
      summs$sex[summs$noYtot >= male.y & summs$noHetXtot <= male.hetX] <- "Male"
    # Call as Female if number Ytot == 0 AND noXtot > 8
      summs$sex[summs$noYtot == female.y & summs$noXtot >= female.Xtot] <- "Female"
    # Also call as Female if number Ytot < 1 AND HetXtot > 2
      summs$sex[summs$noYtot == female.y & summs$noHetXtot >= female.hetXtot] <- "Female"
    # Call as "WARNING" if both male and female (number Ytot > 2 AND HetX > 2)
      summs$sex[summs$noYtot >= warning.noYtot & summs$noHetXtot >= warning.noHetXtot] <- "WARNING"
    # lastly call all sex as "Uncertain" if Marker dropout > 25%
      summs$sex[summs$Marker_dropout > 25] <- "Uncertain"


  # Process the species marker set sp.marker
      # Remove y-chromosome based markers
      if(sum(!is.na(sp.marker))>0 ){
        ### Sex assignment
        summs$noSPHetX1 <- NA
        summs$noSPHetX2 <- NA
        summs$noSPHetXtot <- NA
        summs$noSPX1 <- NA
        summs$noSPX2 <- NA
        summs$noSPXtot <- NA

        sp_marker_pos <- NULL
        if(sum(is.numeric(sp.marker))>0 | sum(sp.marker=="X")>0){
          sp_marker_pos <- which(is.element(map1$V1, sp.marker))
        } else {
          sp_marker_pos <- which(is.element(map1$V2, sp.marker))
        }
        remove_sp <- sort(c(sp_marker_pos*2+5, sp_marker_pos*2+6))

        if(verbose>1){
          cat("Base species-marker analysis on the following markers:\n")
          for(i in 1:length(sp_marker_pos)){
            cat(as.character(map1$V2)[sp_marker_pos[i]],"\n")
          }
          cat("\nConsider the following columns in ped-file:\n")
          cat(remove_sp,"\n")
        }

        for(i in 1:length(unique(nam))){
          # get all genotypes for individual i
          ind <- unique(nam)[i]
          indgen <- ped[ped$V2==ind,]
          genSP <- indgen[,remove_sp]
          genSP[] <- lapply(genSP, factor) 				# 3 new lines - change all variables in genX to factor
          levs <- unique(unlist(lapply(genSP,levels)))		# Get vector of all levels that appear in the data.frame
          genSP <- data.frame(lapply(genSP,factor,levels=levs))	# Set these as the levels for each column

          # count number of of heterozygote genotypes
          SPgeno1 <- 0
          SPgeno1[genSP[1,1] != genSP[1,2]]  <- SPgeno1+1
          SPgeno1[genSP[1,3] != genSP[1,4]]  <- SPgeno1+1
          SPgeno1[genSP[1,5] != genSP[1,6]]  <- SPgeno1+1
          SPgeno1[genSP[1,7] != genSP[1,8]]  <- SPgeno1+1
          SPgeno1[genSP[1,9] != genSP[1,10]]  <- SPgeno1+1
          SPgeno1[genSP[1,11] != genSP[1,12]]  <- SPgeno1+1

          if(replicates){
            SPgeno2 <- 0
          } else {
            SPgeno2 <- NA
          }
          SPgeno2[genSP[2,1] != genSP[2,2]]  <- SPgeno2+1
          SPgeno2[genSP[2,3] != genSP[2,4]]  <- SPgeno2+1
          SPgeno2[genSP[2,5] != genSP[2,6]]  <- SPgeno2+1
          SPgeno2[genSP[2,7] != genSP[2,8]]  <- SPgeno2+1
          SPgeno2[genSP[2,9] != genSP[2,10]]  <- SPgeno2+1
          SPgeno2[genSP[2,11] != genSP[2,12]]  <- SPgeno2+1

          # put these in the right place in table "summs"
          summs$noSPHetX1[summs$Ind==ind] <- SPgeno1
          summs$noSPHetX2[summs$Ind==ind] <- SPgeno2
          summs$noSPHetXtot[summs$Ind==ind] <- sum(SPgeno1,SPgeno2, na.rm=TRUE)

          # count total number of working SP-markers
          SPtgeno1 <- sum(genSP[1,]!=0, na.rm=TRUE)/2
          SPtgeno2 <- sum(genSP[2,]!=0, na.rm=TRUE)/2
          if(!replicates) SPtgeno2 <- NA

          # put these in the right place in table "summs"
          summs$noSPX1[summs$Ind==ind] <- SPtgeno1
          summs$noSPX2[summs$Ind==ind] <- SPtgeno2
          summs$noSPXtot[summs$Ind==ind] <- sum(SPtgeno1,SPtgeno2, na.rm=TRUE)
        }

      }
  }

    ### Export table with summary data for each sample
    ###################################################################
      write.table(summs, file = file.path(dirname,paste0(basename,"_summary_individuals.csv")), quote=FALSE, sep=";", row.names=FALSE, col.names=TRUE)
      summsRE <- summs[summs$categ=="RERUN",]
      write.table(summsRE$Ind, file = file.path(dirname,paste0(basename,"_samples_to_RERUN.txt")), quote=FALSE, sep=";", row.names=FALSE, col.names=FALSE)


    ### FINALY Create consensus PED for all "GOOD" samples AND export #######################################################
    # Ped file specifics:
    # -------------------
    # Family ID (col 1): here sample number
    # Individual ID (col 2): here sample name
    # Paternal ID (col 3): not used
    # Maternal ID (col 4): not used
    # Sex (col 5):1=male; 2=female; other=unknown
    # Condition (col 6): not used

    # Recode so that one locus (both alleles) is in one collumn
    # first take the ped file specific columns:
      ped_wo_y2 <- ped_wo_y[,1:6]
      n_marker <- (ncol(ped_wo_y)-6)/2
     # Add for each marker a factor column for genotype data (without Y-chromosome markers)
      for(k in 1:n_marker){
        ped_wo_y2[,k+6]<- as.factor("NA")
        levels(ped_wo_y2[,k+6])<-c("NA", "0 0",
                                "A A", "A C", "A G", "A T",
                                "C A", "C C", "C G", "C T",
                                "G A", "G C", "G G", "G T",
                                "T A", "T C", "T G", "T T")
      }

    # target column: n+6
    # take column 1: 6+2*(n-1)+1 = 7+2*(n-1)
    # take column 2: 6+2*(n-1)+2 = 8+2*(n-1)
      for(m in 1:nrow(ped_wo_y)){
        for(n in 1:n_marker){
          ped_wo_y2[m,n+6] <- as.factor(paste(ped_wo_y[m,7+2*(n-1)],ped_wo_y[m,8+2*(n-1)]))
        }
      }
      levels(ped_wo_y2$V5) <- c(levels(ped_wo_y2$V5), "0", "1", "2")
      ped_wo_y2[,5][is.na(ped_wo_y2[,5])]<- "0"

    # take out only samples with category "GOOD"
      goodInd <- summs$Ind[summs$categ == "GOOD"]

    # Create new dataframe for consensus sequence of GOOD samples
      GoodPed <- data.frame(matrix(as.factor("0 0"), nrow=length(unique(goodInd)), ncol=ncol(ped_wo_y2)), stringsAsFactors = TRUE)
      GoodPed$X1 <- as.integer(0)
      GoodPed$X3 <- as.integer(0)
      GoodPed$X4 <- as.integer(0)
      GoodPed$X5 <- as.integer(0)
      GoodPed$X6 <- as.integer(0)

    # Check if first and second run of sample have the same genotype, loop over all GOOD samples (j)
      for(j in 1:length(unique(goodInd))){
      # get all genotypes for individual j
        Gind <- unique(goodInd)[j]
        Gindgen <- ped_wo_y2[ped_wo_y2$V2==Gind,]
        Gindgen[,c(7:ncol(Gindgen))][is.na(Gindgen[,c(7:ncol(Gindgen))])]<- "0 0"

      # add factor levels to each of GoodPed factor columns
        levels(GoodPed[,2]) <- unique(c(levels(GoodPed[,2]), levels(Gindgen[,2])))
        for(m in 7:ncol(Gindgen)){
          levels(GoodPed[,m]) <- unique(c(levels(GoodPed[,m]), levels(Gindgen[,m])))
        }

      # check if both replicates are identical:
      # and if so add to consensus sequence
        GoodPed[j,2] <- Gindgen[1,2]
        if(nrow(Gindgen)>1){
          for(k in 7:ncol(Gindgen)){
            if(Gindgen[1,k]==Gindgen[2,k])
              GoodPed[j,k] <- Gindgen[1,k]
          }
        } else {
          for(k in 7:ncol(Gindgen)){
              GoodPed[j,k] <- Gindgen[1,k]
          }
        }
      # add family number
        GoodPed[j,1] <- j
      # add sex info (from summs)
        GoodPed[j,5] <- 0
        if(sexing){
          if(summs$sex[summs$Ind==Gind]=="Male")
            GoodPed[j,5] <- 1
          if(summs$sex[summs$Ind==Gind]=="Female")
            GoodPed[j,5] <- 2
        }
      }

    # export ped file
      ped.filename <- gsub(".ped", ".GOOD.ped", filename)
      write.table(GoodPed, file = file.path(dirname, ped.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

    # Create a database
      if(!is.na(db)){
        if(file.exists(file.path(dirname, db))){
          database <- read.table(file = file.path(dirname, db), sep="\t", header=FALSE)
          if(verbose>0) cat("Database file found, imported", nrow(db), "samples from the existing database\n")
          if(appendSamplesToDB){
            colnames(database) <- colnames(GoodPed)
            if(sum(is.element(GoodPed[,2], database[,2]))>(nrow(GoodPed)/2)) warning("There are lots of samples from current file already in the database, is that expected and wanted?\n")
            rownames(database) <- 1:nrow(database)
            rownames(GoodPed) <- (nrow(database)+1):(nrow(database)+nrow(GoodPed))
            database <- rbind(database, GoodPed)
            rownames(database) <- NULL
            write.table(database, file = file.path(dirname, db), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
          }
        } else {
          if(verbose>0) cat("Provided database file not found, create a new one with existing data!\n")
          database <- GoodPed
          write.table(database, file = file.path(dirname, db), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        }
      }

    # Create .map file without the Y-markers:
      if(!is.na(y.marker)){
        map.filename <- gsub(".ped", ".GOOD.map", filename)
        write.table(map1_wo_y, file = file.path(dirname, map.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      } else {
        map.filename <- gsub(".ped", ".GOOD.map", filename)
        write.table(map1, file = file.path(dirname, map.filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      }
      if(verbose>0) cat("\n ### Estimating errors: DONE! ",date(),"\n","##############################################################\n")

      output <- list(gensim=gensim,
                     summs=summs)
      output
}
