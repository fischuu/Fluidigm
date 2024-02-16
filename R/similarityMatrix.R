#' @title Run the fluidigm analysis script together
#'
#' @description
#' This function is a wrapper for the whole analysis
#'
#' @param file Path to the input file
#' @param mibs.file Path to the mibs input file (only if name differs)
#' @param pairs.file Path to the pairs input file (only if name differs)
#' @param ped.file Path to the pairs input file (only if name differs)
#' @param group Sample identified for sample-wise statistics
#' @param plots logical, shall plots be created?
#' @param similarity  Threshold defining the level of similarity
#' @param verbose Should the output be verbose, logical
#' @param verbosity Level of verbosity, set to higher number for more details
#'
#' #' @details
#' Additional details...
#'
#' @examples
#' \dontrun{
#'   similarityMatrix()
#' }
#'
#' @return A list containing the following elements:
#'         gensim, a matrix indicating if genotypes are called correctly for replicates and/or if genotypes are missing
#'         summs, a matrix with summary statistics
#' @export

similarityMatrix <- function(file=NA, mibs.file=NA, pairs.file=NA, ped.file=NA, group=NA, plots=TRUE, similarity=0.85, verbose=TRUE, verbosity=1){

  if(!verbose & verbosity > 0) verbosity <- 0
  verbose <- verbosity
  if(is.na(file)) stop("ERROR: No input file provided")
  if(is.na(mibs.file)) mibs.file <- paste0(file,".mibs")
  if(is.na(pairs.file)) pairs.file <- paste0(file,".pairs")
  if(is.na(ped.file)) ped.file <- paste0(file,".ped")
  dirname <- dirname(file)
  #####################################################################################################

  # Import similarity matrix
    name <- basename(mibs.file)
    #name1 <- sub("MATCH.mibs","",name)
    mat <- read.table(mibs.file)

  # Import sample names:
    ped <- read.table(ped.file)
    nam <- ped[,2]
    names(mat) <- nam
    if(sum(is.element(names(table(table(nam))), "1")) && length(is.element(names(table(table(nam))), "1"))==1 ){
      row.names(mat) <- nam
    } else {
      row.names(mat) <- paste(nam, rownames(mat),sep="_row")
    }


  # convert from matrix to long format
    dd <- as.matrix(mat)
  # Switch off the warnings, see this discussion on it:
  # https://stackoverflow.com/questions/69666867/constant-warning-message-with-reshapemelt-in-r
    oldw <- getOption("warn")
    options(warn = -1)
      ddl <- reshape::melt(dd)[reshape::melt(upper.tri(dd))$value,]
    options(warn = oldw)
    names(ddl) <- c("sample1", "sample2", "similarity")
    levels(ddl$sample1) <- c(levels(ddl$sample1), levels(ddl$sample2))
    levels(ddl$sample2) <- levels(ddl$sample1)

  # Import data on number of loci for each pairwise comparison
    ant <- read.table(pairs.file)

    ddant <- as.matrix(ant)
    # Switch off the warnings, see this discussion on it:
    # https://stackoverflow.com/questions/69666867/constant-warning-message-with-reshapemelt-in-r
    oldw <- getOption("warn")
    options(warn = -1)
      ddant2 <- reshape::melt(ddant)[reshape::melt(upper.tri(ddant))$value,]
    options(warn = oldw)
    ddl$noLoci <- ddant2$value

  # get the number of loci typed for each sample (diagonal from ddant)
    samples <- data.frame(matrix(NA, nrow=nrow(mat), ncol=2))
    names(samples) <- c("sample", "noLoci")
    samples$sample <- nam
    samples$noLoci <- diag(ddant)

  # check and plot the data
    if(plots){
      fig1.filename <- paste0(file, ".pairwise_similarities.png")
      png(fig1.filename, width=1200, height=800)
        par(mfrow=c(1,2))
        hist(ddl$similarity, breaks=50, main="", xlab="Pairwise similarities")
        hist(ddl$similarity, xlim=c(0.75,1), ylim=c(0,200), breaks=50, main="", xlab="Pairwise similarities (axes restricted)")
      dev.off()
    }

  # output all pairwise similarities over 85%
    ddsim <- ddl[ddl$similarity>=similarity,]
    write.table(ddsim, file=paste0(file,paste0("_similar_",similarity,"_genotypes.rout") ), quote=FALSE, sep=",", row.names=FALSE)

  if(!is.na(group)){
    ### RUN BELOW THIS ONLY IF YOU WANT SUMMARY OUTPUT FOR EACH SAMPLE SEPARATELLY ###********************************************
    ### check each (i) of the new samples individually against both the current run and the Database
    # Take all levels of "nam" beginning with grupp (specified above)
    newsamp <- as.data.frame(grep(group, nam, value=TRUE))
    names(newsamp) <- c("samp")
    levels(newsamp$samp) <- c(levels(newsamp$samp),levels(ddl$sample1))

    # Create output file using a loop over all samples:
    writeLines(paste0("Database matching for samples in file: ", name), paste0(file,"_individual.rout"))
    cat("****************************************************************", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
    cat(" ", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)

    for(i in 1:nrow(newsamp)){
      # take all cases where sample 1 or 2 = current individual
      simil <- ddl[ddl$sample1==newsamp[i,1]|ddl$sample2==newsamp[i,1],]
      # rearange so that focal individual is always in first column (sample1)
      simil$sample2[simil$sample2==newsamp[i,1]] <- simil$sample1[simil$sample2==newsamp[i,1]]
      simil$sample1[simil$sample1==simil$sample2] <- newsamp[i,1]

      sim <- simil[simil$similarity >= 0.95,]
      warn <- simil[simil$similarity >= 0.85,]   # Changed from 0.8 to 0.85 20180711 #
      warn <- warn[warn$similarity < 0.95,]
      if(nrow(warn)>0){
        warn$stars <- "************"	# Added column of stars to warn 20180711 #
      }

      # get info on number of loci for focal sample
      NL <- samples$noLoci[samples$sample==newsamp[i,1]]

      # take rows where match is against own run and against database
      # Add results to output
      cat(paste("Matches for sample: ", newsamp[i,1], sep=""), file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat(paste("Number of loci typed in this sample: ", NL, sep=""), file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat(paste("Matches against same run: ", length(grep(group, sim$sample2, value=TRUE)), sep=""), file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat(paste("Matches against DataBase: ", length(grep(group, sim$sample2, value=TRUE, invert=TRUE)), sep=""), file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat(" ", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      write.table(sim,file=paste0(file,"_individual.rout"), append=TRUE, quote=FALSE, row.names=FALSE, sep=" ", col.names=FALSE)
      cat(" ", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)

      if(nrow(warn)>0){
        cat("WARNING, similarities between 0.85 and 0.95:", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)

        write.table(warn, file=paste0(file,"_individual.rout"), append=TRUE, quote=FALSE, row.names=FALSE, sep=" ", col.names=FALSE)
      }

      cat(" ", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat("-----------------------------------------------------------------", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
      cat(" ", file=paste0(file,"_individual.rout"),sep="\n", append=TRUE)
    }
  }

    if(verbose>0)cat("\n ### Similarity Matrix and groups: DONE! ",date(),"\n","##############################################################\n")
}
