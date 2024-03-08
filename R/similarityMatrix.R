#' @title Calculate the similarity matrix
#'
#' @description
#' Performs pairwise similarity analysis on genotypic data.
#'
#' @param file Input file path. Default: NA.
#' @param mibs.file MIBS input file path. Default: NA.
#' @param pairs.file PAIRS input file path. Default: NA.
#' @param ped.file PED input file path. Default: NA.
#' @param group Sample identifier for statistics. Default: NA.
#' @param plots Should plots be created? Default: TRUE.
#' @param similarity Similarity threshold. Default: 0.85.
#' @param verbose Should output be verbose? Default: TRUE.
#' @param verbosity Verbosity level. Default: 1.
#'
#' @details
#' Reads genotype data, performs pairwise similarity calculations, generates plots, and outputs data for further analysis.
#'
#' @examples
#' \dontrun{
#'      similarityMatrix(file = file.path(outdir, "example_data.csv.GOOD"))
#' }
#'
#' @return Does not return a value. Creates output files in the same directory as the input files.
#' @export

similarityMatrix <- function(file=NA, mibs.file=NA, pairs.file=NA, ped.file=NA,
                             group=NA, plots=TRUE, similarity=0.85,
                             verbose=TRUE, verbosity=1){

  ## Verbose output of input parameters
  #################################################

  if (verbosity >= 2) {
    cat("Input parameters:\n")
    cat("file: ", file, "\n")
    cat("mibs.file: ", mibs.file, "\n")
    cat("pairs.file: ", pairs.file, "\n")
    cat("ped.file: ", ped.file, "\n")
    cat("group: ", group, "\n")
    cat("plots: ", plots, "\n")
    cat("similarity: ", similarity, "\n")
    cat("verbose: ", verbose, "\n")
    cat("verbosity: ", verbosity, "\n")
  }


  ## Input checks
  #################################################

  # Check if 'verbose' is logical
  if(!is.logical(verbose) || length(verbose) != 1){
    stop("'verbose' must be a single logical value.")
  }

  # Check if 'verbosity' is numeric
  if(!is.numeric(verbosity) || length(verbosity) != 1){
    stop("'verbosity' must be a single numeric value.")
  }

  if(!verbose & verbosity > 0) verbosity <- 0
  verbose <- verbosity
  if(is.na(file)) stop("ERROR: No input file provided")
  if(is.na(mibs.file)) mibs.file <- paste0(file,".mibs")
  if(is.na(pairs.file)) pairs.file <- paste0(file,".pairs")
  if(is.na(ped.file)) ped.file <- paste0(file,".ped")
  dirname <- dirname(file)


  # Check if 'file' is a character string
  if(!is.character(file) || length(file) != 1){
    stop("'file' must be a single character string representing the path to the main input file.")
  }

  # Check if 'mibs.file', 'pairs.file', and 'ped.file' are character strings
  if(!is.character(mibs.file) || length(mibs.file) != 1){
    stop("'mibs.file' must be a single character string representing the path to the MIBS input file.")
  }
  if(!is.character(pairs.file) || length(pairs.file) != 1){
    stop("'pairs.file' must be a single character string representing the path to the PAIRS input file.")
  }
  if(!is.character(ped.file) || length(ped.file) != 1){
    stop("'ped.file' must be a single character string representing the path to the PED input file.")
  }

  # Check if 'plots' is logical
  if(!is.logical(plots) || length(plots) != 1){
    stop("'plots' must be a single logical value.")
  }

  # Check if 'similarity' is numeric
  if(!is.numeric(similarity) || length(similarity) != 1){
    stop("'similarity' must be a single numeric value.")
  }

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
    # Suppress the warnings, see this discussion on it:
    # https://stackoverflow.com/questions/69666867/constant-warning-message-with-reshapemelt-in-r
    # Use suppressWarnings instead the restting of options
    ddant2 <- suppressWarnings({
      reshape::melt(ddant)[reshape::melt(upper.tri(ddant))$value,]
    })
    ddl$noLoci <- ddant2$value

  # get the number of loci typed for each sample (diagonal from ddant)
    samples <- data.frame(matrix(NA, nrow=nrow(mat), ncol=2))
    names(samples) <- c("sample", "noLoci")
    samples$sample <- nam
    samples$noLoci <- diag(ddant)

  # check and plot the data
    if(plots){
      fig1.filename <- paste0(file, ".pairwise_similarities.png")

    # Restore old par settings, when function exists
      oldpar <- par(no.readonly = TRUE)
      on.exit({
        oldpar$new <- NULL
        par(oldpar)
      })


      png(fig1.filename, width=1200, height=800)
        par(mfrow=c(1,2))
        hist(ddl$similarity, breaks=50, main="", xlab="Pairwise similarities")
        hist(ddl$similarity, xlim=c(0.75,1), ylim=c(0,200), breaks=50, main="", xlab="Pairwise similarities (axes restricted)")
      dev.off()
    }

  # output all pairwise similarities over 85%
    ddsim <- ddl[ddl$similarity>=similarity,]
    write.table(ddsim, file=paste0(file,paste0(".similar_",similarity,".genotypes.rout") ), quote=FALSE, sep=",", row.names=FALSE)

  if(!is.na(group)){
    ### RUN BELOW THIS ONLY IF YOU WANT SUMMARY OUTPUT FOR EACH SAMPLE SEPARATELLY ###********************************************
    ### check each (i) of the new samples individually against both the current run and the Database
    # Take all levels of "nam" beginning with grupp (specified above)
    newsamp <- as.data.frame(grep(group, nam, value=TRUE))
    names(newsamp) <- c("samp")
    levels(newsamp$samp) <- c(levels(newsamp$samp),levels(ddl$sample1))

    # Create output file using a loop over all samples:
    writeLines(paste0("Database matching for samples in file: ", name), paste0(file,".individual.rout"))
    cat("****************************************************************", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
    cat(" ", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)

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
      cat(paste("Matches for sample: ", newsamp[i,1], sep=""), file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat(paste("Number of loci typed in this sample: ", NL, sep=""), file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat(paste("Matches against same run: ", length(grep(group, sim$sample2, value=TRUE)), sep=""), file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat(paste("Matches against DataBase: ", length(grep(group, sim$sample2, value=TRUE, invert=TRUE)), sep=""), file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat(" ", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      write.table(sim,file=paste0(file,".individual.rout"), append=TRUE, quote=FALSE, row.names=FALSE, sep=" ", col.names=FALSE)
      cat(" ", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)

      if(nrow(warn)>0){
        cat("WARNING, similarities between 0.85 and 0.95:", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)

        write.table(warn, file=paste0(file,".individual.rout"), append=TRUE, quote=FALSE, row.names=FALSE, sep=" ", col.names=FALSE)
      }

      cat(" ", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat("-----------------------------------------------------------------", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
      cat(" ", file=paste0(file,".individual.rout"),sep="\n", append=TRUE)
    }
  }

    if(verbose>0) message("\n ### Similarity Matrix and groups: DONE! ",date(),"\n","##############################################################\n")
}
