# Fluidigm
Fluidigm data analysis r-package, a more detailed description and vignette with
corresponding example data is on its way.

# Installation
You can install the package directly from GitHub using `devtools` like this

```{r}
devtools::install_github("fischuu/fluidigm")
```

# Further requirements
Please note, that you need to have Plink installed on your system and it needs to be available on the PATH.

# Running example

```{r}
# Load the library
  library("Fluidigm")
  
# Set the working directory
  setwd("~/Project/My_fluidigm_project")

# Define the required files
  file <- "Run.csv"
  map <- "Run.map"
  db <- "Run.ped"
  neg_controls=c("STA-blank", "Chipblank")
  
# Define y and x markers
  y.marker <- "markerY1"
  x.marker <- c("markerX1",
                "markerX2",
                "markerX3",
                "markerX4",
                "markerX5",
                "markerX6",
                "markerX7",
                "markerX8",
                "markerX9")

# Run the analysis
  out <- fluidigmAnalysisWrapper(file=file,
                                 db=db,
                                 map=map,
                                 neg_controls=neg_controls,
                                 y.marker = y.marker,
                                 x.marker = x.marker)
  
# You can also run the individual steps  

# First, turn the fluidigm output into Plink format
  fluidigm2PLINK(file=NA, out=NA, map=NA, plots=TRUE, rearrange=TRUE, verbose=TRUE)
```


# Acknowledgements
The pacakge uses also a perl script that was written by Doug Scofield and which is published here:

https://github.com/douglasgscofield/bioinfo/blob/main/scripts/plink-pairwise-loci.pl
