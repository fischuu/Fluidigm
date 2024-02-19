# Fluidigm
Fluidigm data analysis r-package, a more detailed description and vignette with
corresponding example data is on its way.

# Installation
You can install the package directly from GitHub using `devtools` like this

````{r}
devtools::install_github("fischuu/fluidigm")
```

# Further requirements
Please note, that you need to have Plink installed on your system and it needs to be available on the PATH.

# Running example

```{r}
# Load the library
  library("Fluidigm")
  
# Set the working directory
  setwd("~/tmp/fixing_sex_assign_fluiRscript")

# Define the required files
  file <- "Run22_1802150005_JH.csv"
  map <- "Run22_1802150005_JH.map"
  db <- "Run22_1802150005_JH.ped"
  neg_controls=c("STA-blank", "Chipblank")
  
# Define y and x markers
  y.marker <- "DBY7"
  x.marker <- c("BICF2G630532567",
                "BICF2P345621",
                "BICF2P354225",
                "BICF2S23212310",
                "BICF2P10954231",
                "BICF2S23041234",
                "BICF2G630532123",
                "BICF2G63054321",
                "BICF2P1922342")

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

[plink-pairwise-loci.pl](https://github.com/douglasgscofield/bioinfo/blob/main/scripts/plink-pairwise-loci.pl)
