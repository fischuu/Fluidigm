<div align=“right”> <img src=“logo.jpeg” alt=“Fluidigm logo”> </div>
# Fluidigm
Fluidigm is an R package that can analyze single-cell genotyping data from Fluidigm instruments. It can perform various operations such as:

* Converting the raw data into PLINK format
* Estimating genotyping errors
* Calculating pairwise similarities
* Determining pairwise similarity loci
* Generating a similarity matrix

Fluidigm provides a convenient interface to the powerful capabilities of the PLINK software, which is a free, open-source whole genome association analysis toolset. It also provides detailed output files and plots to help you explore and understand your genotypic data.

# Installation
You can install Fluidigm from CRAN with (soon...):

```
install.packages("Fluidigm")
```
Or you can install the development version from GitHub with:

```
# install.packages("devtools")
devtools::install_github("fischuu/Fluidigm")
```

# Further requirements
Please note, that you need to have Plink installed on your system and it needs to be available on the PATH variable.

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

fluidigm2PLINK(...)
estimateErrors(...)
calculatePairwiseSimilarities(...)
getPairwiseSimilarityLoci(...)
similarityMatrix(...)

```

# License
Fluidigm is licensed under the GPL-3 license. See the LICENSE file for more information.

# Citation
If you use Fluidigm in your research, please cite it as follows:

```
citation("Fluidigm")
```

Contact
If you have any questions, suggestions, or feedback, please feel free to contact me. I would love to hear from you and improve the package. 🙌

# Acknowledgements
The package uses also a perl script that was written by Doug Scofield and which is published here:

https://github.com/douglasgscofield/bioinfo/blob/main/scripts/plink-pairwise-loci.pl
