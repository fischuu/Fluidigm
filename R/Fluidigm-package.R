#' Fluidigm
#'
#' @section Title:
#' Comprehensive Analysis of Fluidigm Genotyping Data
#'
#' @section Description:
#' A suite of tools designed to streamline the process of analyzing genotyping data from Fluidigm machines.
#' It includes functions for converting Fluidigm data to 'PLINK' format, estimating errors, calculating pairwise
#' similarities, determining pairwise similarity loci, and generating a similarity matrix.
#'
#' @section Details:
#' The package provides a comprehensive analysis pipeline for Fluidigm genotyping data. It starts by converting
#' the raw data from the Fluidigm machine into a format that can be used with the 'PLINK' software. It then estimates
#' errors in the data, calculates pairwise similarities between genotypes, determines pairwise similarity loci,
#' and generates a similarity matrix. The package is designed to make it easier and more efficient for researchers
#' to extract meaningful insights from their genotyping studies.
#'
#' @section References:
#' \itemize{
#'   \item{Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
#'   \emph{PLINK: a toolset for whole-genome association and population-based linkage analysis.} American Journal of Human Genetics, 81.}
#'   \item{Purcell, Shaun. PLINK. \url{https://zzz.bwh.harvard.edu/plink/}}}
#'
#' @docType package
#' @name Fluidigm
#' @keywords fluidigm genotyping
"_PACKAGE"

## usethis namespace: start
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom stats cor.test
#' @importFrom utils read.table
#' @importFrom utils write.table
## usethis namespace: end
NULL
