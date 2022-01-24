.onAttach <- function(libname, pkgname){
  packageStartupMessage("This is Fluidigm version ", utils::packageVersion("Fluidigm"))
}
