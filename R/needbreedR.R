.onAttach <- function(libname, pkgname)
{ if(!("breedR" %in% loadedNamespaces()))
  { breedR.loaded <- requireNamespace("breedR", quietly=TRUE)
    if (!breedR.loaded)  
    { packageStartupMessage("breedR needs to be loaded if the mixed-model functions are to be used.

breedR is available from Github. Please visit https://github.com/famuvie/breedR.\n")}
  }
}

