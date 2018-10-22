### The first three functions are here for backwards compatibility.
### They are officially "deprecated", which means (a) you should not
### use them in new code because (b) at some point theu will simply
### go away.
tumorGen <- function(...) {
  tumor <- Tumor(...)
  as(tumor, "list")
}

dataGen <- function(tumor, ...) {
  generateTumorData(as(tumor, "Tumor"), ...)
}

runAlg <- function(...) findClones(...)
