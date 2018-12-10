if (!isGeneric("hist"))
  setGeneric("hist",
             function(x, ...) standardGeneric("hist"))

if (!isGeneric("plot"))
  setGeneric("plot",
             function(x, y, ...) standardGeneric("plot"))

if (!isGeneric("summary"))
  setGeneric("summary",
             function(object, ...) standardGeneric("summary"))

