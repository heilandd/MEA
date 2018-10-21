
# Set up new class 

setClass("MEA", slots = c(path="character", MAT = "array", Aligned_Zero="array", Shape="data.frame", Peaks="data.frame", Raster="data.frame"))
setValidity("MEA", function(object) {
  msg <- NULL
  if ( ! is.character(object@path) ){
  msg <- c(msg, "input data must be data.frame")}
  if (is.null(msg)) TRUE
  else msg
})         
setMethod("initialize", signature = "MEA", definition = function(.Object, path){
  require(R.matlab)
  .Object@path <- path
  validObject(.Object)
  return(.Object)})

path="/media/daka/Data Intern/Henrik/ElectroPhys/Tumor Project-MEA Recordings/Records_3/233_1.mat"

object=MEA(path)
