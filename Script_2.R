
#define Files for input 
def=data.frame(dir())
def$def=c("Basline", "Glutamate","Post_Glutatate", "EDTA")
names(def)=c("files", "def")

#shape file

shape=data.frame(chanel=seq(1,50, by=1), tag=rep(50,50), thresold=0)


# Set up new class 

MEA= setClass("MEA", slots = c(path="character", MAT_Basline = "data.frame", MAT_Tr = "data.frame", Aligned_Zero_base="data.frame",Aligned_Zero_TR="data.frame", Shape="data.frame", Peaks="list",IPI="list", Raster="data.frame"))
setValidity("MEA", function(object) {
  msg <- NULL
  if ( ! is.character(object@path) ){
  msg <- c(msg, "input data must be data.frame")}
  if (is.null(msg)) TRUE
  else msg
})         
setMethod("initialize", signature = "MEA", definition = function(.Object, path, def, TR){
  require(R.matlab)
  #look up for basline
  path_baseline=paste(path,"/", as.character(def[def$def=="Basline", ]$files), sep="")
  .Object@MAT_Basline <- as.data.frame(readMat(path_baseline))
  path_TR=paste(path,"/", as.character(def[def$def==TR, ]$files), sep="")
  .Object@MAT_Tr <- as.data.frame(readMat(path_TR))
  validObject(.Object)
  return(.Object)})

path="~/Desktop/Proj_Charlotta/Tumor Project-MEA Recordings/Cell Line 233/20180809"


#Imut a mat file into the MEA class 
object=MEA(path, def, TR="Glutamate")

object@Shape=shape

setGeneric("Align_Zero", function(object,Edges=1000) standardGeneric("Align_Zero"))
setMethod("Align_Zero",signature = "MEA", definition = function(object,Edges=1000){
  
  #Basline
  mat=object@MAT_Basline
  dim(mat)
  #Cut of the Edges
  mat=mat[Edges:nrow(mat)-Edges, ]
  #Normal Values
  mat=mat/1000000
  
  rownames(mat)=mat[,1]
  mat=mat[, 2:ncol(mat)]
  dim(mat)

  #Align to Zero
  mat_zero=mat
  
  for(i in 1:ncol(mat_zero)){
    mean=mean(mat_zero[,i])
    mat_zero[,i]=(mat_zero[,i]/mean)-1
  }
  
  object@Aligned_Zero_base=mat_zero
  
  
  
  #Tretament
  mat=object@MAT_Tr
  dim(mat)
  #Cut of the Edges
  mat=mat[Edges:nrow(mat)-Edges, ]
  #Normal Values
  mat=mat/1000000
  
  rownames(mat)=mat[,1]
  mat=mat[, 2:ncol(mat)]
  dim(mat)
  
  #Align to Zero
  mat_zero=mat
  
  for(i in 1:ncol(mat_zero)){
    mean=mean(mat_zero[,i])
    mat_zero[,i]=(mat_zero[,i]/mean)-1
  }
  
  object@Aligned_Zero_TR=mat_zero
  
  
  
  
  return(object)
  
  })


object=Align_Zero(object,Edges=1000)


setGeneric("Calculate_Threshold", function(object) standardGeneric("Calculate_Threshold"))
setMethod("Calculate_Threshold",signature = "MEA", definition = function(object){
  
  #First get the basline Threshold
  
  mat=object@Aligned_Zero_base
  dim(mat)
  
  threshold=sd(mat[1:1000,1])*3
  
  for(i in 2:ncol(mat)){
    sd=sd(mat[1:1000,i])*3
    threshold=c(threshold,sd)
  }
  
  object@Shape$thresold=threshold
  
  return(object)
  
})


object=Calculate_Threshold(object)


setGeneric("Finde_Peaks", function(object) standardGeneric("Finde_Peaks"))
setMethod("Finde_Peaks",signature = "MEA", definition = function(object){
  require(pracma)
  require(parallel)
  #First get the basline Threshold
  
  mat=object@Aligned_Zero_TR
  dim(mat)
  out=mclapply(1:ncol(mat),mc.cores=4, function(i){
    vector=mat[,i]
    peaks=findpeaks(vector,minpeakheight=object@Shape$thresold[i])
    return(peaks)
    })
  
  object@Peaks=out
  
  return(object)
  
})

object=Finde_Peaks(object)

setGeneric("Find_IPI", function(object,samp_freq) standardGeneric("Find_IPI"))
setMethod("Find_IPI",signature = "MEA", definition = function(object,samp_freq){
  
  time=data.frame(Time=as.numeric(rownames(object@MAT_Tr))/samp_freq, dat_points=seq(1:length(rownames(object@MAT_Tr))))
  dim(time)
  
  
  peaks=object@Peaks
  
  nr=as.numeric(length((peaks)))
  
  out=as.numeric(unlist(lapply(1:nr, function(i){
    time_peaks=time[time$dat_points %in% peaks[[i]][,3], ]$Time
    IPI=time_peaks[2]-time_peaks[1]
    for(ii in 2:length(time_peaks)-1){IPI=c(IPI, (time_peaks[ii+1]-time_peaks[ii]) ) }
    return(mean(IPI))
    })))
  out2=lapply(1:nr, function(i){
    
    time_peaks_x=time[time$dat_points %in% peaks[[i]][,3], ]
    if(nrow(time_peaks_x)>2){
      time_peaks_x=time[time$dat_points %in% peaks[[i]][,3], ]
      time_peaks=time[time$dat_points %in% peaks[[i]][,3], ]$Time
      IPI=time_peaks[2]-time_peaks[1]
      for(ii in 2:length(time_peaks)-1){IPI=c(IPI, (time_peaks[ii+1]-time_peaks[ii]) ) }
      
      time_peaks_x$IPI=as.numeric(IPI)
      
      
      print(i)
      
      if(is.na(out[i])==F){
        if (out[i]<1 & max(peaks[[i]][,1])>100){plot(x=time_peaks_x$IPI, y= peaks[[i]][,1], pch=19, bty="n",
                           main=paste("Channel:",i, sep=" "),
                           xlab="IPI in sec", ylab="Peak Intensity")}
      }
      
      
      
      return(time_peaks_x$IPI)
    } else {return(NA)}
    
    
    
  })
  

  object@IPI=out2
  object@Shape$IPI=out

  
  return(object)
  
})


object=Find_IPI(object,samp_freq=1000)

setGeneric("Plot_Raster", function(object) standardGeneric("Plot_Raster"))
setMethod("Plot_Raster",signature = "MEA", definition = function(object){
  
 
  
  return(object)
  
})


object=Plot_Raster(object)




object@Shape

plot(object@Aligned_Zero_TR[,20], type="l")

library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
k <- kde2d(na.omit(object@Shape$thresold), na.omit(object@Shape$IPI),n=200)
image(k, col=r)
















