
#define Files for input 
def=data.frame(dir())
def$def=c("Basline", "Post_Glutatate", "Glutamate", "EDTA")
names(def)=c("files", "def")

#Example
# def
# files            def
# 1     1.mat        Basline
# 2 233_1.mat Post_Glutatate
# 3     2.mat      Glutamate
# 4     4.mat           EDTA
 



#shape file
#Could also include info of coordinates and so on 
shape=data.frame(chanel=seq(1,50, by=1), tag=rep(50,50), thresold=0)



# Set up new class 

MEA= setClass("MEA", slots = c(path="character", MAT_Basline = "data.frame", MAT_Tr = "data.frame", 
                               baselinefile="data.frame", 
                               Aligned_Zero_base="data.frame",Aligned_Zero_TR="data.frame", 
                               Shape="data.frame", Peaks="list", Raster="data.frame"))



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

path="/media/daka/Data Intern/Henrik/ElectroPhys/Tumor Project-MEA Recordings/Cell Line 233/20180809"


#Imut a mat file into the MEA class 
object=MEA(path, def, TR="Post_Glutatate")

object@Shape=shape

setGeneric("Align_Zero", function(object,STD=10,Edges=1000,Exp_start=3000,Exp_end=40000, Exp_chanel=21,plot=T,n=10000) standardGeneric("Align_Zero"))
setMethod("Align_Zero",signature = "MEA", definition = function(object,STD=10,Edges=1000,Exp_start=1000,Exp_end=200000, Exp_chanel=21,plot=T,n=5000){
  
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
  
  
  
  
  
  mat_zero=mat
  
  #Align to Zero
  
  for(i in 1:ncol(mat_zero)){
    mean=mean(mat_zero[,i])
    mat_zero[,i]=(mat_zero[,i]/mean)-1
  }

  
  
  
  
  #call Basline
  

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
  if(plot==T){
    layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
    plot(mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l", ylab="", xlab="Raw Data")
    
  }
  
  for(i in 1:ncol(mat_zero)){
    mean=mean(mat_zero[,i])
    mat_zero[,i]=(mat_zero[,i]/mean)-1
  }
  
  
  if(plot==T){
    plot(mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l",ylab = "", xlab="Data Aligned to Zero")
    
  }
  
  #Define baseline
  
  
  require(parallel)
  require(pracma)
  baselinefile=as.data.frame(do.call(cbind, mclapply(1:ncol(mat_zero),mc.cores = 8, function(i){
    print(paste("Check for Chanel:",i, table(is.na(t(mat_zero[,i,drop=F])))))
    #mean=IDPmisc::baseli (x=as.numeric(rownames(mat_zero)), y=mat_zero[,i],span=1000, NoXP=1000, maxit = c(1, 1))
    #mean=baseline::baseline.modpolyfit(t(mat_zero[,i, drop=F]),degree = 4, tol = 0.001, rep = 100)
    mean=movavg(as.numeric(mat_zero[,i]), n=n, type=c("t"))
    
    
    return(mean)
  })))
  
  if(plot==T){
    plot(x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), y=mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l",ylab = "", xlab="Add Baseline")
    points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel], col="red")
    
  }
  
  object@baselinefile=baselinefile
  object@Aligned_Zero_TR=mat_zero
  
  mat=object@Aligned_Zero_base
  dim(mat)
  
  threshold=sd(mat[1:1000,1])*STD
  
  for(i in 2:ncol(mat_zero)){
    sd=sd(mat[1:1000,i])*3
    threshold=c(threshold,sd)
  }
  
  object@Shape$thresold=threshold
  
  
  if(plot==T){
    
    plot(x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), y=mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l",ylab = "", xlab="Addaped Threshold")
    points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel]+threshold[Exp_chanel], col="cyan", lwd=2)
    points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel]-threshold[Exp_chanel], col="cyan", lwd=2)
  }
  
  
  return(object)
  
})



object=Align_Zero(object,STD=10,Edges=1000,Exp_start=100000,Exp_end=500000, Exp_chanel=21,plot=T)


setGeneric("plot_chanels", function(object,Exp_start=100000,Exp_end=500000, Exp_chanel=21) standardGeneric("plot_chanels"))
setMethod("plot_chanels",signature = "MEA", definition = function(object,Exp_start=100000,Exp_end=50000, Exp_chanel=21){
  
  mat_zero=object@Aligned_Zero_base
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  plot(mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l", ylab="", xlab="Raw Data")
  
  mat_zero=object@Aligned_Zero_TR
  baselinefile=object@baselinefile
  threshold= object@Shape$thresold
  plot(x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), y=mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l",ylab = "", xlab="Add Baseline")
  points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel], col="red")
  plot(x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), y=mat_zero[Exp_start:Exp_end, Exp_chanel], axes=F, type="l",ylab = "", xlab="Addaped Threshold")
  points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel]+threshold[Exp_chanel], col="cyan", lwd=2)
  points(type="l",x=as.numeric(rownames(mat_zero)[Exp_start:Exp_end]), baselinefile[Exp_start:Exp_end, Exp_chanel]-threshold[Exp_chanel], col="cyan", lwd=2)
  
  #return(object)
  
})

plot_chanels(object, Exp_chanel=26,Exp_start=9000,Exp_end=1199000)


setGeneric("Finde_Peaks", function(object) standardGeneric("Finde_Peaks"))
setMethod("Finde_Peaks",signature = "MEA", definition = function(object){
  require(spectral)
  require(parallel)
  #First get the basline Threshold
  
  mat=object@Aligned_Zero_TR
  dim(mat)
  out=mclapply(1:ncol(mat),mc.cores=14, function(i){
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
  require(spectral)
  require(parallel)
  
  
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
  
  object@Shape$IPI=out
  
  
  return(object)
  
})


object=Find_IPI(object,samp_freq=1000)

saveRDS(object, "object_233.R")

##Go further with the time plot 


#input timings
timings=readRDS("object_233.R")
layout=read.csv("Layout.csv", sep=";")





setGeneric("Plot_Events_Scatter", function(object,RasterP,layout) standardGeneric("Plot_Events_Scatter"))
setMethod("Plot_Events_Scatter",signature = "MEA", definition = function(object,RasterP=T,layout){
  
  #read in events
  layout=layout
  events=object@Peaks
  
  events=events[!sapply(events, is.null)]
  
  
  if (RasterP==T){
    #create raster plot
    
    plot(y=c(1,length(events)), x=c(0,500000), bty="n", type="n")
    
    for(i in 1:length(events)){
      x=as.numeric(na.omit(events[[i]][,2]))
      y=rep(i,length(x))
      print(i)
      points(x=x, y=y, cex=0.8, lwd=0.2, pch="|")
      
      
    }
  }
  
  
  plot(c(-2,11), c(-2,11), bty="n", type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  require(RColorBrewer); rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(length(events))
  r=sample(r)
  
  for (i in 1:length(events)){
    
    #start chanels
    Intensity=(as.numeric(na.omit(events[[i]][,1])))/100
    Time=(as.numeric(na.omit(events[[i]][,2])))/100000000
    #Intensity=rnorm(length(as.numeric(na.omit(events[[i]][,1]))))/1.5
    #Intensity=exp(min(Intensity)/Intensity)/2
    
    #set points
    x=layout$x[i]; y=layout$y[i]
    
    
    for (j in 1:length(Intensity)){
      print(j)
      
      #define Angle  
      
      alpha=c((((1/360)*length(Time))*j)*Time[j])
      
      
      
      
      xx=sin(alpha)*Time[j]
      
      yy=cos(alpha)*Time[j]
      
      points(x+(xx*(1/Intensity[j])),y+(yy*(1/Intensity[j])), pch=19, cex=0.2, col=r[i])
      
    }
    
  }
  
  points(x=layout$x, y=layout$y, pch=19, cex=1.5, col="black")
  
  
  
  #Add Channels
  
  y_cor=seq(10,-2, length.out=length(events))
  for (i in 1:length(events)) {
    
    y=y_cor[i]
    x=-2
    
    polygon(x=c(x,x+0.2), y=c(y,y), border =r[i], lwd=4)
    text(x=x+0.8, y=y, labels = paste("Channel: ",i,sep="" ), cex = 0.5)
    
  }
  
  
  # Compute Connections
  
  #sorted DB with all events
  DB_Events=as.data.frame(do.call(rbind, lapply(1:length(events), function(i){out=data.frame(Event=as.numeric(na.omit(events[[i]][,2])), Channel=paste("Channel_", i , sep=""))})))
  #Add coordinates to each event
  layout$Channel_new=paste("Channel_",layout$Channel,sep="")
  
  DB_Events$x=1
  DB_Events$y=1
  
  for(i in 1:nrow(DB_Events)){DB_Events[i,]$x=layout[layout$Channel_new==DB_Events[i,]$Channel, ]$x}
  for(i in 1:nrow(DB_Events)){DB_Events[i,]$y=layout[layout$Channel_new==DB_Events[i,]$Channel, ]$y}
  DB_Events=DB_Events[order(DB_Events$Event, decreasing = F), ]
  connection=na.omit(as.data.frame(do.call(rbind, lapply(2:nrow(DB_Events),function(i){
    
    if(DB_Events[i,]$Event-DB_Events[i-1,]$Event<1 ){
      from=DB_Events[i-1,]$Channel
      to=DB_Events[i,]$Channel
      
      if(from!=to){out=data.frame(from=from, to=to)}else{print("No Connection")}
    }else{print("No Connection")}
    
    
  }))))
  
  
  #count the number of connections
  
  from_x=connection[!duplicated(connection$from), ]$from
  connection_counts=as.data.frame(do.call(rbind, lapply(1:length(from_x),function(i){
    
    a_x=connection[connection$from==from_x[i], ]
    
    to_x=as.character(a_x[!duplicated(a_x$to), ]$to)
    
    out_x=as.data.frame(do.call(rbind, lapply(1:length(to_x),function(j){
      print(to_x[j])
      
      a_xx=as.numeric(nrow(a_x[a_x$to==to_x[j], ]))
      
      out_xx=data.frame(from=from_x[i], to=to_x[j], Nr_con=a_xx)
      
      
      
    })))
    
    
    
  })))
  connection_counts=connection_counts[order(connection_counts$Nr_con, decreasing = T), ]
  
  
  
  #remove only one connection
  connection_counts=connection_counts[connection_counts$Nr_con>4, ]
  
  
  #plot(x=layout$x, y=layout$y, pch=19, cex=1.5, col="black", bty="n", yaxt="n", xaxt="n", xlab="", ylab="")
  for(i in 1:nrow(connection_counts)){
    xx=c(layout[layout$Channel_new==connection_counts[i,]$from, ]$x,layout[layout$Channel_new==connection_counts[i,]$to, ]$x)
    yy=c(layout[layout$Channel_new==connection_counts[i,]$from, ]$y,layout[layout$Channel_new==connection_counts[i,]$to, ]$y)
    
    polygon(x=xx,y=yy, lwd=connection_counts[i, ]$Nr_con/15)
    
  }
  
  
  #object@connections
  
  
  
  
  return(object)
  
})



Plot_Events_Scatter(object,RasterP=T,layout)



#Full Pipe:

path="/media/daka/Data Intern/Henrik/ElectroPhys/Tumor Project-MEA Recordings/Cell Line 233/20180809"


#Imut a mat file into the MEA class 
object=MEA(path, def, TR="Post_Glutatate")
object@Shape=shape
object=Align_Zero(object,STD=10,Edges=1000,Exp_start=100000,Exp_end=500000, Exp_chanel=21,plot=T)
#plot_chanels(object, Exp_chanel=26,Exp_start=9000,Exp_end=1199000)
object=Finde_Peaks(object)
object=Find_IPI(object,samp_freq=1000)

#saveRDS(object, "object_233.R")
#layout=read.csv("Layout.csv", sep=";")

Plot_Events_Scatter(object,RasterP=T,layout)


