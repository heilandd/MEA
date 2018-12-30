# R_MEA: An R-Tool for Comprehensive Analysis of Multi-Electrode Array's 
 

The pipeline was designed to integrate following aspects:

Reproducibility: Analysis needs to be easily reproduced by external researchers.

Easy-to-Use: The pipeline was designed to be user-friendly and applicable for non-expert users.

Compatible: The pipeline should be feasible for multiple designs of MEA layouts



### Description

The full pipeline could be used by following commands:


```
# Define your data structures with ctr and treatment conditions

#Example (shape file):

#       files               def
# 1     filename1.mat        Basline
# 2     filename2.mat        TR_X
# 3     filename3.mat        TR_Y
# 4     filena44e.mat        TR_Z



# load MEA and shape file (experimental shape) data into a MEA object

object=MEA(path, def, TR="Post_Glutatate")
object@Shape=shape

#Define a dynamic basline based on your "Basline" input
object=Align_Zero(object,STD=10,Edges=1000,Exp_start=100000,Exp_end=500000, Exp_chanel=21,plot=T)

object=Finde_Peaks(object)
object=Find_IPI(object,samp_freq=1000)

#optional
saveRDS(object, "object_233.R")

#Layout file defines your MEA layout
layout=read.csv("Layout.csv", sep=";")

#Analysis and plotting of the raster- or a scatterplot
Plot_Events_Scatter(object,RasterP=T,layout)

```


### Authors

D. H. Heiland & K. Joseph, Translational Research Group, Medcal-Center Freiburg, University of Freiburg

