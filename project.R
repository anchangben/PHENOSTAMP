
source("stamp_compute.R")

## main

libs <- c( "graph", "fastcluster", "party",  
          "plot3D", "MASS", "RColorBrewer", "flowCore", "cluster",
          "plotly", "bigvis", "tripack", "deldir", "sp", "Rtsne")

lapply(libs, library, character.only = TRUE)

outdir <- "./output/"
dir.create(outdir)

load("data/newdat_8clusters.rdata")
load("data/newdat2_8clusters.rdata")
dataset = asinh(newdat)

##load tsne output (tsnedat)
##load group or cluster labels (groups)
##load voronoi boundaries (vor)
## load neural network (nn1)
load("data/tsnedat.rdata")
load("data/groups.rdata")
load("data/vor.rdata")
load("data/nn1.rdata")

Timepoint <- newdat2[, "Timepoint"]
antibody1=colnames(newdat2)

### plot EMT map for complete data ##
window = c(min(tsnedat[, 1])-5, max(tsnedat[, 1])+5, min(tsnedat[, 2])-5, max(tsnedat[, 2])+5)
maptsne(tsnedat,asinhp=NULL,newdat=dataset,newdat2,antibody1,nn=nn1,outputDir='./output', vor=vor,file1="EMT timecourse data.tiff",sample="EMT timecourse data", window=window)

### plot EMT map for each time or condition ##
for (time in unique(Timepoint)){
  maptsne(tsnedat,
          asinhp=NULL,newdat[newdat2[, "Timepoint"] == time, ],
          newdat2[newdat2[, "Timepoint"] == time, ],NULL,
          nn=nn1,outputDir='./output', vor=vor,file1=paste("EMT_timecourse_", time, ".tiff"),
          sample="EMT timecourse data", window=window)
}
