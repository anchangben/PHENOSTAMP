
source("ccast_compute.R")

## main

libs <- c("spade", "graph", "fastcluster", "party",  
          "plot3D", "MASS", "RColorBrewer", "flowCore", "cluster",
          "plotly", "bigvis", "tripack", "deldir", "sp", "Rtsne")

lapply(libs, library, character.only = TRUE)

outdir <- "./output/"
dir.create(outdir)
dir.create(paste(outdir,"First_pass",sep="/"))
dir.create(paste(outdir,"Second_pass",sep="/"))

load("data/newdat.rdata")
load("data/newdat2.rdata")
dataset = asinh(newdat)

print("Starting CCAST tree clustering starting with a hierarchical clustering of the data")
print("This can take several minutes...")
first_res <- ccast_tree_twopass(dataset, outputDir=paste(outdir,"First_pass",sep="/"),
                            k=8, groups = NULL, ylabel="Vimentin", s.plot=TRUE, s.save=FALSE,
                            max_depth=5, min_bucket=1000, subanalysis=FALSE)

print("Fist CCAST tree done")
final_dataset <- first_res$Finaldata[,-dim(first_res$Finaldata)[2]]

finalclusters <- ccast_tree_analyze_and_clean(tree=first_res$finaltree, DD=as.data.frame(final_dataset), s.save=TRUE, outputDir=paste(outdir,"Second_pass",sep="/"))

first_pass_groups <- finalclusters$newclusters

final_res <- ccast_tree_twopass(as.data.frame(final_dataset),outputDir=paste(outdir,"Second_pass",sep="/"),k=length(unique(first_pass_groups)),
                                ylabel="Vimentin", groups=first_pass_groups, s.save=FALSE, s.plot=TRUE, subanalysis=FALSE,max_depth=5, min_bucket=5)

groups <- Predict(final_res$finaltree, as.data.frame(dataset))

Timepoint <- newdat2[, "Timepoint"]
print("CCAST finalized")

set.seed(123)
print("Starting TSNE computation")
print("This can take up to an hour...")
tsnedat <- Rtsne(asinh(newdat),perplexity = 30)$Y

optimalbinsize <- optimalbinsize2(groups, Timepoint, as.data.frame(tsnedat), 15)
vor <- voronoi_mapping(groups, tsnedat,optimalbinsize)

antibody1=colnames(newdat2)

rr=trainmodel(tsnedat,newdat,trainprop=0.9,hidden=11,maxit=2000)
nn1=rr[[1]]

window = c(min(tsnedat[, 1])-5, max(tsnedat[, 1])+5, min(tsnedat[, 2])-5, max(tsnedat[, 2])+5)
maptsne(tsnedat,asinhp=NULL,newdat,newdat2,antibody1,nn=nn1,outputDir='./output', vor=vor,file1="EMT timecourse data.tiff",sample="EMT timecourse data", window=window)
for (time in unique(Timepoint)){
  maptsne(tsnedat[newdat2[, "Timepoint"] == time, ],
          asinhp=NULL,newdat[newdat2[, "Timepoint"] == time, ],
          newdat2[newdat2[, "Timepoint"] == time, ],NULL,
          nn=nn1,outputDir='./output', vor=vor,file1=paste("EMT_timecourse_", time, ".tiff"),
          sample="EMT timecourse data", window=window)
}