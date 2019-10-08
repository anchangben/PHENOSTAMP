# PHENOSTAMP

This repository contains the code for analyzing single cell data using PHENOtypic STAte MaP (PHENOSTAMP) algorithm (Karacosta et al. accepted by Nature Communication 2019). This algorithm combines the Clustering Classification and sorting Tree(CCAST) algorithm (Anchang et al. 2014) with an optimal 2D visualization of projected single-cell data using a neural network.

It contains a dataset representing the protein expression levels for cells undergoing an Epithelial Mesenchymal Transition (EMT), a biological process exploited by cancer cells for invasion and metastasis. This dataset is separated in two parts as follows :
+ newdat.rdata contains 6 markers selected for their relevance in the EMT process
+ newdat2.rdata contains the whole raw dataset with many more markers

The codes needed to generate the results from new unprocessed data are accessible using the main.r script provided.

## Installation
PHENOSTAMP relies on R (>= 3.5.0) and other R libraries (See main.R). The current implementation below runs on Mac OS X 10.7 and above operating systems. Instructions to setup the macOS toolchain for compiling used in the 3.5.* and above series of R  can be found in  https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/. A complete installation of PHENOSTAMP requires that all dependencies are installed.

A complete installation tool is available. All you have to do is run :
```{R}
source("setup.R")
followed by 
source("main.R") to run the complete analysis.
```
The complete installation script can be long to run. 

If you have a new EMT data already processed as a data matrix r object and you want to project on the map from Karacosta et al. 2019, use code from "project.R" instead.

## Questions
Any questions can be addressed to anb28636@gmail.com
