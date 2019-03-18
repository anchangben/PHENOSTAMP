# STAMP

This repository contains the code for analyzing single cell data using STATE MAP (STAMP) algorithm (Karacosta et al. under review). This algorithm is an extension of the Clustering Classification and sorting Tree(CCAST) algorithm (Anchang et al. 2014) with an optimal 2D visualization of projected single-cell data using a neural network.

It contains a dataset representing the protein expression levels for cells undergoing an Epithelial Mesenchymal Transition (EMT), a biological process exploited by cancer cells for invasion and metastasis. This dataset is separated in two parts as follows :
+ newdat.rdata contains 6 markers selected for their relevance in the EMT process
+ newdat2.rdata contains the whole raw dataset

The code needed to generate the results in the paper are accessible using the main.r script provided.

## Installation

A complete installation tool is available. All you have to do is run :
```{R}
source("setup.R")
followed by 
source("main.R")
```

The complete installation script can be long to run (between half an hour and an hour), but is complete.
