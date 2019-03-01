# CCAST

This repository contains the code for analyzing single cell data with the CCAST algorithm.

It contains a dataset representing the protein expression levels for cells undergoing an EMT transition process. This dataset is separated in two parts as follows :
+ newdat.rdata contains 6 markers selected for their relevance in the EMT process
+ newdat2.rdata contains the whole raw dataset

The code needed to generate the results in the paper are accessible using the main.r script provided. 

## Installation

A complete installation tool is available. All you have to do is run :
```{R}
source("setup.R")
```

The complete installation script can be long to run (between half an hour and an hour), but is complete.
