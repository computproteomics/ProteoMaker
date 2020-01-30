# PhosFake

This repository contains the scripts necessary for the generation of an in-silico bottom-up phosphoproteomics data set. 

All the parameters that are used to generate the data are listed in `Parameter.R`. The script that runs the entire pipeline is `RunAll.R`.

The pipeline is described in the figure `PhosFakeLayout.svg` and can be described as follows:

1) Generation of ground truth data at the proteoform level (`01_GenerateGroundTruth.R`):
