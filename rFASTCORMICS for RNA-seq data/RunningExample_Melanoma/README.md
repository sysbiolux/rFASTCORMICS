# Simple running example rFASTCORMICS
## Metabolic network modelling of A375 melanoma and Mold melanocyte (control) cell lines

Paper: [Bintener et al, 2021](https://pubmed.ncbi.nlm.nih.gov/37495601/)

Data: [EGAS00001006463](https://ega-archive.org/datasets/EGAD00001009089)
The consensus models can be downloaded from [consensus models](https://github.com/sysbiolux/MelanomaPaper).
The FPKM can be obtained from the authors on request.

In this example we reconstruct Metabolic Networks Models of A375 and Mold consensus models, as well as a sample specific model for A375 sample1.
This will demonstrate the model building and analysis pipeline with short run time (approx. 20 minutes on local desktop).

driverModBuild_melanoma: Performs some quality control and PCA of the RNAseq data, generates medium constraint Metabolic Network Models and checks for biomass production of the sample specific A375 model.

driverAnalysis_melanoma: Generates some model stats and filters specific metabolites and genes, identifies pathways with different reaction presence in the consensus models
and performs drug screening on the A375 consensus model.
