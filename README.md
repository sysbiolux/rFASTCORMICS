# INTRODUCTION

The FASTCORE algorithm family is a collection of model-building algorithms that allow for the reconsctruction of context-specific metabolic models based on a generic genome-scale metabolic reconstruction and some input data.

For more details on the FASTCORE algorithm family see:

	- FASTCORE (Vlassis et al., 2014, PloS Computational Biology)
	- FASTCORMICS (Pacheco et al., 2015, BMC Genomics)
	- Benchmarking (Pacheco et al., 2016, Frontiers in Physiology)
	- rFASTCORMICS (Pacheco et al., 2019, EBioMedicine)



#  for rFASTCORMICS

Similar to FASTCORMICS but takes RNA-seq data as input, preferably FPKM transformed.

## PREREQUISITES

SOFTWARE

	- Matlab 2013 or higher
		- compatible IBM CPLEX installation
		- Statistics and Machine Learning Toolbox
		- Curve Fitting Toolbox 
		- COBRA Toolbox (https://opencobra.github.io/cobratoolbox/latest/installation.html )
		
DATA (provided in the [exampleData folder](https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/exampleData)

	- RNA-seq data (FPKM transformed)
		- colnames:		1xC cell with the sample names (size C)
		- rownames:		Rx1 cell with the gene identifiers (size R)
		- fpkm: 		RxC double matrix or table containing the fpkm values
	- model:		(consistent) genome-scale metabolic reconstruction in the COBRA format, i.e. Recon 2.04 (from https://vmh.uni.lu/#downloadview )
	- dico:			table which contains corresponding gene identifier information. Needed to map the rownames to the genes in the model. Can be manually assembled in https://www.ensembl.org/biomart/martview and imported into Matlab.
	- medium:		[optional] defines metabolites in the growth medium of cells to constrain the model, see example medium_example.mat
	
	
## USAGE

An [example script](https://github.com/sysbiolux/rFASTCORMICS/blob/master/rFASTCORMICS%20for%20RNA-seq%20data/fastcormicsRNAseq_example_v4.mlx) for the creation of context-specific models based on two samples from the TCGA data is provided ([from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944)).
Note that the FPKM values have been used.

The two included samples are TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07.
In the provided example, we will create single models for each sample, as well as one consensus model from both samples. 
For the latter, a consensus proportion needs to be decided. 
The default consensus is 0.9, meaning a reaction is only considered present in the final model, if it can be derived from at least 90% from the input data.

An [additional script](https://github.com/sysbiolux/rFASTCORMICS/blob/master/rFASTCORMICS%20for%20RNA-seq%20data/book_example_4_publish_v2.mlx) is also provided that contains more samples as well as some visualization methods and the drug target prediction workflow.

## Troubleshooting

For more information and possible troubleshooting, see the [original publication](https://doi.org/10.1016/j.ebiom.2019.04.046) and below.


> "The input was too complicated or too big for MATLAB to parse"

For some newer versions of Matlab you might encounter this error. Please use *feature astheightlimit 2000* in the beginning of your script.


#  for FASTCORMICS

## PREREQUISITES

SOFTWARE

	- Matlab 2013 or higher
		- compatible IBM CPLEX installation
		- COBRA Toolbox (https://opencobra.github.io/cobratoolbox/latest/installation.html )
	- R (optional: R Studio) with 
		- BiocManager
		- affy
		- frma
		- corresponding BARCODE vectors
			- hugene.1.0.st.v1frmavecs or hgu133afrmavecs or hgu133plus2frmavecs or hgu133a2frmavecs
		
DATA (provided in the [exampleData folder](https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/exampleData)

	- raw microarray data (.CEL files)
		- colnames:		1xC cell with the sample names (size C)
		- rownames:		Rx1 cell with the gene identifiers (size R)
		- barcode: 		RxC double matrix or table containing the barcode-transformed values
	- model:		(consistent) genome-scale metabolic reconstruction in the COBRA format, i.e. Recon 2.04 (from https://vmh.uni.lu/#downloadview )
	- dico:			table which contains corresponding gene identifier information. Needed to map the rownames to the genes in the model. Can be manually assembled in https://www.ensembl.org/biomart/martview and imported into Matlab.
	- medium:		[optional] defines metabolites in the growth medium of cells to constrain the model, see example medium_example.mat
	
	
## USAGE

The script for FASTCORMICS is divided into two parts:

### Barcode-transformation of microarray data in R

In the provided [example script](https://github.com/sysbiolux/rFASTCORMICS/blob/master/FASTCORMICS%20for%20microarray%20data/Barcode%20for%20FASTCORMICS.R), the .CEL files from a microarray experiment are read using the affy package. Then the data is fRMA-normalized followed by the BARCODE transformation of retrieving z-scores. BARCODE 3 has been [published](https://doi.org/10.1093/nar/gkt1204).

### Discretiation and model building in Matlab

The z-scores obtained from BARCODE are read into Matlab and discretized as follows

## Troubleshooting

> my models reconstructed with FASTCORMICS are very small

Please check the number of expressed genes after the discretization. If the number is low < 30 %, you can try changing the expression threshold to 3 instead of 5.


#  for FASTCORE
under construction

# ABOUT
	

## rFASTCORMICS

FASTCORMICS RNA-seq (c) was published in "Identifying and targeting cancer-specific metabolism with network-based drug target prediction".

Pacheco, M. P., Bintener, T., Ternes, D., Kulms, D., Haan, S., Letellier, E., & Sauter, T. (2019). Identifying and targeting cancer-specific metabolism with network-based drug target prediction. EBioMedicine, 43, 98-106.

https://doi.org/10.1016/j.ebiom.2019.04.046
https://www.sciencedirect.com/science/article/pii/S2352396419302853


## FASTCORMICS

FASTCORMICS (c) was published in "Integrated metabolic modelling reveals cell-type specific epigenetic control points of the macrophage metabolic network".

Pacheco, M.P., John, E., Kaoma, T., Heinäniemi, M., Nicot, N., Vallar, L., Bueb, J.L., Sinkkonen, L. and Sauter, T. (2015). Integrated metabolic modelling reveals cell-type specific epigenetic control points of the macrophage metabolic network. BMC genomics, 16(1), 1-24.

https://doi.org/10.1186/s12864-015-1984-4
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1984-4

## FASTCORE

FASTCORE (c) was published in "Fast reconstruction of compact context-specific metabolic network models".

Vlassis, N., Pacheco, M. P., & Sauter, T. (2014). Fast reconstruction of compact context-specific metabolic network models. PLoS Comput Biol, 10(1), e1003424.

https://doi.org/10.1371/journal.pcbi.1003424
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424	
	
	
	
	
Tamara Bintener

