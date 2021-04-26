INTRODUCTION
==========================

The FASTCORE algorithm family is a collection of model-building algorithms that allow for the reconsctruction of context-specific metabolic models based on a generic genome-scale metabolic reconstruction and some input data.

For more details on the FASTCORE algorithm family see:

	- FASTCORE (Vlassis et al., 2014, PloS Computational Biology)
	- FASTCORMICS (Pacheco et al., 2015, BMC Genomics)
	- Benchmarking (Pacheco et al., 2016, Frontiers in Physiology)
	- rFASTCORMICS (Pacheco et al., 2019, EBioMedicine)

PREREQUISITES
==========================

##  for rFASTCORMICS

SOFTWARE

	- Matlab 2013 or higher
		- compatible IBM CPLEX installation
		- Statistics and Machine Learning Toolbox
		- Curve Fitting Toolbox 
		- COBRA Toolbox (https://opencobra.github.io/cobratoolbox/latest/installation.html )
		
DATA

	- RNA-seq data (FPKM transformed)
		- colnames:		1xC cell with the sample names (size C)
		- rownames:		Rx1 cell with the gene identifiers (size R)
		- fpkm: 		RxC double matrix or table containing the fpkm values
	- model:		(consistent) genome-scale metabolic reconstruction in the COBRA format, i.e. Recon 2.04 (from https://vmh.uni.lu/#downloadview )
	- dico:			table which contains corresponding gene identifier information. Needed to map the rownames to the genes in the model. Can be manually assembled in https://www.ensembl.org/biomart/martview and imported into Matlab.
	- medium:		[optional] defines metabolites in the growth medium of cells to constrain the model, see example medium_example.mat
	
	
USAGE
==========================
An example script for the creation of context-specific models based on two samples from the TCGA data is provided (from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944 ).
Note that, for the publication, the FPKM values have been used.

The two included samples are TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07.
In the provided example, we will create single models for each sample, as well as one consensus model from both samples. 
For the latter, a consensus proportion needs to be decided. 
The default consensus is 0.9, meaning a reaction is only considered present in the final model, if it can be derived from at least 90% from the input data.


SCRIPT	fastcormicsRNAseq_example_v4
==========================
First, the RNA-seq data is loaded.
Then, the fpkm data is log2-transformed and discretized (for details, see Material and Methods from Pacheco et al., 2018).

Finally, the context-specific models are reconstructed taking into consideration the optional_settings such as 

	- unpenalized systems
	- function to be forced into the model (biomass_reraction and DM_atp_c)
	- not_medium_constrained
	- medium constraints
	
and the mandatory settings:

	- epsilon
	- consensus proportion
	- dico


##  for rFASTCORMICS
under construction

##  for FASTCORE
under construction

ABOUT
==========================	

# rFASTCORMICS

FASTCORMICS RNA-seq (c) was published in "Identifying and targeting cancer-specific metabolism with network-based drug target prediction".

Pacheco, M. P., Bintener, T., Ternes, D., Kulms, D., Haan, S., Letellier, E., & Sauter, T. (2019). Identifying and targeting cancer-specific metabolism with network-based drug target prediction. EBioMedicine, 43, 98-106.

https://doi.org/10.1016/j.ebiom.2019.04.046
https://www.sciencedirect.com/science/article/pii/S2352396419302853


# FASTCORMICS

FASTCORMICS (c) was published in "Integrated metabolic modelling reveals cell-type specific epigenetic control points of the macrophage metabolic network".

Pacheco, M.P., John, E., Kaoma, T., Heinäniemi, M., Nicot, N., Vallar, L., Bueb, J.L., Sinkkonen, L. and Sauter, T. (2015). Integrated metabolic modelling reveals cell-type specific epigenetic control points of the macrophage metabolic network. BMC genomics, 16(1), 1-24.

https://doi.org/10.1186/s12864-015-1984-4
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1984-4

# FASTCORE

FASTCORE (c) was published in "Fast reconstruction of compact context-specific metabolic network models".

Vlassis, N., Pacheco, M. P., & Sauter, T. (2014). Fast reconstruction of compact context-specific metabolic network models. PLoS Comput Biol, 10(1), e1003424.

https://doi.org/10.1371/journal.pcbi.1003424
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424	
	
	
	
	
Tamara Bintener

