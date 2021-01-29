# rFASTCORMICS


FASTCORMICS RNA-seq (c) Pacheco et al., 2018 02/07/2018

INTRODUCTION
==========================

FASTCORMICS RNA-seq is an automated workflow to create a context specific model from the input data and a genome scale reconstruction.

For more details on the FASTCORE algorithm family see:
	- FASTCORE (Vlassis et al., 2014, PloS Computational Biology)
	- FASTCORMICS (Pacheco et al., 2015, BMC Genomics)
	- Benchmarking (Pacheco et al., 2016, Frontiers in Physiology)

PREREQUISITES
==========================

SOFTWARE
	- Matlab 2013 or higher
		- compatible IBM CPLEX installation, (please use this link to check which version to use https://www.ibm.com/software/reports/compatibility/clarity/softwarePrereqsMatrix.html )
		- Statistics and Machine Learning Toolbox
		- Curve Fitting Toolbox 
		- COBRA Toolbox (https://opencobra.github.io/cobratoolbox/latest/installation.html )
		
DATA
	- RNA-seq data (FPKM transformed)
		- colnames:		1xC cell with the sample names (size C)
		- rownames:		Rx1 cell with the gene identifiers (size R)
		- fpkm: 		RxC double matrix containing the fpkm values
	- model:		(consistent) genome-scale metabolic reconstruction in the COBRA format, i.e. Recon 2.04 (from https://vmh.uni.lu/#downloadview )
	- dico.mat		table which contains corresponding gene identifier information. Needed to map the rownames to the genes in the model
	- medium.mat:	[optional] defines metabolites in the growth medium of cells to constrain the model, see example medium_example.mat
	
	
USAGE
==========================
An example for the creation of context-specific models based on two samples from the TCGA data is provided (from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944 ).
Note that, for the publication, the FPKM values have been used.

The two included samples are TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07.
In the provided example, we will create single models for each sample, as well as one generic model from both samples. 
For the latter, a consensus proportion needs to be decided. 
The default consensus is 0.9, meaning a reaction is only considered present in the final model, if it can be derived from at least 90% from the input data.


SCRIPT	
==========================
First, the RNA-seq data is loaded.
The figure flag (figflag) can be set to 1 (show figures) or 0(do not show figures). If set to 1, the density plot with hybrid curves will be saved for each sample in the Figures folder.
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
	
	
	
Tamara Bintener

