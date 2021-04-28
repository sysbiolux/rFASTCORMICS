# INTRODUCTION

The FASTCORE algorithm family is a collection of model-building algorithms that allow for the reconsctruction of context-specific metabolic models based on a generic genome-scale metabolic reconstruction and some input data.

For more details on the FASTCORE algorithm family see:

	- FASTCORE (Vlassis et al., 2014, PloS Computational Biology)
	- FASTCORMICS (Pacheco et al., 2015, BMC Genomics)
	- Benchmarking (Pacheco et al., 2016, Frontiers in Physiology)
	- rFASTCORMICS (Pacheco et al., 2019, EBioMedicine)

Last major update was on April 26, 2021.

<!------------------------------------------------------------------------------------------------>
#  for rFASTCORMICS

rFASTCORMICS is an algorithm for the reconstruction of context-specific metabolic models based on RNA-seq data. It is similar to FASTCORMICS for microarray data but takes RNA-seq data as input, preferably FPKM transformed.

## PREREQUISITES

SOFTWARE

	- Matlab 2013 or higher
		- compatible IBM ILOG CPLEX installation
		- Statistics and Machine Learning Toolbox
		- Curve Fitting Toolbox 
		- COBRA Toolbox (https://opencobra.github.io/cobratoolbox/latest/installation.html )
		
DATA (provided in the [exampleData folder](https://github.com/sysbiolux/rFASTCORMICS/tree/master/rFASTCORMICS%20for%20RNA-seq%20data/exampleData)

	- RNA-seq data (FPKM or TPM transformed)
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

## INPUTS

The code to run rFASTCORMICS looks as follows:

``` Matlab
[model, A_final] = fastcormics_RNAseq(model, data, rownames, dico, biomass_rxn, ...
			already_mapped_tag, consensus_proportion, epsilon, optional_settings)
```


**Required inputs:**
Input | Explanation
------------ | -------------
model | (consistent) genome-scale metabolic reconstruction in the COBRA format, i.e. [Recon 2.04](https://vmh.uni.lu/#downloadview)
data | discretized experimental data (1 for expressed, -1 for not expressed, and 0 for unknown expression genes)
rownames | cell aray with the gene IDs from the experiment
dico | table that contains corresponding gene identifier information. Can also be a matrix. Needed to map the rownames to the genes in the model. Can be manually assembled in [Biomart](https://www.ensembl.org/biomart/martview) and imported into Matlab as a text-only table.


**Optional inputs:**
Input | Explanation | Default
------------ | -------------| ------------- 
biomass_rxn | name of the biomass reaction in the model if present (check model.rxns). Setting this variable will always enable the biomass to carry a flux. Some examples are: biomass_reaction, biomass_components, biomass_human,... | ''
already_mapped_tag | 1, if the data was already to the model.rxns in this case data p = n  and  0, if the data has to be mapped  using the GPR rules of the model	|	0
consensus_proportion |	gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples |	0.9
epsilon | to avoid small number errors	|	1e-4
optional_settings | a structure with the following variables:	|	''
/ | unpenalizedSystems: 	|	
/ | medium: medium composition, defines metabolites in the growth medium of cells to constrain the model |	
/ | not_medium_constrained: 	|	
/ | func: reaction(s) forced to be present in the model	|	





<!------------------------------------------------------------------------------------------------>
#  for FASTCORMICS

FASTCORMICS was designed to reconstruct context-specific metabolic models based on microarray data.

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

For R, the [example script](https://github.com/sysbiolux/rFASTCORMICS/blob/master/FASTCORMICS%20for%20microarray%20data/Barcode%20for%20FASTCORMICS.R) will automatically download and install the required packages as well as perform the BARCODE transformation of the data. 
This output of the examples script will be saved in 4 separate files in your working directory and will be used for the reconstruction of context specific models in the ['FASTCORMICS example script](https://github.com/sysbiolux/rFASTCORMICS/blob/master/FASTCORMICS%20for%20microarray%20data/FASTCORMICS_Example.mlx):

	- barcode.txt		BARCODE transofrmed data, based on the .CEL input
	- colnames.txt		sample names of the .CEL input
	- frma.txt		fRMA-normalized data
	- rownames.txt		probe IDs from the .CEL files
		
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

When using different data, make sure to use the correct BARCODE vectors for your platform/GeneChip:

GeneChip | Platform | BARCODE vector
------------ | ------------- | -------------
U133A	| GPL96	| hgu133a2frmavecs
U133 plus 2.0	| GPL570	| hgu133plus2frmavecs
U133A 2.0	| GPL571	| hgu133afrmavecs
Human Gene 1.0 ST	| GPL6244	| hugene.1.0.st.v1frmavecs


### Discretiation and model building in Matlab

The z-scores obtained from BARCODE are read into Matlab and discretized as follows

	- not expressed: 	z-score <= 0
	- unknown expression: 	0 < z-score < 5
	- expressed: 		5 < z-score

The discretized values are then used to find active reactions in the model based on the GPR-rules, see [original publication](https://doi.org/10.1186/s12864-015-1984-4) for a full explanation.




<!------------------------------------------------------------------------------------------------>
#  for FASTCORE

FASTCORE can be used to reconstruct a context-specific metabolic model based on a list of reactions, called core reactions, that are known to take place in the context of interest.

In the [short provided example](https://github.com/sysbiolux/rFASTCORMICS/blob/master/FASTCORE/FASTCORE_example.mlx), we use Recon 1 and a previously compiled list of liver reactions to reconstruct a liver core model.


<!------------------------------------------------------------------------------------------------>

# Troubleshooting

For more information and possible troubleshooting, see the [original publication](https://doi.org/10.1016/j.ebiom.2019.04.046) and below.

> My FPKM density plots look very different and are not processed correctly during the discretization step, i.e. the peaks are not correctly determined.

If you observe 2 peaks in your FPKM density plot and the left peak is higher than the rightmost peak, you can try to use the discretize_FPKM_skewed function instead. Otherwise, make sure that no altering pre-filtering steps of unexpressed genes has been performed.

> "The input was too complicated or too big for MATLAB to parse"

For some newer versions of Matlab you might encounter this error. Please use *feature astheightlimit 2000* in the beginning of your script.

> I do not have a biomass reaction in my model. What can I do?

You can use [], to omit any input for fastcormics_RNAseq, or define biomass_rxns = {''}.

> I get the warning that optional settings are not set even though I did. I also get errors.

Please note that the inputs of fastcormics_RNAseq have changed recently. Please re-check the inputs, you might have forgotten to define the biomass_rxns input.

> My models reconstructed with FASTCORMICS are very small

Please check the number of expressed genes after the discretization. If the number is low < 20 %, you can try changing the expression threshold to 3 instead of 5.



<!------------------------------------------------------------------------------------------------>


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
	
