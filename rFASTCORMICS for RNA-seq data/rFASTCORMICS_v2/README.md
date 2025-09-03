## What is rFASTCORMICS_v2?
rFASTCORMICS is an algorithm for the reconstruction of context-specific metabolic models based on RNA-seq data that has first been published in 2019 (https://doi.org/10.1016/j.ebiom.2019.04.046). 
Here, we provide an updated version of rFASTCORMICS, called rFASTCORMICS_v2.

## What will you find on this repository?
On this repository, you will find:
- A `README.md` file, containing a description of this repository, and some tutorials.
- A *scripts* folder, containing the scripts of rFASTCORMICS_v2, and the ones of FASTCORE. It is important to use the functions from this FASTCORE folder instead of the one included in the COBRA Toolbox, as some changes were also made on these functions and that they are called in rFASTCORMICS_v2. The two folders can directly be associated with to the COBRA Toolbox.
- A *BCRA_example* folder, containing materials and a script to run an example using rFASTCORMICS_v2. The `test_script.m` file can be run directly. It will load the `workspace_example` and the `optionalSettings` files. The `finalWorkspace` file contains all the variables after the script has finished running.

## What are the main changes between rFASTCORMICS and rFASTCORMICS_v2?
1. Instead of not penalizing the addition of transporters in the model, we removed the nonPen argument and rather decided to remove transporters from the core set of reactions. In this way, we make sure only the necessary transporters will be included in the model.
2. The FASTCORE algorithm is run only once instead of twice, reducing the time of execution of rFASTCORMICS_v2.

## Installation
rFASTCORMICS_v2 requires MATLAB 2019 or a later version and works with the IBM CPLEX or Gurobi solvers.<br> The user also needs to install the COBRA Toolbox. As mentioned previously, both FASTCORE and rFASTCORMICS_v2 folders in the *scripts* folders sould be downloaded. The FASTCORE folder already included in the COBRA Toolbox should be replaced by the one provided here.

## How to use rFastcormics_v2?
Here is the workflow of rFastcormics_v2:

<img src="WorkflowrFastcormics_v2.png" width="550">

**Usage**

`[contextSpecificModel, retainedRxns] = rFastcormics4cobra_v2(model, discretized, rownames, dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, fillingMediumFlag, adaptiveScalingFlag)`

| **Inputs** | **Description** |
|--------------|-----------------|
| `model`      | metabolic model (following fields are required: S, lb, ub, rxns, metFormulas) |
| `discretized` | matrix with the discretized value (1, 0, or -1) for each gene and sample (as many rows as the number of genes, as many columns as the number of samples) |
| `rownames` | cell array with the gene IDs (as many rows as the number of genes in `discretized`) |
| `dico` | table with two columns, first for gene IDs in the `model` format, second for gene IDs in the `discretized` format |
| **Optional inputs** | 
| `consensusProportion` | the rate of samples that have to express or not to express a gene for the gene to be considered expressed or not (default 0.9) |
| `epsilon` | smallest flux that is considered nonzero (default getCobraSolverParams('LP', 'feasTol')*100 = 1e-4) |
| `optionalSettings.func` | cell array of reaction abbreviations that should carry a flux <br> It is recommended to put the objective function of your model to ensure its preservation in the context-specific model |
| `optionalSettings.medium` | cell array of metabolite abbreviations that are present in the growth medium of the cells and that will be used to constrain the model |
| `optionalSettings.notMediumConstrained` | ??? |
| `biomassReactionName` | ??? |
| `fillingMediumFlag` | fill the medium with supplementary reactions in case the provided medium is not sufficient to fulfill the objective function<br> 1 for active (default), 0 for inactive |
| `adaptiveScalingFlag` | adaptive scaling of the flux values (see LP10)<br>0 for inactive (default), 1 for active |
| **Outputs** |
| contextSpecificModel | context-specific model, reduced to the retained reactions and associated genes |
| Afinal | indices in `model` of the retained reactions |
