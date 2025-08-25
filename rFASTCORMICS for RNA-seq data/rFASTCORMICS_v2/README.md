## What is rFASTCORMICS_v2?
rFASTCORMICS is an algorithm for the reconstruction of context-specific metabolic models based on RNA-seq data that has first been published in 2019 (https://doi.org/10.1016/j.ebiom.2019.04.046). 
Here, we provide an updated version of rFASTCORMICS, called rFASTCORMICS_v2.

## What will you find on this repository?
On this repository, you will find:
- A README.md file, containing a description of this repository, and some tutorials.
- A *scripts* folder, containing the scripts of rFASTCORMICS_v2, and the ones of FASTCORE. It is important to use the functions from this FASTCORE folder instead of the one included in the COBRA Toolbox, as some changes were also made on these functions and that they are called in rFASTCORMICS_v2. The two folders can directly be joined to the COBRA Toolbox.
- A *BCRA example* folder, containing material to run an example using rFASTCORMICS_v2.

## What are the main changes between rFASTCORMICS and rFASTCORMICS_v2?
1. Instead of not penalizing the addition of transporters in the model, we removed the nonPen argument and rather decided to remove transporters from the core set of reactions. In this way, we make sure only the necessary transporters will be included in the model.
2. The FASTCORE algorithm is run only once instead of twice, reducing the time of execution of rFASTCORMICS_v2.

## How to use rFastcormics_v2?
Here is the workflow of rFastcormics_v2:
<img src="WorkflowrFastcormics_v2.png">
