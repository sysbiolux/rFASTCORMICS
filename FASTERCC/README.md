# FASTERCC

⚡ **FASTERCC** is an accelerated implementation of flux consistency testing designed for large-scale metabolic network models, including single-cell–derived reconstructions.

---

## PUBLICATION

**FASTERCC: Accelerating Flux Consistency Testing and Context-Specific Reconstruction for Large-Scale Metabolic Network Models**  
Maria Pacheco, Evelyn Gonzalez, Thomas Sauter  
DOI: https://doi.org/10.64898/2026.03.19.712885  

### Status
This work is currently available as a **preprint**.

---

## ABSTRACT

The increase in size of metabolic network models especially with the advent of single-cell data calls for scalable reconstruction and analysis tools. Such models, often used for drug discovery and the analysis of microbial communities rely on consistency testing and reconstruction algorithms such as FASTCORE and FASTCC. However, with models nowadays comprising hundreds of thousands of reactions, the running times of such algorithms increased from few minutes to hours or days even with high performance computing. Experiments that require multiple reconstructions, such as parameter tuning or cross-validation, are practically infeasible in very large networks. 
Here we introduce a new version of FASTCC that leverages structural information for removing type I and II dead-ends, the orientation of reversible reactions and correcting the reversibility of reactions that are structurally incapable of carrying flux in both directions prior to any feasibility tests. These improvements reduce drastically the running time of FASTERCC by a median 20-fold speedup in comparison to FASTCC for networks with a larger number of block reactions. The model cleaning performed by FASTERCC also reduces the computational time of downstream analyses, notably of FASTCORE up to 50%.


**Key results:**

- Median **20× speedup** compared to FASTCC  
- Up to **50% reduction in FASTCORE runtime**

---

## INTRODUCTION

FASTERCC extends and optimizes the FASTCORE family by focusing on efficient flux consistency testing in large and highly constrained genome-scale metabolic models (GEMs). The core function, `fasterCC`, integrates structural preprocessing with LP-based consistency testing to identify and remove blocked reactions.

Applying physiological constraints—such as medium composition or gene expression filtering—often disrupts network connectivity, leading to flux-inconsistent models. FASTERCC addresses this by simplifying the model prior to and during optimization. It corrects reaction directionality, removes dead ends, reduces structurally constrained reversible reactions, and iteratively refines the network.

An optional reversibility testing step further improves model quality by reorienting reactions and tightening bounds for reactions that operate only in one direction. This results in a cleaner and more standardized model, significantly improving the performance of downstream algorithms such as FASTCORE.

---

## METHOD OVERVIEW

FASTERCC improves flux consistency testing by combining structural preprocessing with efficient optimization.

The workflow includes:

- Removal of type I and II dead-end metabolites  
- Correction of reaction directionality  
- Conversion of structurally constrained reversible reactions into irreversible ones  
- Iterative flux consistency testing (batch and one-by-one modes)  
- Optional post-processing to refine reaction reversibility  

These steps reduce the number of blocked and reversible reactions, minimizing costly feasibility checks and improving scalability for large models.

---

## PREREQUISITES

### Software

- MATLAB (2013 or higher)  
- COBRA Toolbox  
  https://opencobra.github.io/cobratoolbox/latest/installation.html  
- Compatible LP solver (e.g., IBM ILOG CPLEX)

---

### Required Functions

Ensure the following functions are in your MATLAB path:

- `fasterCC.m`
- `fixIrr_rFASTCORMICS.m`
- `restoreBounds.m`
- `findDeadEndsFastbox.m`
- `structureAnalyseFastbox.m`
- `LP3_rFASTCORMICS.m`
- `LP7_rFASTCORMICS.m`
- `LP9_rFASTCORMICS.m`

---

### COBRA Toolbox Functions Used

- `removeRxns`

---

### Data

- A COBRA-format metabolic model  
- Example large-scale model: `Recon3D.mat`

---

## USAGE

Basic usage in MATLAB:

```matlab
[A, model_consistent] = fasterCC(model, epsilon, printLevel, revopt);
```

---

## DEMO

A full working example is provided in the `code` folder:

```matlab
code/run_fasterCC_Recon3D_demo.m
```

---
