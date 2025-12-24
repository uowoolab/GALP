# Machine Learned Probability Distributions

This document describes the optional machine learned probability distribution workflow implemented in GALP. This method was developed by Jake Burner and is described in detail in the following preprint:

[J. Burner *et al.*, 2025](https://chemrxiv.org/engage/chemrxiv/article-details/68bf29fd23be8e43d6c7af1f)

The purpose of this workflow is to bypass explicit GCMC simulations by using a trained machine learning model to predict adsorbate probability distributions at specified temperature and pressure conditions. The predicted probability grids are compatible with GALP and can be used to identify binding sites in the same manner as grids generated from GCMC simulations.


## Method Overview

In this workflow, binding site maxima extracted with GALP from GCMC generated probability grids were used during training to identify relevant regions of configuration space. The resulting model predicts complete adsorbate probability distributions directly, including the energetic differences arising from temperature and pressure conditions.

When this option is used, no GCMC simulation is required. GALP operates directly on the machine learned probability grids produced by the model.


## Requirements and Assumptions

### Crystallographic Symmetry

As with all GALP workflows, input structures **must be provided in P1 symmetry**. This requirement applies to both the structural CIF file and the force field definitions used in the energy evaluation step.

### Additional Dependencies

In addition to the core GALP dependencies, the machine learned probability distribution workflow requires the following Python packages:

- **DeepAPD**  
  Provides the trained neural network model used to generate adsorbate probability distributions. Installation instructions are available at:  
  https://github.com/uowoolab/DeepAPD

- **ASE** (Atomic Simulation Environment), version **3.26.0**

ASE can be installed using either pip or conda:

```bash
pip install ase==3.26.0
```
or
```bash
conda install -c conda-forge ase=3.26.0
```

## Input Files

To run GALP using machine learned probability distributions, the following files are required.

### Mandatory

- `FIELD`  
  Contains force field definitions and the reference guest molecular geometry. The format is identical to that used for FastMC workflows. The guest must be specified explicitly under the `&guest` section. Currently supported guests are:
  - Xe  
    `&guest Xe`
  - CH<sub>4</sub>  
    `&guest O`

- `<Structure>.cif`  
  CIF file containing the host structure for which adsorbate probability distributions and binding sites will be generated. The structure must be in P1 symmetry.

No grid folding or fold factor needs to be specified for this workflow. The machine learning model automatically generates probability distributions at the unit cell defined by the input structure.

## GALP Input Configuration

In the `GALP.inp` file, the GCMC engine must be set to **FASTMC** in the general input section, as the current implementation relies on the FastMC style `FIELD` file format.

In addition, the user must explicitly request the machine learned probability distribution option and specify the temperature and pressure conditions. At present, the following guest and condition combinations are supported:

- Xe at 1 bar
- CH<sub>4</sub> at 1 bar
- CH<sub>4</sub> at 65 bar

Users must specify the desired settings in the `ML APD Generation` section of the `GALP.inp` input file. Requests outside of these supported conditions will result in an error.

## Example Workflows

Example workflows demonstrating the use of machine learned probability distributions with GALP are provided in the repository:

- Example: [sbmof-1](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/ml_apd/sbmof-1)

The example includes all required input files, the corresponding GALP.inp, and the resulting GALP outputs. [Runnable version](https://github.com/uowoolab/GALP/tree/main/examples/runnable_examples/ml_apd/sbmof-1).

## Notes and Limitations

The machine learned probability distribution workflow is intended as an alternative to GCMC for supported guest and condition combinations. The quality of the predicted probability grids reflects the training domain of the underlying model. Unsupported guests or thermodynamic conditions will raise an error at runtime.
