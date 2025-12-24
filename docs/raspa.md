# RASPA3 Workflow

This document describes how to run GALP using adsorbate probability distributions generated with RASPA3.

GALP supports RASPA3 version **3.0.18 or newer**, which introduced equitable binning and site resolved adsorption probability distributions. These features are required for compatibility with GALP and are not available in earlier RASPA3 releases.

RASPA3 is an open source Monte Carlo simulation package for adsorption and diffusion in porous materials. The source code, documentation, and examples are available at:

https://github.com/iRASPA/RASPA3

## Required RASPA3 Options

GALP relies on two RASPA3 input options introduced in version 3.0.18:

- `"DensityGridBinning"`
- `"DensityGridPseudoAtomsList"`

These options enable equitable binning and site specific probability density output, respectively. Detailed descriptions of these options can be found in the RASPA3 command manual:

[commands.md](https://github.com/iRASPA/RASPA3/blob/main/docs/manual/commands.md#density-grids)

Example input files demonstrating the required density grid setup are provided in the RASPA3 documentation:

[examples_density_grids.md](https://github.com/iRASPA/RASPA3/blob/main/docs/manual/examples_density_grids.md)  
[density_grids_examples](https://github.com/iRASPA/RASPA3/tree/main/examples/density_grids)

## Requirements and Assumptions

### Crystallographic Symmetry

All RASPA3 based workflow used with GALP **must employ structures in P1 symmetry**.

In addition, when preparing RASPA3 input files, the `force_field.json` file must also reflect a P1 representation, where each atom in the structure is explicitly defined. This is required to ensure a one to one mapping between probability grid sites and Cartesian coordinates.

If the structure is not in P1, GALP will generate the `<Guest>_binding_sites.cif` file but will fail during the DL_POLY based energy evaluation step, and the file `<Guest>_galp_binding_sites.xyz`, containing binding energies and relative occupancies, will not be produced.

To convert structures to P1 symmetry, you may use:
- [CSD-cleaner](https://github.com/uowoolab/CSD-cleaner)
- [pymatgen](https://github.com/materialsproject/pymatgen)
- [ASE](https://github.com/DeepChoudhuri/Atomic-Simulation-Environment) 

### Atomic Charges

For the MD based energy evaluation step in GALP, partial charges must be specified for each atom in the host structure within `force_field.json`. These charges are used to compute Coulombic interactions during the DL_POLY calculations.

If charges are not provided, all framework charges default to **0.0**, and binding energies will be computed using Lennard Jones interactions only. In this case, the `<Guest>_binding_sites.cif` file is unaffected, but the reported binding energies will not include electrostatic contributions.

## Input Files

To run GALP using RASPA3 generated data, the following files are required.

### Mandatory

- `simulation.json`  
  Contains system level information such as cutoff distances, grid factors, and the definition of host and guest components.

- `force_field.json`  
  Defines force field parameters and atom labels used for mapping probability grids to atomic sites. Must correspond to a P1 representation of the structure.

- `<Guest>.json`  
  Contains force field parameters for the guest molecule. For binary or multicomponent simulations, multiple guest JSON files must be provided.

- `Prob_Guest_X_Site_Y.cube` or `Prob_Guest_X_Site_Y_folded.cube`  
  Adsorbate probability distributions grids for each guest and site. RASPA3 generates these in a proprietary cube format, which must be converted to a pymatgen compatible format prior to running GALP.   

  A simple command line conversion tool is available here:  
  https://github.com/uowoolab/raspa3-cube-to-pymatgen  
  This conversion step can be easily integrated into high throughput workflows.

## GALP Input Configuration

In the `GALP.inp` file, the GCMC engine must be set to **RASPA** in the general input section. No additional engine specific configuration is required beyond the presence of the converted probability grid files.


## Folded vs Unfolded Probability Grids

When unfolded probability grids are used, the number of identified binding sites scales with the folding factor. While this may be useful for diagnostic purposes, it significantly increases the runtime of GALP and the number of configurations passed to the DL_POLY energy evaluation step.

For large scale studies, folded probability grids are recommended unless explicitly required.

## Example Workflows

Example RASPA3 based workflows are provided in the repository:

- Example 1: [calf-20 binary](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/raspa/calf-20_binary)
- Example 2: [calf-20](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/raspa/calf-20_unary)
- Example 3: [cpl-1](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/raspa/cpl-1)

Each example includes the required RASPA3 output files, a corresponding `GALP.inp`, and the resulting GALP outputs. [Runnable versions](https://github.com/uowoolab/GALP/tree/main/examples/runnable_examples/raspa) of these examples will be provided for users who wish to reproduce the calculations.

