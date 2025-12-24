# FastMC Workflow

This document describes how to run GALP using adsorbate probability distributions generated with FastMC.

FastMC v1.4.0 is a DL_POLY Classic based Grand Canonical Monte Carlo code written in Fortran and optimized for adsorption simulations in porous crystalline materials. The input format closely follows DL_POLY conventions, with `FIELD` and `CONTROL` files analogous to those used for molecular dynamics simulations. Detailed documentation and example cases for FastMC are available in the FastMC repository:

https://github.com/uowoolab/FastMC-1.4.0

## Requirements and Assumptions

### Crystallographic Symmetry

All FastMC based workflows used with GALP **must employ structures in P1 symmetry**.

GALP assumes a one to one mapping between probability grids and Cartesian coordinates. If the input structure is not in P1, GALP will successfully generate the `<Guest>_binding_sites.cif` file but will fail during the DL_POLY based energy evaluation step. In this case, the file `<Guest>_galp_binding_sites.xyz`, which contains binding energies and relative occupancies, will not be produced.

To convert structures to P1 symmetry, you may use:
- [CSD-cleaner](https://github.com/uowoolab/CSD-cleaner)
- [pymatgen](https://github.com/materialsproject/pymatgen)
- [ASE](https://github.com/DeepChoudhuri/Atomic-Simulation-Environment) 

## Input Files

To run GALP with FastMC generated data, the following files are required:

### Mandatory
- `FIELD`  
  Contains force field definitions and guest molecular geometry.
- `Prob_Guest_X_Site_Y.cube` or `Prob_Guest_X_Site_Y_folded.cube`  
  Volumetric adsorbate probability distributions files. A separate file must be provided for each guest and each site present in the system, following the naming convention shown above. If the expected filename is unclear, running GALP once will report the required guest and site specific filenames in `galp.log`, which can be used as a reference to ensure correct naming.

### Optional but Recommended
- `CONTROL`  
  If provided, GALP will automatically parse:
  - the MD cutoff value
  - the GCMC grid factor used for probability grid folding

If the CONTROL file is not present, GALP handles the MD cutoff and grid spacing independently. If the MD cutoff is not specified in GALP.inp, it defaults to 12.5 Å. If the grid spacing is not specified, grid folding defaults to treating the supercell as a single unit cell (1 x 1 x 1). When grid folding is requested under these conditions, the unfolded supercell is therefore used directly.


## GALP Input Configuration

In the `GALP.inp` file, the GCMC engine must be set to FastMC in the general input section.

Only the `FIELD` and the adsorbate probability distributions files are strictly required. However, for high throughput studies, providing the `CONTROL` file is strongly recommended to ensure consistent cutoff handling and grid folding behavior across systems.

## Folded vs Unfolded Probability Grids

When unfolded probability grids are used, the number of identified binding sites scales with the folding factor. While this may be useful for diagnostic purposes, it significantly increases the runtime of GALP and the number of configurations passed to the DL_POLY energy evaluation step.

For large scale studies, folded probability grids are recommended unless explicitly required.

## Example Workflows

Example FastMC based workflows are provided in the repository:

- Example 1: [cald-15](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/fastmc/calf-15)
- Example 2: [calf-16](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/fastmc/calf-16)
- Example 3: [calf-20](https://github.com/uowoolab/GALP/tree/main/examples/example_outputs/fastmc/calf-20)

Each example includes the required FastMC output files, a corresponding GALP.inp, and the resulting GALP outputs. [Runnable versions](https://github.com/uowoolab/GALP/tree/main/examples/runnable_examples/fastmc) of these examples are provided for users who wish to reproduce the calculations.
