# GALP Input File Reference (`GALP.inp`)

This document provides a detailed reference for the `GALP.inp` input file. Each keyword is described in terms of its purpose, effect on the workflow, and relationship to the generated output files. This document is intended as a technical reference rather than a step-by-step guide.

## Contents

1. [General Input](#general-input)
2. [Binding Site Algorithm](#binding-site-algorithm)
3. [MD Parameters](#md-parameters)
4. [GCMC Analysis Parameters](#gcmc-analysis-parameters)
5. [Machine Learned APD Generation](#machine-learned-apd-generation)
6. [Post Process Cleaning](#post-process-cleaning)

## General Input

### Software Choice

**Input line:**  
`Software Choice: [FASTMC, RASPA]`

**Type:**  
String (required)

**Description:**  
Selects the GCMC engine used to generate the probability distributions. This option determines which parsing routines are invoked, as the input formats and metadata differ substantially between FastMC and RASPA.

---

### Input Files Path

**Input line:**  
`Input Files Path:`

**Type:**  
String (path, optional)

**Description:**  
Specifies the directory containing the input files to be processed. This may be an absolute or relative path. If left blank, GALP defaults to using the current working directory.

### Maximum Number of Binding Sites

**Input line:**  
`Maximum Number of Binding Sites (integer, optional):`

**Type:**  
Integer (optional, default: 100)

**Description:**  
Sets the maximum number of binding sites retained in the output files. After candidate sites are identified and ranked by relative occupancy, only the top *n* sites are kept.

This option is useful when working with large supercells or unfolded probability grids, where many low-occupancy sites may be identified.

### Maximum Binding Sites Energy

**Input line:**  
`Maximum Binding Sites Energy (integer, optional, default: -100 kJ/mol):`

**Type:**  
Integer (optional, default: −100 kJ/mol)

**Description:**  
Prunes binding sites whose binding energy is greater than −*n* kJ/mol. By default, sites weaker than −100 kJ/mol are discarded. Positive binding energies are always removed.

Energies above this threshold are generally unrealistic in practice, particularly for chemisorbed systems.

### Generate CIF for High Occupancy Sites

**Input line:**  
`Generate CIF for high occupancy sites (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, writes a CIF file containing the extracted local maxima prior to molecular fitting or pruning. This file is useful for diagnosing the quality of the probability distributions.

When enabled, the file `<Guest>_local_maxima.cif` is generated.

---

## Binding Site Algorithm

### Algorithm Selection

**Input line:**  
`Algorithm Selection: [0: Legacy GALA (vectorial), 1: RMSD] (default: 1):`

**Type:**  
Integer (default: `1`)

**Description:**  
Selects the binding site fitting algorithm. The legacy vectorial method corresponds to the original implementation. The RMSD based algorithm is more robust and has been validated across a wider range of guest molecules, making it the recommended option.

### Overlap Tolerance

**Input line:**  
`Overlap Tolerance (float, default: 0.35):`

**Type:**  
Float (default: `0.35`)

**Description:**  
Controls how candidate maxima are grouped prior to fitting.

For the legacy algorithm, this value acts as an absolute cutoff when comparing bond vectors. For the RMSD based algorithm, it is used as a preprocessing filter to reduce the number of candidate combinations passed to the RMSD fitting step.

### RMSD Cutoff

**Input line:**  
`RMSD Cutoff (float, default: 0.1):`

**Type:**  
Float (default: `0.1`)

**Description:**  
Sets the rotational tolerance for molecular fitting using the RMSD based algorithm. Smaller values enforce stricter matching between the reference geometry and candidate configurations.

This option has no effect when the legacy algorithm is selected.

### Exclude Hydrogen in Fitting Algorithm

**Input line:**  
`Exclude hydrogen in fitting algorithm (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, hydrogen atoms are excluded during molecular fitting. This option is useful for organic guests where hydrogen Lennard Jones parameters are often set to zero or where coarse-grained models are employed.

---

## MD Parameters

### MD Program Executable

**Input line:**  
`MD Program Executable:`

**Type:**  
String (path)

**Description:**  
Path to the DL_POLY Classic executable used for binding energy evaluation and optional binding site optimization.

### Optimize Binding Sites (OBS)

**Input line:**  
`Optimize Binding Sites (OBS) [T/F] (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, performs a molecular dynamics optimization of each binding site using DL_POLY. This refines guest positions prior to energy evaluation.

When enabled, the file `<Guest>_binding_sites_optimized.cif` is generated.

### Optimization Step

**Input line:**  
`Optimization Step (default: 1000):`

**Type:**  
Integer (default: `1000`, required if OBS is `T`)

**Description:**  
Number of MD steps used during optimization. Increasing this value beyond approximately 10 000 generally yields diminishing returns.

### Timestep

**Input line:**  
`Timestep (ps):`

**Type:**  
Float  
Default: `0.0` if OBS is `F`, `0.001` if OBS is `T`

**Description:**  
MD timestep used during optimization. This parameter is sensitive. Excessively large timesteps may cause DL_POLY failures, while very small values may result in negligible relaxation.

### Cutoff

**Input line:**  
`Cutoff (Angstrom) (default: 12.5):`

**Type:**  
Float (Å, default: `12.5`)

**Description:**  
Nonbonded interaction cutoff used for Lennard Jones and Coulombic interactions during MD evaluation.

### Dedicated CPUs for MD Simulation

**Input line:**  
`Dedicated CPUs for MD simulation (default: 1):`

**Type:**  
Integer (default: `1`)

**Description:**  
Number of CPU cores allocated for multiprocessing of DL_POLY calculations. Increasing this value can significantly reduce runtime for systems with many binding sites.

---

## GCMC Analysis Parameters

### Sigma

**Input line:**  
`Sigma (default: 0.4):`

**Type:**  
Float (default: `0.4`)

**Description:**  
Controls the width of Gaussian smoothing applied to the probability density grids.

### Radius

**Input line:**  
`Radius (default: 0.45):`

**Type:**  
Float (default: `0.45`)

**Description:**  
Defines the exclusion radius around each identified maximum. Maxima within this radius of a higher-occupancy maximum are removed to reduce redundancy.

### Cutoff

**Input line:**  
`Cutoff (default: 0.1):`

**Type:**  
Float (default: `0.1`)

**Description:**  
Relative occupancy threshold used to prune low-occupancy maxima prior to molecular fitting.

### Write Folded

**Input line:**  
`Write Folded [T/F] (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, folded probability grids are written to cube files ending with `_folded.cube`.

### Grid Factors

**Input line:**  
`Grid Factors: <nx ny nz>`

**Type:**  
Integer array (default: `1 1 1`)

**Description:**  
Specifies folding factors applied to the probability grids. Required if folding is enabled or unfolded grids are provided. These values may be inferred automatically from engine-specific input files when available.

### Write Smoothed

**Input line:**  
`Write Smoothed [T/F] (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, smoothed unit cell probability grids are written to files ending with `_smooth.cube`.

---

## Machine Learned APD Generation

### Generate Adsorbate Probability Distribution

**Input line:**  
`Generate adsorbate probability distribution [T/F] (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
Enables the machine learned probability distribution workflow using DeepAPD. Only a limited set of guests and thermodynamic conditions are currently supported.

### Temperature and Pressure

**Input line:**  
`Temperature and pressure (T, p) corresponding to ML APD:`

**Type:**  
Tuple `(T, p)`

**Description:**  
Thermodynamic conditions used for machine learned probability generation.

**Supported combinations:**
- Xe at 1 bar  
- CH₄ at 1 bar  
- CH₄ at 65 bar  

---

## Post Process Cleaning

### Remove DL_POLY_BS Directory

**Input line:**  
`Remove DL_POLY_BS directory [T/F] (default: F):`

**Type:**  
Boolean (`T` / `F`, default: `F`)

**Description:**  
If enabled, removes the `DL_POLY_BS` directory after successful execution to reduce disk usage.
