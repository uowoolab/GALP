# Usage

This document describes how to run GALP for both small scale analyses and high throughput workflows. GALP can be executed directly within a working directory or launched from an external location, making it suitable for local testing as well as large scale batch processing on HPC systems.

The same execution logic applies to both FastMC and RASPA based workflows. Engine specific requirements are documented separately.

## Small Scale Execution (Manual)

This mode is intended for testing, debugging, or analyzing a small number of systems.

### Preparation

Ensure that the following files are available:

- `GALP.py`
- `GALP.inp`

Both files may reside in the same directory, or `GALP.py` may be located elsewhere and called via a relative or absolute path.

<p align="center">
  <img src="https://raw.githubusercontent.com/uowoolab/GALP/main/images/Small_Scalle.png" alt="Small scale workflow" width="100%">
</p>

#### Execution Steps
1. Open the `GALP.inp` file.
2. Find the comment line `>> Absolute path to directory with input files (Inputs + Probability plots) <<`.
3. Below this comment, insert the absolute path of the directory where you plan to run `GALP.py`.
   - Example: To execute GALP in the `MOF1` directory, update the path to the absolute path where `MOF1` is located.
4. Save and close the `GALP.inp` file.
5. In the terminal, navigate to the directory containing both `GALP.inp` and `GALP.py`.
6. Execute the command:
   ```bash
   python GALP.py
   ```
7. GALP will process using the specified directory settings.

---

### 2. High Throughput (Using sbatch)

#### Preparation
For extensive analysis, utilize sbatch to enable parallel execution and optimize computation time. In this method, leave the absolute path in the `GALP.inp` file empty.

<p align="center">
  <img src="https://github.com/uowoolab/GALP/blob/main/images/Large_Scale.png" alt="Large scale workflow" width="100%">
</p>

#### Execution Steps
1. Ensure the `GALP.py` script is accessible.
2. In the terminal, change to the directory intended for processing.
   - Example: To process the `MOF1` directory, point your terminal to the `MOF1` directory.
3. Execute the command, substituting `path/to/code/` with the actual path leading to `GALP.py`:
   ```bash
   python ../path/to/code/GALP.py
   ```

4. Since the absolute path was left blank in `GALP.inp`, the script defaults to the current working directory for execution.
5. For multiple directories like `MOF1, MOF2, MOF3`, etc., you may set up a `.sh` script to initiate separate jobs for each, utilizing the power of parallel processing.

### Important Tips
1. Always confirm your paths and directories are orderly to avoid errors and promote an efficient execution.
2. Monitor the output and logs to track the processing status and identify any potential issues requiring troubleshooting.
