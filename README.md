# conn_def_polyNet
This project investigates how local connectivity defects affect building block mobility in metallo-supramolecular polymer networks. By introducing specific flaws into model star polymer networks, we aim to uncover the mechanisms behind defect-accelerated relaxation to design tunable transient networks

## Workflow & Directory Structure

The project code is organized sequentially to match the simulation pipeline:

- `src/01_generation/`: Contains `generate.py` to build the initial network configurations (`init_file.gsd`).
- `src/02_simulation/`: Contains `run_sim.py` to run MD simulations in HOOMD-blue and output trajectories (`raw_data.gsd`).
- `src/03_analysis/`: Contains `analyze.py` to calculate network properties (MSD, $R_g$, diffusion $D$, bond dissociation $c(t)$, Rouse modes, etc.) and generate plots.

## Setup & Installation

To install the required dependencies (HOOMD-blue, GSD, NumPy, SciPy, Matplotlib), use the provided Conda environment file:

```bash
conda env create -f environment.yml
conda activate polymer-networks
