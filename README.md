# conn_def_polyNet
This project investigates how local connectivity defects affect building block mobility in metallo-supramolecular polymer networks. By introducing specific flaws into model star polymer networks, we aim to uncover the mechanisms behind defect-accelerated relaxation to design tunable transient networks

## Project Details

SFB1552

## Detailed Workflow & Code Description

### 1. Network Generation (`src/01_generation/`)
- `generate_heteroleptic.py` & `generate_homoleptic.py`: These scripts define the topology of your star polymer networks. 
  - **Key Variables:** Core functionality ($f$), arm length ($N_{\text{arm}}$), and defect fraction ($p_{\text{defect}}$).
  - **What the code does:** It algorithmically places coordinates for star polymer cores and linear arm monomers in a periodic simulation box, intentionally leaving a percentage of functional groups unconnected to simulate local connectivity defects. It writes the topology out as an initial `.gsd` file.

### 2. Molecular Dynamics Simulation (`src/02_simulation/`)
- `run_sim.py`: Controls the HOOMD-blue simulation pipeline.
  - **Key Variables:** Time-step ($\Delta t$), total integration steps, temperature ($T$), and friction coefficient ($\gamma$).
  - **What the code does:** It initializes the system from the generated `.gsd` file, applies a Weeks-Chandler-Andersen (WCA) potential for excluded volume interactions, sets up a Langevin thermostat for NVT dynamics, runs an energy minimization, and then executes production runs to output trajectory frames.

### 3. Data Analysis (`src/03_analysis/`)
- `analyze.py`: Processes the raw trajectory data.
  - **What the code does:** It loops through the trajectory file (`raw_data.gsd`) to calculate:
    - Mean Squared Displacement (MSD) to determine building block mobility.
    - Diffusion coefficient ($D$) via the Einstein relation.
    - Bond dissociation correlation function $c(t)$ to monitor the lifetime of the supramolecular links.

## Setup & Installation

To install the required dependencies (HOOMD-blue, GSD, NumPy, SciPy, Matplotlib), use the provided Conda environment file:

```bash
conda env create -f environment.yml
conda activate polymer-networks
```
## Publication

If you use this code or find our work helpful, please cite our paper:

* **Title:** Effects of Connectivity Defects on the Structure and Dynamics of Star Polymer Networks
* **Authors:** Sayam Bandyopadhyay, Sebastian Seiffert, and Arash Nikoubashman
* **Journal:** *Macromolecules*, 2026, 59 (3), 1191–1200
* **DOI:** [10.1021/acs.macromol.5c02272](https://doi.org/10.1021/acs.macromol.5c02272)

### BibTeX
```bibtex
@article{doi:10.1021/acs.macromol.5c02272,
  author  = {Bandyopadhyay, Sayam and Seiffert, Sebastian and Nikoubashman, Arash},
  title   = {Effects of Connectivity Defects on the Structure and Dynamics of Star Polymer Networks},
  journal = {Macromolecules},
  volume  = {59},
  number  = {3},
  pages   = {1191-1200},
  year    = {2026},
  doi     = {10.1021/acs.macromol.5c02272}
}
```

## Funding & Acknowledgments

This work is part of the Collaborative Research Center **SFB 1552: Defects and Defect Control in Soft Matter**, funded by the German Research Foundation (DFG).

* **Project ID:** German Research Foundation (DFG) - Project number: 465145163
* **Speaker:** Professor Dr.-Ing. Sebastian Seiffert (Johannes Gutenberg University Mainz)
* **Funding Period:** Since 2023

For more information about the overall initiative and its sub-projects, visit the [SFB 1552 Project Page](https://gepris.dfg.de/gepris/projekt/465145163).
