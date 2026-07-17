# conn_def_polyNet
This project investigates how local connectivity defects affect building block mobility in metallo-supramolecular polymer networks. By introducing specific flaws into model star polymer networks, we aim to uncover the mechanisms behind defect-accelerated relaxation to design tunable transient networks

## Workflow & Directory Structure

The project code is organized sequentially to match the simulation pipeline:

- `src/01_generation/`: Contains two files `generate_heteroleptic.py` and `generate_homoleptic.py` to build the initial network configurations for the two systems (`init_file.gsd`).
- `src/02_simulation/`: Contains `run_sim.py` to run MD simulations in HOOMD-blue and output trajectories (`raw_data.gsd`).
- `src/03_analysis/`: Contains `analyze.py` to calculate network properties (MSD, $R_{\text{g}}$, diffusion $D$, bond dissociation $c(t)$, Rouse modes, etc.) and generate plots.

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
