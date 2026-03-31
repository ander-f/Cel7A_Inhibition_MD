# Molecular Basis of Cel7A Inhibition by Glucose, Mannose, and Xylose

This repository contains the simulation parameters, analysis scripts, and structural data for the study of **TrCel7A** inhibition (both Carbohydrate-Binding Module (CBM) and the Catalytic Domain (CD)) by monosaccharides (glucose, mannose, and xylose).

## 📂 Repository Structure

* `ff/`: Force field parameters for the protein and sugars (Glucose, Mannose, and Xylose).
* `mdp_files/`: GROMACS input parameters for energy minimization, equilibration (NVT/NPT), and production runs.
* `pdbs/`: Starting structures, selected solvated systems, and representative conformational ensembles for visualization.
* `scripts/`: 
    * `cm_*.jl`: Julia scripts for solvation analysis using [**ComplexMixtures.jl**](https://github.com/m3g/ComplexMixtures.jl).
    * `Project.toml` & `Manifest.toml`: Julia environment configuration for full reproducibility.
    * `packmolinputcreator_*.jl`: Julia scripts for packmol input generation.
    * `plot_fig*.jl`: Julia scripts to generate the figures presented in the manuscript.
    * `run_folder.sh`: Shell script for directory generation.
* `analyses/`: Processed data files.
* `figures/`: Final plots of the results (PNG/SVG).

## 🛠 Reproducibility

### 1. Molecular Dynamics
Simulations were performed using GROMACS 2023. Production parameters for both CD and CBM systems are provided in the `mdp_files/` directory. Starting configurations for each component of the studied systems are located in `pdbs/`.

### 2. Analysis (Julia)
To reproduce the solvation and structural analyses, you will need [Julia](https://julialang.org/) installed.

1.  **Navigate to the scripts directory:**
    ```bash
    cd scripts
    ```
2.  **Activate and instantiate the environment:**
    ```julia
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    ```
    *This command automatically installs all required packages (ComplexMixtures, Plots, etc.) at the specific versions used in this study.*

3.  **Run a specific script (e.g., Plot Figure 1):**
    ```bash
    julia --project=. plot_fig1.jl
    ```
