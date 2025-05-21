# GenMC_Fit ReadMe

## Overview
GenMC_Fit is the second component of the GenMC-MA toolkit, designed to parameterize lattice models from density functional theory (DFT) datasets. It operates as part of the broader GenMC-MA workflow, which is used to model the thermodynamic, magnetic, and structural properties of magnetic alloys and compounds. This tool specifically takes DFT-derived datasets and calculates effective cluster interaction (ECI) terms or other model parameters through a regression approach. These parameters are then used to predict material properties and support Monte Carlo simulations.

---

## Key Features
- Parameterizes lattice models, including cluster expansions, Ising models, and Potts models (More to come).
- Supports multiple linear regression techniques: Ridge, Lasso, Least Squares, and Elastic Net (More to come).
- Handles symmetry-based clustering and decorations for complex lattice systems.
- Outputs parameterized models for use with GenMC_Run.

---

## Prerequisites
- **Python 3.7+**
- Required Python libraries: `numpy`, `pymatgen`, `json`, `yaml`
- DFT datasets with corresponding POSCAR, OUTCAR, and CONTCAR files
- An input parameter file (`param_in`)

---

## Installation
1. Clone the repository containing GenMC_Fit:
   ```bash
   git clone <repository_url>
   ```
2. Install required Python dependencies:
   ```bash
   pip install numpy pymatgen yaml
   ```

---

## Input Files
### 1. **`param_in` File**
This file defines configuration options for GenMC_Fit. Below is an example structure:

```yaml
lat_in: 'POSCAR'  # Path to the lattice structure file
data_file: 'data_file.json'  # Path to compacted DFT dataset
clust_in: 'cluster_in.json'  # Path to cluster definition file
species: ['Ni', 'Mn', 'In']  # Atomic species in the system
fit_ridge: true  # Use Ridge regression for parameter fitting
fit_lasso: false  # Use LASSO regression
fit_eln: false  # Use Elastic Net regression
rescale_enrg: false  # Energy rescaling option
do_fit: true  # Perform regression fitting
do_count: false  # Count clusters (set to true if counts are not precomputed)
```

### 2. **Cluster Definition File (`cluster_in.json`)**
This file defines the motifs and types of clusters to be analyzed. An a simple cluster_in file could look like this:

```json
{
  "List": [
    [ [0, 0, 0], [1], [0] ],
    [ [ [0, 0, 0], [2.6, 0.0, 0.0] ], [1], [1] ]
  ]
}
```

### 3. **Defining Clusters**
Clusters in GenMC_Fit are defined by their geometric motifs and decorations. A motif represents the spatial arrangement of atoms, while a decoration specifies how atomic species are assigned to each position in the motif. For example:

- **Motifs:** These are the unique geometric arrangements of atoms (e.g., a single site, pairs of sites, triangles).
- **Decorations:** These define the specific atomic species (or spins) assigned to each site in the motif.

#### Example:
A motif with three sites located at [0,0,0], [2.6,0,0], and [0,2.6,0] can have different decorations based on species, such as:
- [Ni, Mn, In]
- [Ni, Ni, Mn]

Clusters are specified in the `cluster_in.json` file using the following format:

```json
{
  "List": [
    [ [0, 0, 0], [1], [0] ],  // Single site motif
    [ [ [0, 0, 0], [2.6, 0.0, 0.0] ], [1], [1] ],  // Pair motif
    [ [ [0, 0, 0], [2.6, 0.0, 0.0], [0.0, 2.6, 0.0] ], [1], [2] ]  // Triangle motif
  ]
}
```
- The first list defines the geometric positions.
- The second list defines lattice scaling factors or symmetry flags.
- The third list identifies the cluster type (e.g., Ising (1), Potts (2), or cluster expansion (0) ).

---

## Running GenMC_Fit
1. Prepare the required input files (`param_in`, `cluster_in`, etc.).
2. Run the script using the following command:

   ```bash
   python main.py
   ```

3. The script performs two main tasks depending on the `param_in` settings:
   - **Cluster Counting:** When `do_count: true`, cluster occurrences are computed for the input dataset.
   - **Model Fitting:** When `do_fit: true`, the regression models are applied to fit the lattice model parameters.

---

## Outputs
1. **`CLUSTERS` File**: Contains the parameterized lattice model.
2. **Cluster Counts**: Written to `count_out` if `do_count` is enabled. Counting clusters for each DFT simulation is the most time consuming step in GenMC_Fit. The `count_out` file saves a version of counted clusters for each simulation which can be loaded when trying diffrent regression methods on the same set of clusters. It is highly recomended to 

---

## Example Workflow
### Step 1: Cluster Counting
Update `param_in`:
```yaml
# main
do_count: True                          # define whether to count or not
do_fit: True                           # define whether to fit or not
lat_in: 'POSCAR'                        # name of input lattice file (can be any pymatgen readable format)
data_file: 'NiMnIn_data.txt'                # name of input structure metadata
clust_in: 'cluster_in'                  # name of input cluster file
species: ['Ni', 'Mn', 'In']             # list of input species
fit_lasso: True                         # define fitting method
fit_ridge: False
fit_eln: False
rescale_enrg: False                      # define fitting with rescaled energy or raw energy

# fit
sample_times: 1000                      # bootstrap sampling times
sample_ratio: 0.8                       # bootstrap sampling size over input data size
kfold: 10                               # number of data splitting folds
alpha_range: [-6, 2]                    # range of alpha (need to test carefully)
l1_ratio: [.4, .5, .6, .7, .9]          # range of l1_ratio in ElasticNet (need to test)
convergence: 1e-5                       # tolerance for the coefficient optimization in Lasso/ElasticNet
```
Run the script:
```bash
python main.py
```
Output: `count_out` (contains computed cluster counts).

### Step 2: Try a New Regression
Update `param_in`:
Change `do_count` to False.
```yaml
# main
do_count: False                          # define whether to count or not
do_fit: True                           # define whether to fit or not
lat_in: 'POSCAR'                        # name of input lattice file (can be any pymatgen readable format)
data_file: 'NiMnIn_data.txt'                # name of input structure metadata
clust_in: 'cluster_in'                  # name of input cluster file
species: ['Ni', 'Mn', 'In']             # list of input species
fit_lasso: False                         # define fitting method
fit_ridge: True
fit_eln: False
rescale_enrg: False                      # define fitting with rescaled energy or raw energy

# fit
sample_times: 1000                      # bootstrap sampling times
sample_ratio: 0.8                       # bootstrap sampling size over input data size
kfold: 10                               # number of data splitting folds
alpha_range: [-6, 2]                    # range of alpha (need to test carefully)
l1_ratio: [.4, .5, .6, .7, .9]          # range of l1_ratio in ElasticNet (need to test)
convergence: 1e-5                       # tolerance for the coefficient optimization in Lasso/ElasticNet
```
Run the script:
```bash
python main.py
```
Output: `eci_out` (contains fitted parameters).
---

