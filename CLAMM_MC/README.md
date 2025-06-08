# CLAMM_MC: Monte Carlo Engine for Predicting Alloy Properties

**CLAMM_MC** is a C++ Monte Carlo simulation engine designed to predict thermodynamic, magnetic, and structural properties of complex alloy systems and magnetic materials. It is one of three core components of the CLAMM toolkit, which integrates Density Functional Theory (DFT) data, cluster expansion (CE), and spin-lattice models into a streamlined computational pipeline.

This package is part of the larger CLAMM framework:
- `CLAMM_Prep` (Python): Prepares DFT simulation data
- `CLAMM_Fit` (Python): Fits lattice and spin-lattice models
- `CLAMM_MC` (C++): Runs Monte Carlo simulations on fitted models

---

## 🔧 Build Instructions

### Prerequisites
To compile and run `CLAMM_MC`, you need:

- C++ compiler with C++11 or later support
- CMake ≥ 3.10
 
CLAMM_MC has been tested on Mac Linux and Windows.

### Build Steps
Assuming you have already cloned the main CLAMM repository, navigate the CLAMM_MC folder and create a build directory:
```
cd ./CLAMM_MC
# Create a build directory
mkdir build
cd build
```
Inside the build directory  configure the project using cmake:
```
cmake ..
```
And build it with:
```
cmake --build .
```
The executable `CLAMM_MC.exe` will be generated in the `CLAMM_MC/bin/Debug` directory.
For convenience we will copy the .exe to the Sample_Input folder in preperation for running the tutorial:
```
cd ../
cp ./bin/Debug/CLAMM_MC.exe ./Sample_Input
```

---

## 🚀 Usage

### Required Files
To run a simulation with CLAMM_MC, you need the folloing input files:

1. `POSCAR` — Defines the unit or supercell geometry and sublattice constraints
2. `INPUT` — Parameter file specifying the simulation type and conditions
3. `CLUSTERS` — Output from CLAMM_Fit containing model parameters
4. (Optional) `SPIN_STATES`, `SRO_DEFINITION`, etc., depending on your simulation

### Example Run
Once the nessisary input files are supplied simply run:
```bash
./CLAMM_MC
```
Make sure the required input files are in the same directory or modify the `INPUT` file to provide the full paths.

---

## 🧪 Tutorial Example
### Step 1: Input File Prep
Prepaire all input files. Example input files are included in the Sample_Input folder.
### Prepare Your POSCAR
CLAMM_MC uses a POSCAR-like file to define the simulation cell. The only differnece between CLAMM and VASP POSCAR files lies in the implementation of sub-lattice information. In order to define any sub-lattices CLAMM_MC requires that each set of atomic coordinates is followed by the atomic spieces that are allowed to occupy that, or an equivilant site.
```text
Ni2MnIn
1.0
6.01 0.00 0.00
0.00 6.01 0.00
0.00 0.00 6.01
Ni Mn In
8 4 4
Direct
0.25 0.25 0.75 Ni
...
0.50 0.50 0.50 Mn In
```

### Write an INPUT File
The INPUT file is where the user defines all simulation specific paramiters. An example INPUT file is shown below:
```text
ALGO = 3
STRUCTURE = POSCAR
USE_POSCAR = TRUE
SHAPE = 2 2 2
ATOM_NUMBS = 1728 864 864
SPECIES = Ni Mn In
SIM_TYPE = DEFAULT
SPIN_INIT = FM
TA_PASSES = 100
EQ_PASSES = 300
START_TEMP = 4000
END_TEMP = 0
TEMP_STEP = 40
USE_STATES = TRUE
WRITE_CONTCARS = FALSE
```
The effect of each of these lines is shown in the following table:

| Line | Flag | Purpose |
|------|------|---------|
| 1 | `ALGO = 3` | Selects the simulation algorithm:<br>- `3` enables a full Monte Carlo simulation allowing **both atomic and magnetic configurations** to change. |
| 2 | `STRUCTURE = POSCAR` | Specifies the filename for the unit cell or initial configuration. Must be a modified VASP-style `POSCAR`. |
| 3 | `USE_POSCAR = TRUE` | Indicates that the program should use the `POSCAR` file as input (as opposed to generating one). |
| 4 | `SHAPE = 2 2 2` | Sets the demensions in x, y, z to build a **supercell** from the unit cell (e.g., 2×2×2). |
| 5 | `ATOM_NUMBS = 1728 864 864` | Target **number of atoms of each species** in the supercell, matching order in `SPECIES`. |
| 6 | `SPECIES = Ni Mn In` | Lists the **chemical elements** (species) used in the simulation. Order matters and must match POSCAR. |
| 7 | `SIM_TYPE = DEFAULT` | Placeholder for selecting simulation type (currently always set to `DEFAULT`; may allow future extensions). |
| 8 | `SPIN_INIT = FM` | Sets the **initial magnetic configuration**:<br>- `FM`: Ferromagnetic (all spins aligned)<br>- `AFM`: Antiferromagnetic<br>- `RAND`: Random spins |
| 9 | `TA_PASSES = 100` | Number of **thermalization (equilibration)** passes performed before measuring observables. |
| 10 | `EQ_PASSES = 300` | Number of **passes used to average physical quantities** (after thermalization). |
| 11 | `START_TEMP = 4000` | Starting **temperature** for simulated annealing (in K or arbitrary units). |
| 12 | `END_TEMP = 0` | Ending **temperature** for simulated annealing. |
| 13 | `TEMP_STEP = 40` | **Temperature step** size for simulated annealing. |
| 14 | `SRO_TARGET = 1.0` | (Only used when `ALGO = -1`) Target value for **short-range order** parameter (Warren-Cowley). |
| 15 | `USE_STATES = TRUE` | If `TRUE`, spins are **restricted to user-defined allowed values** (e.g., from `SPIN_STATES` file). |
| 16 | `WRITE_CONTCARS = FALSE` | If `TRUE`, will output **CONTCAR-style snapshots** after each temperature step. If `FALSE`, only the final state is saved. |

Options for the ALGO Flag:
| ALGO | Description |
|------|-------------|
| -2   | SRO generation with 3-body Ising support |
| -1   | Targeted short-range order (SRO) generation |
| 0    | Evaluate Hamiltonian only |
| 1    | Magnetic degrees of freedom only |
| 2    | Atomic degrees of freedom only |
| 3    | Magnetic + atomic degrees of freedom |
| 4    | Magnetic + atomic with vacancies/interstitials |


### Provide SPIN_STATES
```text
-1 0 1
-1 0 1
0
```

### Step 2: Run the Simulation
Yep. Just two steps. Its just that easy. Run the simulation with:
```bash
./CLAMM_MC
```

The as the simulation runs, it will create the following output files:
- `OUTPUT` — thermodynamic data (energy, magnetization, heat capacity, etc.)
- (optional) `CONTCAR_xxx` files if enabled

---

## 📁 File Summary

| File         | Purpose                                                           |
|--------------|-------------------------------------------------------------------|
| `POSCAR`     | Lattice structure and sublattice constraints                      |
| `INPUT`      | Main simulation parameters                                        |
| `CLUSTERS`   | Fitted cluster model parameters                                   |
| `SPIN_STATES`| Allowed spin states per element                                   |
| `OUTPUT`     | Simulation order paramiters and thermodynamic results             |
| `OUTPUT_SRO` | Short-range order parameter results (optional)                    |
| `CONTCAR_*`  | Snapshot of simulation configurations                             |

---

