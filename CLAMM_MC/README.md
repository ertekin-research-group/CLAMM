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
- A Linux system (tested on Ubuntu)

### Build Steps

```bash
# Clone this repository or unpack the source code
cd CLAMM_MC

# Create a build directory
mkdir build && cd build

# Configure the project using CMake
cmake ..

# Compile the code
make
```

The executable `CLAMM_MC` will be generated in the `build` directory.

---

## 🚀 Usage

### Required Files
To run a simulation with CLAMM_MC, you need:

1. `POSCAR` — Defines the unit or supercell geometry and sublattice constraints
2. `INPUT` — Parameter file specifying the simulation type and conditions
3. `CLUSTERS` — Output from CLAMM_Fit containing model parameters
4. (Optional) `SPIN_STATES`, `SRO_DEFINITION`, etc., depending on your simulation

### Example Run

```bash
./CLAMM_MC
```

Make sure the required input files are in the same directory or modify the code to provide full paths.

---

## 🧪 Tutorial Example

### Step 1: Prepare Your POSCAR
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

### Step 2: Write an INPUT File
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

### Step 3: Provide SPIN_STATES
```text
-1 0 1
-1 0 1
0
```

### Step 4: Run the Simulation
```bash
./CLAMM_MC
```

The code will output:
- `OUTPUT` — thermodynamic data (energy, magnetization, heat capacity, etc.)
- (optional) `CONTCAR_xxx` files if enabled

---

## 🧠 Algorithms (ALGO Flag)

| ALGO | Description |
|------|-------------|
| -2   | SRO generation with 3-body Ising support |
| -1   | Targeted short-range order (SRO) generation |
| 0    | Evaluate Hamiltonian only |
| 1    | Magnetic degrees of freedom only |
| 2    | Atomic degrees of freedom only |
| 3    | Magnetic + atomic degrees of freedom |
| 4    | Magnetic + atomic with vacancies/interstitials |

---

## 📁 File Summary

| File         | Purpose                                      |
|--------------|----------------------------------------------|
| `POSCAR`     | Lattice structure and sublattice constraints |
| `INPUT`      | Main simulation parameters                   |
| `CLUSTERS`   | Fitted cluster model parameters              |
| `SPIN_STATES`| Allowed spin states per element              |
| `OUTPUT`     | Simulation thermodynamic results             |
| `OUTPUT_SRO` | Short-range order parameter results (optional)|
| `CONTCAR_*`  | Snapshot of simulation configurations         |

---

## 📚 Citation

If you use CLAMM_MC or the CLAMM toolkit in your research, please cite:

> Brian J. Blankenau, Tianyu Su, Namhoon Kim, Elif Ertekin, "CLAMM: a generalized spin cluster expansion–Monte Carlo toolkit for Alloys and Magnetic Materials", *Computer Physics Communications*, 2025.

---

## 📜 License

CLAMM_MC is licensed under the MIT License. See `LICENSE` for details.

---

## 🤝 Contributions

We welcome contributions! If you would like to add new algorithms, improve performance, or extend CLAMM_MC's functionality, feel free to fork this repository and submit a pull request.

---

## 🔗 Project Homepage

The full CLAMM toolkit (Prep, Fit, MC) and latest tutorials are hosted on:

**[https://github.com/ertekin-research-group/clamm](https://github.com/ertekin-research-group/clamm)**
