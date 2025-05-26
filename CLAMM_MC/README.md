# CLAMM_MC

**CLAMM_MC** is the C++ Monte Carlo engine of the CLAMM toolkit — a Generalized Spin Cluster Expansion–Monte Carlo toolkit for Alloys and Magnetic Materials. This component performs Monte Carlo simulations using models built from DFT data, enabling predictions of thermodynamic, magnetic, and structural properties of complex alloys.

## Features

- Supports Ising, Potts, and decorated cluster expansion models
- Enables prediction of phase diagrams, order parameters, and SRO
- Efficient and modular C++ implementation
- Designed for finite-temperature simulations of magnetic and atomic disorder

## Directory Structure

```
CLAMM_MC/
├── CMakeLists.txt       # Build configuration
├── src/                 # Source files for CLAMM_MC
│   ├── algo*.cpp/h      # Algorithm-specific MC engines
│   ├── main.cpp         # Main program entry point
│   ├── session.cpp/h    # Input/session management
│   ├── sim_cell.cpp/h   # Simulation cell and lattice structure
│   └── ...              # Supporting utilities
```

## Building the Code

### Prerequisites

- CMake (>=3.10)
- C++17 compatible compiler (e.g., `g++`, `clang++`)

### Build Instructions

```bash
# Navigate to the source directory
cd CLAMM_MC

# Create a build directory and configure
mkdir build && cd build
cmake ..

# Compile
make
```

The executable will be generated in the `build/` directory (e.g., `CLAMM_MC`).

## How to Use

CLAMM_MC requires the following inputs:

1. A parameterized model file (`CLUSTERS`)
2. A VASP-style `POSCAR` file with added sublattice labels
3. An `INPUT` file specifying simulation parameters

### Example Workflow

```bash
# Step 1: Place input files in the working directory
cp ../examples/POSCAR .
cp ../examples/CLUSTERS .
cp ../examples/INPUT .

# Step 2: Run the Monte Carlo simulation
./CLAMM_MC
```

### Sample `INPUT` file

```ini
ALGO = 3
STRUCTURE = POSCAR
USE_POSCAR = TRUE
SHAPE = 1 1 1
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

Refer to the [CLAMM Documentation](https://github.com/ertekin-research-group/xxx) for more examples and detailed parameter descriptions.

## License

This project is licensed under the MIT License.

## Citation

If you use CLAMM in your work, please cite:

> Blankenau, B.J., Su, T., Kim, N., Ertekin, E. (2025). *CLAMM: A generalized spin cluster expansion–Monte Carlo toolkit for Alloys and Magnetic Materials.* Computer Physics Communications (in review).
