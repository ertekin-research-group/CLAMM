![title](media/Clamm1.png)
# Wellcome to CLAMM: A Cluster Expansion and Monte Carlo toolkit for Alloys and Magnetic Materials 

**CLAMM** is an open-source software suite for simulating thermodynamic, magnetic, and structural properties of complex alloys and magnetic materials using lattice models informed by DFT data. It provides tools for data preparation, model fitting, and Monte Carlo simulation through its three main components:

- **CLAMM Prep** – Python tool for converting VASP-based DFT data into a compact model-ready format.
- **CLAMM Fit** – Python tool for parameterizing cluster expansion and spin-lattice models.
- **CLAMM MC** – C++ Monte Carlo engine for performing simulations with the parametrized models.
---
  
## Key Features
- **Cluster expansion** and **spin-lattice models** for configurational and magnetic modeling.
- Support for **multi-sublattice systems** and **user-defined motifs**.
- Generation of **special quasi-random structures (SQS)** and **configurations with target short-range order (SRO)**.
- Monte Carlo algorithms with flexible temperature schedules and spin/atomic configuration support.
- Interpretable models using the **decorated cluster expansion (DCE)** formalism.
---

## 📦 Project Structure
The CLAMM repository is devided into seperate folders for each tool
```
CLAMM/
├── CLAMM_Prep/         # CLAMM Prep (Python)
├── CLAMM_Fit/          # CLAMM Fit (Python)
├── CLAMM_MC/           # CLAMM MC (C++)
```
In addition to the code, each folder contains all sample files needed for running a short tutorial. Detailed can be found in the files README.

---

## 🔧 Installation Instructions

### Python Tools (Prep and Fit)
To install all python dependencies, navigate to the main CLAMM directory and run:
```bash
pip install -r requirements.txt
```

### C++ Tool (CLAMM_MC)
To build the c++ code nessisary to perform Monte Carlo simulations, navigate to the main CLAMM directory and run: 
```bash
cd CLAMM_MC
mkdir build && cd build
cmake ..
make
```
This will generates a `CLAMM_MC` executable in the `build` directory.
---

## 🚀 Example Workflow

### 1. Prepare DFT Data
Use `CLAMM Prep` to collate VASP simulation data:
```bash
python clamm_prep.py
```
### 2. Fit a Model
Use `CLAMM Fit` to parameterize cluster models:
```bash
python clamm_fit.py
```
### 3. Run Monte Carlo
Use `CLAMM_MC` to simulate behavior at different temperatures:
```bash
./CLAMM_MC
```
---

## 🤝 Contributions

We welcome contributions! If you would like to add new algorithms, improve performance, or extend CLAMM's functionality, feel free to fork this repository and submit a pull request.

---

## 📘 Documentation

Each module has detailed documentation and examples:
- `CLAMM_Prep/README.md` – DFT data parsing and formatting
- `CLAMM_Fit/README.md` – Model definition and fitting
- `CLAMM_MC/README.md` – Lattice Model initalization and Monte Carlo simulation
A detailed explenation of all CLAMM tools can be found at ... (the link to the paper when published)
---

## 🔗 Repository

The latest version is hosted at:

**[https://github.com/ertekin-research-group/clamm](https://github.com/ertekin-research-group/clamm)**

---
## 📜 License

CLAMM is released under the MIT License.

---

## 📚 Citation

If you use CLAMM in your work, please cite:

> Blankenau, B.J., Su, T., Kim, N., Ertekin, E. (2025). *CLAMM: A generalized spin cluster expansion–Monte Carlo toolkit for Alloys and Magnetic Materials.*(in review).
