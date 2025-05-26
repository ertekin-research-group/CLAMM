![title](media/Clamm1.png)
# Wellcome to CLAMM
A Cluster Expansion and Monte Carlo toolkit for Alloys and Magnetic Materials 

CLAMM consists of three main tools:
- `CLAMM_Prep` (Python): Prepares DFT simulation data
- `CLAMM_Fit` (Python): Fits lattice and spin-lattice models
- `CLAMM_MC` (C++): Runs Monte Carlo simulations on fitted models



**CLAMM** is an open-source software suite for simulating thermodynamic, magnetic, and structural properties of complex alloys and magnetic materials using lattice models informed by DFT data. It provides tools for data preparation, model fitting, and Monte Carlo simulation through its three main components:

- **CLAMM Prep** – Python tool for converting VASP-based DFT data into a compact model-ready format.
- **CLAMM Fit** – Python tool for parameterizing cluster expansion and spin-lattice models.
- **CLAMM MC** – C++ Monte Carlo engine for performing simulations with the parametrized models.

---

## 🌐 Overview

CLAMM allows researchers to:
- Generate special quasi-random structures (SQS)
- Explore short-range and long-range ordering
- Predict magnetic and phase transitions
- Simulate alloy stability across temperature and composition

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

```bash
cd prep
pip install -r requirements.txt

cd ../fit
pip install -r requirements.txt
```

Requirements include:
- `numpy`, `pymatgen`, `scikit-learn`, `spglib`

### C++ Tool (CLAMM_MC)

```bash
cd mc
mkdir build && cd build
cmake ..
make
```

Generates `CLAMM_MC` executable in the `build` directory.

---

## 🚀 Example Workflow

### 1. Prepare DFT Data
Use `CLAMM Prep` to collate VASP simulation data:
```bash
python clamm_prep.py param_in.yaml
```

### 2. Fit a Model
Use `CLAMM Fit` to parameterize cluster models:
```bash
python clamm_fit.py param_in.yaml
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
- `prep/README.md` – DFT data parsing and formatting
- `fit/README.md` – Model definition and fitting
- `mc/README.md` – Monte Carlo flags and usage

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
