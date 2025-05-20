# GenMC_Preprocess

GenMC_Preprocess is a component of the Generalized Monte Carlo Toolkit for Magnetic Alloys (GenMC-MA). It is a Python-based preprocessing tool that facilitates the construction of lattice models by collating, verifying, and formatting data generated from DFT (Density Functional Theory) simulations. The processed data can then be used to parameterize models for simulating magnetic and alloy configurations. Currently, GenMC_Preprocess only supports VASP input files for the DFT dataset.

---

## Features

- **Collates DFT simulation results**: Automatically gathers results of DFT simulations into a single data file.
- **Formats data for lattice models**: Transforms DFT outputs into a compact format suitable for lattice model parameterization
- **Validates simulation results**: Identifies and flags incomplete or corrupt DFT results.
- **Analyze magnetic moment distribution**: Finds the distribution of magnetic moments for each constituent atomic species.

---

## Workflow

The general workflow for using GenMC_Preprocess involves the following steps:
1. Create a dataset from DFT simulations. 
    - Any model created with GenMC-MA is only as good as its dataset. In order to ensure the best possible mode, the dataset should include a diverse set of atomic and magnetic configurations. It is recommended that, for each unique POSCAR file in a dataset should be initialized with as many symmetrically unique configurations of magnetic moments as is practical.

    - The directory structure of the dataset should be organized as follows:
```plaintext
    |
    |_ Root Directory (Name)
	|
	|_ Any sub-directory organization (Optional)
		|
		|_ VASP Sumulation-1 Results (Name)
		|	|
		|	|_ POSCAR
		|	|_ CONTCAR
		|	|_ OUTCAR
		|
		|_ VASP Sumulation-2 Results (Name)
		|	|
		|	|_ POSCAR
		|	|_ CONTCAR
		|	|_ OUTCAR
		|...
   ```
   Each simulation results folder MUST contain at least a POSCAR, CONTCAR, and OUTCAR file. 
   Any additional VASP input or output files are welcome to remain in the folder but will be ignored. 
	
2. Format the "param_in" file.
   - An example "param_in" file is shown below:

	```plaintext
	root_dir: '/Users/Desktop/AlloyData'              # the top directory of your VASP data
	do_compile_data: True                             # Choose whether or not to write the dataset to a single file
	write_data: 'data_mag'                            # Path and name of your output dataset file
	species: ['Ni', 'Mn', 'In']                       # Species in sequence as they appear in your POSCAR file
	read_mag: True                                    # Flag to read magnetic moment from OUTCAR
	use_spin_tol: True                                # Choose whether display spin as integers or total magnetic moment
	spin_tol: [[0.2] , [3.0 , 2.0] , [ -1]]           # Spin tolerance here in the same sequence of species
	do_mag_distrib: True                              # Choose whether or not to calculate magnetic moment distribution
	write_mag_distrib: 'mag_distrib'                  # Path and name of your magnetic distribution output file
	```
   - The "spin_tol" option requires additional explanation. When the "use_spin_tol" flag is set to “True”, a spin tolerance method is used. In this case, magnetic moments read from the OUTCAR file are assigned integer 
     values such as −1, 0, 1 or −1, 1. The mapping between the OUTCAR moments and the integer values and is provided by "spin_tol" using a tolarance vector. Each element in the vector represents how the tolarance is 
     applied to the corisponding atomic species. In the example "param_in" file, the tolarance vector indicates that there are two allowable spin magnitudes for Ni atoms:
     - If ∥OUTCAR moment∥ > 0.2 then ∥spin∥ is set to 1 and its sign is set to sgn(OUTCAR moment)
     - Otherwise ∥spin∥ is set to 0
      
     The second element in the vector indicates that there are three allowable magnitudes for spins assigned to Mn atoms:
     - If ∥OUTCAR moment∥ > 3.0 then ∥spin∥ is set to 2 and its sign is set to sgn(OUTCAR moment)
     - If ∥OUTCAR moment∥ < 3.0 and > 2.0 then ∥spin∥ is set to 1 and its sign is set to sgn(OUTCAR moment)
     - Otherwise ∥spin∥ is set to 0.
    
     For the third element in the vector a special value is used, -1. This indicates that, for all atoms of this type (in this case In) all spins are set to 0.0.

     __Choosing the Spin Tolerance__
     It is important to carefully choose the "spin_tol" values in order to achieve a suitable model. If the "do_mag_distrib" flag is set to "True", GenMC_Preprocess will return the magnetic moments of each atom organized by atomic species. This data can be plotted as a histogram to view the distribution of magnetic moments. It is recommended that the chosen tolerance vector reflect the moment distribution in order to achieve the best model.
     

3. Run the GenMC_Preprocess Python script. This will creat a dataset file at the user defined location in a format that is readable by GenMC_Fit, the next tool in the GenMC-MA workflow. The dataset file should look something like this:

		# Ni Mn In
		8 4 4 \Ni2MnIn\Mart\B0 -88.5581 5.2 5.2 7.8 1.5708 1.5708 1.5708
		5.2 0 0
		0 5.2 0
		0 0 7.8
		0 0 1.0 0.25 0.25 0.25
		1 0 1.0 0.75 0.25 0.25
		2 0 1.0 0.25 0.75 0.25
		3 0 1.0 0.75 0.75 0.25
		4 0 1.0 0.25 0.25 0.75
		5 0 1.0 0.75 0.25 0.75
		6 0 1.0 0.25 0.75 0.75
		7 0 1.0 0.75 0.75 0.75
		8 1 1.0 0.5 0.0 0.0
		9 1 1.0 0.0 0.5 0.0
		10 1 1.0 0.5 0.0 0.5
		11 1 1.0 0.0 0.5 0.5
		12 2 0 0.0 0.0 0.0
		13 2 0 0.5 0.5 0.0
		14 2 0 0.0 0.0 0.5
		15 2 0 0.5 0.5 0.5

The first line lists the atomic species present in the dataset in the same order as in the POSCAR files. The second line lists the number of each atom of each species, the file path the data was found at, the DFT total energy in eV, and the lattice constants. The next three lines are the lattice vectors. The remaining lines show the atom index, spin/moment, and fractional coordinates of each atom. 

Feel Free to try and run GenMC_Preprocess on the included sample data for NiMnIn. 


