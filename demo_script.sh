# Demo file to demonstrate the input and outputs of this code.
# The outputs of the commands listed below are in the Examples folder

# Section ONE
# Equilibrium Intrinsic Defect in NiAl (mp-1487)

# Generate Intrinsic Defect Structure
gen_def_structure --mpid mp-1487 --mapi_key <put_your_mapi_key_here>
# Notes: To obtain your mapi key, refer to "API keys" section at 
# https://materialsproject.org/open
# The command creates a folder mp-1487 with 5 subfolders containing
# input files for first principles DFT calculations with VASP. If you
# have VASP installed and set the POTCAR directory in pymatgen, 
# POTCAR files are also generated.

# After successfully running the DFT calculations in the five 
# subfolders of mp-1487, execute 
gen_def_energy --mpid mp-1487 --mapi_key <put_your_mapi_key_here>
# Notes: This should create mp-1487_raw_defect_energy.json file

# After mp-1487_raw_defect_energy.json is ready
gen_def_profile --mpid mp-1487 --temp 1000
# Output Files: 1) mp-1487_chem_pot.dat 2) mp-1487_def_concentration.dat 3) mp-1487_def_energy.dat 


# Section TWO
# Solute Mo in NiAl (mp-1487)

# Generate Solute Defect Structure
gen_sol_pref_structure --mpid mp-1487 --solute Mo --mapi_key <put_your_mapi_key_here>

# After DFT calculations in two folders are done
gen_sol_def_energy --mpid mp-1487 --solute Mo --mapi_key <put_your_mapi_key_here>

# After mp-1487_solute-Mo_raw_defect_energy.json is ready
gen_sol_site_pref --mpid mp-1487 --solute Mo --temp 1000
# Output File: 1) mp-1487_solute-Mo_def_concentration.dat
