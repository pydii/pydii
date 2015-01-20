=====
pydii
=====

Equilibrium defect concentrations and solute site preference in intermetallic compounds

Requirements
------------
pydii requires pymatgen, sympy and monty packages and optionally VASP for full functionality. 

Source Code
------------
Obtaining
~~~~~~~~~
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pydii source code using the command::

    git clone https://github.com/pydii/pydii.git
    
Code Description
~~~~~~~~~~~
PyDII source code is in the folder ROOT_FOLDR/pydii. The ROOT_FOLDR 
is also named pydii. The root folder additional files for setup and one Examples 
folder that can be used to test run the code. 

Within the pydii folder, the essential functionality is implemented in 
dilute_solution_model.py. Scripts to generate input files for VASP calculations,
parsing vasp runs, and generate defect concentrations are in the scripts subfolder.
The scripts 
#.  generate VASP input files for defect supercell calculations,
#.  parse the output files from VASP calculations, and
#.  generate defect concentration, defect formation energy and chemical potential data.
To test the programs fully, VASP is required. If VASP installation is not found, 
POTCAR file is not generated in step 1 and step 2 can not be performed. In case, 
the user uses alternative density functional packages for first principles 
calculations, the input to step 3 should be preparared in the json format given 
in examples/NiAl_mp-1487/mp-1487_raw_defect_energy.json. Contact the authors for 
any addition details.

Installation
------------
#.  Navigate to pydii root directory::

    cd pydii

#.  Install the code, using the command::

    python setup.py install

    The command tries to obtain the required packages and their dependencies 
    and install them automatically. Access to root may be needed if 
    ``virtualenv`` is not used.

#. The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYDII_ROOTDIR

    where PYDII_ROOTDIR is your choice of directory. In UNIX/Linux environments, 
    add PYDII_ROOTDIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYDII_ROOTDIR
    export PYTHONPATH=$PYTHONPATH:PYDII_ROOTDIR    

Materials Project API key
-------------------------



Examples
--------

From the pydii, go to examples folder by typing::

    cd examples

The examples folder contains two subfolders for Al3V and NiAl intermetallic systems. For a description of
the compounds and how to generate the defect concentration profiles, refer to the manuscript. 



