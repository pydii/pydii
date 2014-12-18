=====
pydii
=====

Equilibrium defect concentrations and solute site preference in intermetallic compounds

Installation
------------
pydii requires pymatgen, sympy and monty packages. 

Source Code
------------
Obtaining
~~~~~~~~~
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pydii source code using the command::

    git clone https://github.com/pydii/pydii.git
Description
~~~~~~~~~~~
PyDII source code is in the folder ROOT_FOLDR/pydii. The ROOT_FOLDR 
is also named pydii. The root folder additional files for setup and one Examples 
folder that can be used to test run the code. 

Within the pydii folder, the essential functionality is implemented in 
dilute_solution_model.py. Scripts to generate input files for VASP calculations,
parsing vasp runs, and generate defect concentrations are in the scripts subfolder.

Installation
------------
#. Navigate to pydii root directory::

    cd pydii

#. Install the code, using the command::

    python setup.py install

The command tries to obtain the required packages and their dependencies and install them automatically.
Access to root may be needed if ``virtualenv`` is not used.

Examples
--------

From the pydii, go to examples folder by typing::

    cd examples

The examples folder contains two subfolders for Al3V and NiAl intermetallic systems. For a description of
the compounds and how to generate the defect concentration profiles, refer to the manuscript. 



