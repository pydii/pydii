=====
PyDII
=====

Equilibrium defect concentrations and solute site preference in intermetallic compounds.

Requirements
------------
PyDII requires Python 3.4+, `Pymatgen <https://pymatgen.org/>`_, `SymPy <https://www.sympy.org/en/index.html>`_, and `Monty <https://guide.materialsvirtuallab.org/monty/>`_ packages and optionally `VASP <https://www.vasp.at>`_ for full functionality. 

Source Code
------------
Obtaining
~~~~~~~~~
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_ if it's not already packaged with your system.

#. Download the PyDII source code using the command::

    git clone https://github.com/pydii/pydii.git

Code Description
~~~~~~~~~~~
The PyDII source code is in the folder ``ROOT_FOLDR/pydii``. The ``ROOT_FOLDR``
is also named ``pydii``. The root folder contains additional files for setup and one
``examples`` folder that can be used to test the code.

Within the ``pydii`` folder, the essential functionality is implemented in
``dilute_solution_model.py``. Scripts to generate input files for VASP calculations,
parse VASP runs, and generate defect concentrations are in the ``scripts`` subfolder.
The scripts

#.  Generate VASP input files for defect supercell calculations.
#.  Parse the output files from VASP calculations.
#.  Generate defect concentration, defect formation energy, and chemical potential data.

To test the programs fully, VASP is required. If a VASP installation is not found,
a POTCAR file is not generated in step 1 and step 2 cannot be performed. In case
the user uses alternative density functional theory packages for first-principles
calculations, the input to step 3 should be prepared in the JSON format given
in ``examples/NiAl_mp-1487/mp-1487_raw_defect_energy.json``. Contact the authors
for any additional details.

Installation
------------
#. Navigate to the PyDII root directory ::

    cd pydii

#. Install the code, using the command ::

    python setup.py install

The command tries to obtain the required packages and their dependencies
and install them automatically. Access to root may be needed if
``virtualenv`` is not used. Alternatively you can try ::

   pip install . --user

#. The package can be installed at non-standard locations using the command ::

    python setup.py install --prefix PYDII_ROOTDIR

where ``PYDII_ROOTDIR`` is your choice of directory. In UNIX/Linux environments,
add ``PYDII_ROOTDIR`` to ``PATH`` and ``PYTHONPATH`` variables by the following
commands, which you can add to your ``.bashrc`` file ::

    export PATH=$PATH:PYDII_ROOTDIR
    export PYTHONPATH=$PYTHONPATH:PYDII_ROOTDIR

Materials Project API key
-------------------------
PyDII makes use of the Materials Project (MP) database to get the crystal structure of
the intermetallic. To access the MP database, a MP API key is required. Please refer
to https://materialsproject.org/open to obtain your MP API key.


Examples
--------

From the PyDII root directory, go to the ``examples`` folder by typing::

    cd examples

The ``examples`` folder contains two subfolders for Al3V and NiAl intermetallic systems.
For a description of the compounds and how to generate the defect concentration profiles,
refer to the manuscript: 
`Ding, H. et al, Computer Physics Communications, 193, (2015) <https://www.sciencedirect.com/science/article/pii/S0010465515001149>`_.

Version 2.0
----------

Starting with version 1.9 (pre-release version of 2.0), the signature of the parser data is changing. 
So data files produced with PyDII codes with different versions may not work. 
Please report any bugs on the Issues tracker.

Optionally, instead of specifying ``mpid``, users can use their own structures to run the PyDII calculations, 
as long as the structure is `compatible with Pymatgen <https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure.from_file>`_. 
Use the ``-h`` flag for more information.
