import os
import glob

from setuptools import setup, find_packages

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
        name="pydii",
        packages=find_packages(),
        version="1.9.0",
        install_requires=["pymatgen>=2020.1", "sympy>=1.0", "numpy>=1.13"],
        author="Bharat Medasani, Enze Chen",
        author_email="mbkumar@gmail.com, chenze@berkeley.edu",
        url="http://github.com/mbkumar/pydii",
        description="PyDII is a python tool to evaluate point defect "
                  "concentrations and micro-alloy solute site preference "
                  "in intermetallics.",
        long_description=readme(),
        classifiers=[
            "Programming Language :: Python :: 3.4",
            "Development Status :: 1 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Software Development :: Libraries :: Python Modules"
            ],
        license="MIT",
        scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*")),
        test_suite='nose.collector',
        tests_require=['nose'],
        entry_points={
            'console_scripts':[
                'gen_def_structure=pydii.scripts.gen_def_structure:'
                               'im_vac_antisite_def_struct_gen',
                'gen_def_energy=pydii.scripts.gen_def_energy:'
                               'im_vac_antisite_def_energy_parse',
                'gen_def_profile=pydii.scripts.gen_def_profile:'
                               'im_vac_antisite_def_profile',
                'gen_sol_pref_structure=pydii.scripts.gen_def_structure:'
                               'im_sol_sub_def_struct_gen',
                'gen_sol_def_energy=pydii.scripts.gen_def_energy:'
                               'im_sol_sub_def_energy_parse',
                'gen_sol_site_pref=pydii.scripts.gen_def_profile:'
                               'im_sol_sub_def_profile']}
)

