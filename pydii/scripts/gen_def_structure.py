#!/usr/bin/env python
"""
This file reads arguments and generate vacancy and antisite structures 
in intermetallics.
"""

__author__ = "Bharat Medasani, Enze Chen"
__data__  = "Apr 28, 2020"

import os
import sys
from argparse import ArgumentParser 
import json

from pymatgen.io.vasp.sets import MPMetalRelaxSet
from pymatgen.ext.matproj import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site
from pymatgen.analysis.defects.core import Vacancy
from pymatgen.io.vasp.inputs import Kpoints


def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1./3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(mpid, mapi_key, cellmax):
    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return

    # Get conventional and primitive structures from the Materials Project DB
    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_ratio = int(conv_struct_sites / prim_struct_sites)

    # Default VASP settings
    blk_vasp_incar_param = {'IBRION':-1, 'EDIFF':1e-5, 'EDIFFG':0.001, 'NSW':0,}
    def_vasp_incar_param = {'ISIF':2, 'EDIFF':1e-6, 'EDIFFG':0.001,}
    kpoint_den = 6000
    
    # Create bulk structure and associated VASP files
    sc_scale = get_sc_scale(inp_struct=struct, final_site_no=cellmax)
    blk_sc = struct.copy()
    blk_sc.make_supercell(scaling_matrix=sc_scale)
    site_no = blk_sc.num_sites

    # Rescale if needed
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        blk_sc = struct.copy()
        blk_sc.make_supercell(scaling_matrix=sc_scale)
    blk_str_sites = set(blk_sc.sites)
    custom_kpoints = Kpoints.automatic_density(blk_sc, kppa=kpoint_den)
    mpvis = MPMetalRelaxSet(blk_sc, user_incar_settings=blk_vasp_incar_param,
                        user_kpoints_settings=custom_kpoints)
    ptcr_flag = True
    try:
        potcar = mpvis.potcar
    except:
        print ("VASP POTCAR folder not detected.\n" \
               "Only INCAR, POSCAR, KPOINTS are generated.\n" \
               "If you have VASP installed on this system, \n" \
               "refer to pymatgen documentation for configuring the settings.")
        ptcr_flag = False
    fin_dir = os.path.join(mpid, 'bulk')
    mpvis.write_input(fin_dir)

    # Create each defect structure and associated VASP files
    for i, site in enumerate(struct.sites):
        vac = Vacancy(structure=struct, defect_site=site)
        vac_sc = vac.generate_defect_structure(supercell=sc_scale)
        site_no = vac_sc.num_sites
        
        # Get vacancy site information
        vac_str_sites = set(vac_sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = next(iter(vac_sites))
        site_mult = int(vac.get_multiplicity() / conv_prim_ratio)
        vac_site_specie = vac_site.specie
        vac_symbol = vac_site_specie.symbol

        custom_kpoints = Kpoints.automatic_density(vac_sc, kppa=kpoint_den)
        mpvis = MPMetalRelaxSet(vac_sc, user_incar_settings=def_vasp_incar_param,
                        user_kpoints_settings=custom_kpoints)
        vac_dir = 'vacancy_{}_mult-{}_sitespecie-{}'.format(
                    str(i+1), site_mult, vac_symbol)
        fin_dir = os.path.join(mpid, vac_dir)
        mpvis.write_input(fin_dir)

        # Antisites generation at the vacancy site
        struct_species = blk_sc.species
        for specie in set(struct_species) - set([vac_site_specie]):
            specie_symbol = specie.symbol
            anti_sc = vac_sc.copy()
            anti_sc.append(specie, vac_site.frac_coords)
            mpvis = MPMetalRelaxSet(anti_sc, user_incar_settings=def_vasp_incar_param,
                        user_kpoints_settings=custom_kpoints)
            anti_dir = 'antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                    str(i+1), site_mult, vac_symbol, specie_symbol)
            fin_dir = os.path.join(mpid, anti_dir)
            mpvis.write_input(fin_dir)


def substitute_def_struct_gen(mpid, solute, mapi_key, cellmax):
    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return
    if not solute:
        print ("============\nERROR: Provide solute atom\n============")
        return
    
    # Get conventional and primitive structures from the Materials Project DB
    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_ratio = int(conv_struct_sites / prim_struct_sites)

    # Default VASP settings
    def_vasp_incar_param = {'ISIF':2, 'EDIFF':1e-6, 'EDIFFG':0.001,}
    kpoint_den = 6000
    
    # Create each substitutional defect structure and associated VASP files
    sc_scale = get_sc_scale(inp_struct=struct, final_site_no=cellmax)
    blk_sc = struct.copy()
    blk_sc.make_supercell(scaling_matrix=sc_scale)
    site_no = blk_sc.num_sites

    # Rescale if needed
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        blk_sc = struct.copy()
        blk_sc.make_supercell(scaling_matrix=sc_scale)

    blk_str_sites = set(blk_sc.sites)
    for i, site in enumerate(struct.sites):
        vac = Vacancy(structure=struct, defect_site=site)
        vac_sc = vac.generate_defect_structure(supercell=sc_scale)
        site_no = vac_sc.num_sites
        
        # Get vacancy site information
        vac_str_sites = set(vac_sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = next(iter(vac_sites))
        vac_specie = vac_site.specie.symbol
        site_mult = int(vac.get_multiplicity() / conv_prim_ratio)

        # Solute substitution defect generation at the vacancy site
        solute_struct = vac_sc.copy()
        solute_struct.append(solute, vac_site.frac_coords)
        custom_kpoints = Kpoints.automatic_density(solute_struct, kppa=kpoint_den)
        mpvis = MPMetalRelaxSet(solute_struct, user_incar_settings=def_vasp_incar_param,
                              user_kpoints_settings=custom_kpoints)
        
        # Check if POTCAR file can be generated
        ptcr_flag = True
        try:
            potcar = mpvis.potcar
        except:
            print ("VASP POTCAR folder not detected.\n" \
                   "Only INCAR, POSCAR, KPOINTS are generated.\n" \
                   "If you have VASP installed on this system, \n" \
                   "refer to pymatgen documentation for configuring the settings.")
            ptcr_flag = False

        # Generate VASP directory
        sub_def_dir ='solute_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                str(i+1), site_mult, vac_specie, solute)
        fin_dir = os.path.join(mpid, sub_def_dir)
        mpvis.write_input(fin_dir)


def im_vac_antisite_def_struct_gen():
    m_description = 'Command to generate vacancy and antisite defect ' \
                    'structures for intermetallics.' 

    parser = ArgumentParser(description=m_description)

    parser.add_argument("--mpid",
            type=str.lower,
            help="Materials Project id of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument("--mapi_key",
            default = None,
            help="Your Materials Project REST API key.\n" \
                 "For more info, please refer to " \
                 "www.materialsproject.org/opne")

    parser.add_argument("--cellmax",
            type=int,
            default=128,
            help="Maximum number of atoms in supercell.\n" \
                 "The default is 128\n" \
                 "Keep in mind the number of atoms in the supercell" \
                 "may vary from the provided number including the default.")
    
    args = parser.parse_args()
    vac_antisite_def_struct_gen(args.mpid, args.mapi_key, args.cellmax)


def im_sol_sub_def_struct_gen():
    m_description = 'Command to generate solute substitution defect ' \
                    'structures for intermetallics.' 

    parser = ArgumentParser(description=m_description)

    parser.add_argument("--mpid",
            type=str.lower,
            help="Materials Project id of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument("--solute", 
            type=str,
            help="Solute Element")

    parser.add_argument("--mapi_key",
            default = None,
            help="Your Materials Project REST API key.\n" \
                 "For more info, please refer to " \
                 "www.materialsproject.org/opne")

    parser.add_argument("--cellmax",
            type=int,
            default=128,
            help="Maximum number of atoms in supercell.\n" \
                 "The default is 128\n" \
                 "Keep in mind the number of atoms in the supercell" \
                 "may vary from the provided number including the default.")
    
    args = parser.parse_args()
    substitute_def_struct_gen(args.mpid, args.solute, args.mapi_key, args.cellmax)


if __name__ == '__main__':
    #im_vac_antisite_def_struct_gen()
    im_sol_sub_def_struct_gen()
