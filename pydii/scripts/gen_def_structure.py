#!/usr/bin/env python
"""
This file reads arguments and generate vacancy and antisite structures 
in intermetallics.
"""

from __future__ import division

__author__ = "Bharat Medasani"
__data__  = "Sept 14, 2014"

import os
import sys
from argparse import ArgumentParser 
import json

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.ext.matproj import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site
#from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import Vacancy, Substitution, Interstitial
from pymatgen.io.vasp.inputs import Kpoints


def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(mpid, mapi_key, cellmax):
    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return

    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)

    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites/prim_struct_sites)
    sc_scale = get_sc_scale(inp_struct=struct, final_site_no=cellmax)

    mpvis = MVLRelax52Set(struct)

    # Begin defaults: All default settings.
    blk_vasp_incar_param = {'IBRION':-1,'EDIFF':1e-4,'EDIFFG':0.001,'NSW':0,}
    def_vasp_incar_param = {'ISIF':2,'NELM':99,'IBRION':2,'EDIFF':1e-6, 
                            'EDIFFG':0.001,'NSW':40,}
    kpoint_den = 6000
    # End defaults
    
    ptcr_flag = True
    try:
        potcar = mpvis.get_potcar(struct)
    except:
        print ("VASP POTCAR folder not detected.\n" \
              "Only INCAR, POSCAR, KPOINTS are generated.\n" \
              "If you have VASP installed on this system, \n" \
              "refer to pymatgen documentation for configuring the settings.")
        ptcr_flag = False


    vac = Vacancy(struct, {})
    scs = vac.make_supercells_with_defects(sc_scale)
    site_no = scs[0].num_sites
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        scs = vac.make_supercells_with_defects(sc_scale)

    for i in range(len(scs)):
        sc = scs[i]
        poscar = mpvis.get_poscar(sc)
        kpoints = Kpoints.automatic_density(sc,kpoint_den)
        incar = mpvis.get_incar(sc)
        if ptcr_flag:
            potcar = mpvis.get_potcar(sc)

        interdir = mpid
        if not i:
            fin_dir = os.path.join(interdir,'bulk')
            try:
                os.makedirs(fin_dir)
            except:
                pass
            incar.update(blk_vasp_incar_param)
            incar.write_file(os.path.join(fin_dir,'INCAR'))
            poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            if ptcr_flag:
                potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        else:
            blk_str_sites = set(scs[0].sites)
            vac_str_sites = set(sc.sites)
            vac_sites = blk_str_sites - vac_str_sites
            vac_site = list(vac_sites)[0]
            site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
            vac_site_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol

            vac_dir ='vacancy_{}_mult-{}_sitespecie-{}'.format(str(i),
                    site_mult, vac_symbol)
            fin_dir = os.path.join(interdir,vac_dir)
            try:
                os.makedirs(fin_dir)
            except:
                pass
            incar.update(def_vasp_incar_param)
            poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            incar.write_file(os.path.join(fin_dir,'INCAR'))
            if ptcr_flag:
                potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))

            # Antisite generation at all vacancy sites
            struct_species = scs[0].types_of_specie
            for specie in set(struct_species)-set([vac_site_specie]):
                subspecie_symbol = specie.symbol
                anti_struct = sc.copy()
                anti_struct.append(specie, vac_site.frac_coords)
                poscar = mpvis.get_poscar(anti_struct)
                incar = mpvis.get_incar(anti_struct)
                incar.update(def_vasp_incar_param)
                as_dir ='antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                        str(i), site_mult, vac_symbol, subspecie_symbol)
                fin_dir = os.path.join(interdir,as_dir)
                try:
                    os.makedirs(fin_dir)
                except:
                    pass
                poscar.write_file(os.path.join(fin_dir,'POSCAR'))
                incar.write_file(os.path.join(fin_dir,'INCAR'))
                if ptcr_flag:
                        potcar.write_file(os.path.join(fin_dir,'POTCAR'))
                kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))


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
    def_vasp_incar_param = {'ISIF':2, 'ISMEAR':1, 'SIGMA':0.1, 'NSW':60, 
                            'EDIFF':1e-5, 'EDIFFG':0.001,}
    kpoint_den = 6000
    
    # Create each substitutional defect structure and VASP files
    sc_scale = get_sc_scale(inp_struct=struct, final_site_no=cellmax)
    blk_sc = struct.copy()
    blk_sc.make_supercell(scaling_matrix=sc_scale)
    blk_str_sites = set(blk_sc.sites)
    for i, site in enumerate(struct.sites):
        vac = Vacancy(structure=struct, defect_site=site)
        vac_sc = vac.generate_defect_structure(supercell=sc_scale)
        site_no = vac_sc.num_sites

        # Rescale if needed
        if site_no > cellmax:
            max_sc_dim = max(sc_scale)
            i = sc_scale.index(max_sc_dim)
            sc_scale[i] -= 1
            vac_sc = vac.generate_defect_structure(supercell=sc_scale)
        
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
        mpvis = MPRelaxSet(solute_struct, user_incar_settings=def_vasp_incar_param,
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
    im_vac_antisite_def_struct_gen()
    #im_sol_sub_def_struct_gen()


    '''
    vac = Vacancy(structure=struct, defect_site=struct.sites[0])
    sc_scale = get_sc_scale(inp_struct=struct, final_site_no=cellmax)
    scs = vac.generate_defect_structure(supercell=sc_scale)
    site_no = scs.num_sites
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        scs = vac.generate_defect_structure(supercell=sc_scale)

    interdir = mpid
    blk_str_sites = set(scs.sites)
    for i in range(1, len(scs)):
        sc = scs[i]
        vac_str_sites = set(sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = list(vac_sites)[0]
        site_mult = int(vac.get_defectsite_multiplicity(i-1) / conv_prim_rat)
        vac_site_specie = vac_site.specie
        vac_specie = vac_site.specie.symbol

        # Solute substitution defect generation at all vacancy sites
        struct_species = scs[0].types_of_specie
        solute_struct = sc.copy()
        solute_struct.append(solute, vac_site.frac_coords)
        '''
    '''
        incar = mpvis.incar
        incar.update(def_vasp_incar_param)
        poscar = mpvis.poscar
        kpoints = Kpoints.automatic_density(solute_struct, kpoint_den)
        if ptcr_flag:
            potcar = mpvis.potcar

        try:
            os.makedirs(fin_dir)
        except:
            pass
        poscar.write_file(os.path.join(fin_dir,'POSCAR'))
        incar.write_file(os.path.join(fin_dir,'INCAR'))
        kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        if ptcr_flag:
            potcar.write_file(os.path.join(fin_dir,'POTCAR'))
    '''
