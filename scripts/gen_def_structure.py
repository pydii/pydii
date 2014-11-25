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

from pymatgen.io.vaspio_set import MPGGAVaspInputSet
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure import Structure
from pymatgen.serializers.json_coders import pmg_load, pmg_dump
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.io.vaspio.vasp_input import Poscar, Kpoints


def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(mpid, mapi_key, cellmax):
    print (mpid, mapi_key, cellmax)
    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return

    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)

    print (struct.formula)
    sc_scale = get_sc_scale(struct,cellmax)

    mpvis = MPGGAVaspInputSet()

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

    print (ptcr_flag)

    vac = Vacancy(struct, {}, {})
    scs = vac.make_supercells_with_defects(sc_scale)
    site_no = scs[0].num_sites
    if site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        scs = vac.make_supercells_with_defects(sc_scale)
    print ("No. of atoms in supercell", scs[0].num_sites)

    for i in range(len(scs)):
        sc = scs[i]
        poscar = mpvis.get_poscar(sc)
        kpoints = Kpoints.automatic_density(sc,kpoint_den)
        incar = mpvis.get_incar(sc)
        if ptcr_flag:
            potcar = mpvis.get_potcar(sc)

        #print sc
        interdir = mpid
        if not i:
            fin_dir = os.path.join(interdir,'bulk')
            try:
                os.makedirs(fin_dir)
            except:
                #print 'Failed creating ',fin_dir
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
            site_mult = vac.get_defectsite_multiplicity(i-1)
            vac_site_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol

            vac_dir ='vacancy_{}_mult-{}_sitespecie-{}'.format(str(i),
                    site_mult, vac_symbol)
            fin_dir = os.path.join(interdir,vac_dir)
            try:
                os.makedirs(fin_dir)
            except:
                #print 'Failed creating ',fin_dir
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
                    #print 'Failed creating ',fin_dir
                    pass
                poscar.write_file(os.path.join(fin_dir,'POSCAR'))
                incar.write_file(os.path.join(fin_dir,'INCAR'))
                if ptcr_flag:
                        potcar.write_file(os.path.join(fin_dir,'POTCAR'))
                kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))

def substitute_def_struct_gen(mpid, solute, mapi_key, cellmax):
    print (mpid, solute, mapi_key, cellmax)
    if not mpid:
        print ("============\nERROR: Provide an mpid\n============")
        return
    if not solute:
        print ("============\nERROR: Provide solute atom\n============")
        return

    if not mapi_key:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            struct = mp.get_structure_by_material_id(mpid)

    print (struct.formula)

    mpvis = MPGGAVaspInputSet()

    # Begin defaults: All default settings.
    blk_vasp_incar_param = {'IBRION':-1,'EDIFF':1e-4,'EDIFFG':0.001,'NSW':0,}
    def_vasp_incar_param = {'ISIF':2,'NELM':99,'IBRION':2,'EDIFF':1e-6, 
                            'EDIFFG':0.001,'NSW':40,}
    kpoint_den = 6000
    # End defaults
    
    # Check if POTCAR file can be geneated
    ptcr_flag = True
    try:
        potcar = mpvis.get_potcar(struct)
    except:
        print ("VASP POTCAR folder not detected.\n" \
              "Only INCAR, POSCAR, KPOINTS are generated.\n" \
              "If you have VASP installed on this system, \n" \
              "refer to pymatgen documentation for configuring the settings.")
        ptcr_flag = False

    print (ptcr_flag)

    vac = Vacancy(struct, {}, {})
    sc_scale = get_sc_scale(struct,cellmax)
    scs = vac.make_supercells_with_defects(sc_scale)
    site_no = scs[0].num_sites
    if site_no > cellmax:
            max_sc_dim = max(sc_scale)
            i = sc_scale.index(max_sc_dim)
            sc_scale[i] -= 1
            scs = vac.make_supercells_with_defects(sc_scale)
    print len(scs)

    interdir = mpid
    blk_str_sites = set(scs[0].sites)
    for i in range(1,len(scs)):
        sc = scs[i]
        vac_str_sites = set(sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = list(vac_sites)[0]
        site_mult = vac.get_defectsite_multiplicity(i-1)
        vac_site_specie = vac_site.specie
        vac_specie = vac_site.specie.symbol

        # Solute substitution defect generation at all vacancy sites
        struct_species = scs[0].types_of_specie
        solute_struct = sc.copy()
        solute_struct.append(solute, vac_site.frac_coords)

        incar = mpvis.get_incar(solute_struct)
        incar.update(def_vasp_incar_param)
        poscar = mpvis.get_poscar(solute_struct)
        kpoints = Kpoints.automatic_density(solute_struct,kpoint_den)
        if ptcr_flag:
            potcar = mpvis.get_potcar(solute_struct)

        sub_def_dir ='solute_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                str(i), site_mult, vac_specie, solute)
        fin_dir = os.path.join(interdir,sub_def_dir)
        try:
            os.makedirs(fin_dir)
        except:
            pass
        poscar.write_file(os.path.join(fin_dir,'POSCAR'))
        incar.write_file(os.path.join(fin_dir,'INCAR'))
        kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        if ptcr_flag:
            potcar.write_file(os.path.join(fin_dir,'POTCAR'))

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
    print args
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

    parser.add_argument("--solute", help="Solute Element")

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
    print args
    substitute_def_struct_gen(args.mpid,args.solute,args.mapi_key,args.cellmax)

if __name__ == '__main__':
    im_vac_antisite_def_struct_gen()
    #im_sol_sub_def_struct_gen()

