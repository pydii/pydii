#!/usr/bin/env python
"""
This file reads arguments and generate vacancy and antisite structures
in intermetallics.
"""

__author__ = "Bharat Medasani, Enze Chen"
__data__  = "Aug 10, 2020"

import os
from argparse import ArgumentParser

from pymatgen.io.vasp.sets import MPMetalRelaxSet
from pymatgen.ext.matproj import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import Vacancy
from pymatgen.io.vasp.inputs import Kpoints


def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1./3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(mpid, mapi_key, cellmax, struct_file=None):
    if not mpid and not struct_file:
        print ("============\nERROR: Provide an mpid\n============")
        return

    # Get primitive structure from the Materials Project DB
    if not struct_file:
        if not mapi_key:
            with MPRester() as mp:
                struct = mp.get_structure_by_material_id(mpid)
        else:
            with MPRester(mapi_key) as mp:
                struct = mp.get_structure_by_material_id(mpid)
    else:
        struct = Structure.from_file(struct_file)

    sga = SpacegroupAnalyzer(struct)
    prim_struct = sga.find_primitive()
    #prim_struct_sites = len(prim_struct.sites)
    #conv_struct = sga.get_conventional_standard_structure()
    #conv_struct_sites = len(conv_struct.sites)
    #conv_prim_ratio = int(conv_struct_sites / prim_struct_sites)

    # Default VASP settings
    def_vasp_incar_param = {'ISIF':2, 'EDIFF':1e-6, 'EDIFFG':0.001,}
    kpoint_den = 15000

    # Create bulk structure and associated VASP files
    sc_scale = get_sc_scale(inp_struct=prim_struct, final_site_no=cellmax)
    blk_sc = prim_struct.copy()
    blk_sc.make_supercell(scaling_matrix=sc_scale)
    site_no = blk_sc.num_sites

    # Rescale if needed
    while site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        blk_sc = prim_struct.copy()
        blk_sc.make_supercell(scaling_matrix=sc_scale)
        site_no = blk_sc.num_sites
    
    blk_str_sites = set(blk_sc.sites)
    custom_kpoints = Kpoints.automatic_density(blk_sc, kppa=kpoint_den)
    mpvis = MPMetalRelaxSet(blk_sc, user_incar_settings=def_vasp_incar_param,
                            user_kpoints_settings=custom_kpoints)

    if mpid:
        root_fldr = mpid
    else:
        root_fldr = struct.composition.reduced_formula

    fin_dir = os.path.join(root_fldr, 'bulk')
    mpvis.write_input(fin_dir)
    if not mpid:    # write the input structure if mpid is not used
        struct.to(fmt='poscar', filename=os.path.join(fin_dir, 'POSCAR.uc'))

    # Create each defect structure and associated VASP files
    # First find all unique defect sites
    periodic_struct = sga.get_symmetrized_structure()
    unique_sites = list(set([periodic_struct.find_equivalent_sites(site)[0] \
                             for site in periodic_struct.sites]))
    temp_struct = Structure.from_sites(sorted(unique_sites))
    prim_struct2 = SpacegroupAnalyzer(temp_struct).find_primitive()
    prim_struct2.lattice = prim_struct.lattice  # a little hacky
    for i, site in enumerate(prim_struct2.sites):
        vac = Vacancy(structure=prim_struct, defect_site=site)
        vac_sc = vac.generate_defect_structure(supercell=sc_scale)

        # Get vacancy site information
        vac_str_sites = set(vac_sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = next(iter(vac_sites))
        site_mult = vac.get_multiplicity()
        vac_site_specie = vac_site.specie
        vac_symbol = vac_site_specie.symbol

        custom_kpoints = Kpoints.automatic_density(vac_sc, kppa=kpoint_den)
        mpvis = MPMetalRelaxSet(vac_sc,
                                user_incar_settings=def_vasp_incar_param,
                                user_kpoints_settings=custom_kpoints)
        vac_dir = 'vacancy_{}_mult-{}_sitespecie-{}'.format(
                    str(i+1), site_mult, vac_symbol)
        fin_dir = os.path.join(root_fldr, vac_dir)
        mpvis.write_input(fin_dir)

        # Antisites generation at the vacancy site
        struct_species = blk_sc.species
        for specie in set(struct_species) - set([vac_site_specie]):
            specie_symbol = specie.symbol
            anti_sc = vac_sc.copy()
            anti_sc.append(specie, vac_site.frac_coords)
            mpvis = MPMetalRelaxSet(anti_sc,
                                    user_incar_settings=def_vasp_incar_param,
                                    user_kpoints_settings=custom_kpoints)
            anti_dir = 'antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                        str(i+1), site_mult, vac_symbol, specie_symbol)
            fin_dir = os.path.join(root_fldr, anti_dir)
            mpvis.write_input(fin_dir)


def substitute_def_struct_gen(mpid, solute, mapi_key, cellmax,
                              struct_file=None):
    if not mpid and not struct_file:
        print ("============\nERROR: Provide an mpid\n============")
        return
    if not solute:
        print ("============\nERROR: Provide solute atom\n============")
        return

    # Get primitive structure from the Materials Project DB
    if not struct_file:
        if not mapi_key:
            with MPRester() as mp:
                struct = mp.get_structure_by_material_id(mpid)
        else:
            with MPRester(mapi_key) as mp:
                struct = mp.get_structure_by_material_id(mpid)
    else:
        struct = Structure.from_file(struct_file)

    if mpid:
        root_fldr = mpid
    else:
        root_fldr = struct.composition.reduced_formula

    sga = SpacegroupAnalyzer(struct)
    prim_struct = sga.find_primitive()
    #prim_struct_sites = len(prim_struct.sites)
    #conv_struct = sga.get_conventional_standard_structure()
    #conv_struct_sites = len(conv_struct.sites)
    #conv_prim_ratio = int(conv_struct_sites / prim_struct_sites)

    # Default VASP settings
    def_vasp_incar_param = {'ISIF':2, 'EDIFF':1e-6, 'EDIFFG':0.001,}
    kpoint_den = 15000

    # Create each substitutional defect structure and associated VASP files
    sc_scale = get_sc_scale(inp_struct=prim_struct, final_site_no=cellmax)
    blk_sc = prim_struct.copy()
    blk_sc.make_supercell(scaling_matrix=sc_scale)
    site_no = blk_sc.num_sites

    # Rescale if needed
    while site_no > cellmax:
        max_sc_dim = max(sc_scale)
        i = sc_scale.index(max_sc_dim)
        sc_scale[i] -= 1
        blk_sc = prim_struct.copy()
        blk_sc.make_supercell(scaling_matrix=sc_scale)
        site_no = blk_sc.num_sites

    # Create solute structures at vacancy sites
    # First find all unique defect sites
    blk_str_sites = set(blk_sc.sites)
    periodic_struct = sga.get_symmetrized_structure()
    unique_sites = list(set([periodic_struct.find_equivalent_sites(site)[0] \
                             for site in periodic_struct.sites]))
    temp_struct = Structure.from_sites(sorted(unique_sites))
    prim_struct2 = SpacegroupAnalyzer(temp_struct).find_primitive()
    prim_struct2.lattice = prim_struct.lattice  # a little hacky
    for i, site in enumerate(prim_struct2.sites):
        vac = Vacancy(structure=prim_struct, defect_site=site)
        vac_sc = vac.generate_defect_structure(supercell=sc_scale)

        # Get vacancy site information
        vac_str_sites = set(vac_sc.sites)
        vac_sites = blk_str_sites - vac_str_sites
        vac_site = next(iter(vac_sites))
        vac_specie = vac_site.specie.symbol
        site_mult = vac.get_multiplicity()

        # Solute substitution defect generation at the vacancy site
        solute_struct = vac_sc.copy()
        solute_struct.append(solute, vac_site.frac_coords)
        custom_kpoints = Kpoints.automatic_density(solute_struct,
                                                   kppa=kpoint_den)
        mpvis = MPMetalRelaxSet(solute_struct,
                                user_incar_settings=def_vasp_incar_param,
                                user_kpoints_settings=custom_kpoints)

        # Generate VASP directory
        sub_def_dir ='solute_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                str(i+1), site_mult, vac_specie, solute)
        fin_dir = os.path.join(root_fldr, sub_def_dir)
        mpvis.write_input(fin_dir)


def im_vac_antisite_def_struct_gen():
    m_description = 'Command to generate vacancy and antisite defect ' \
                    'structures for intermetallics.'

    parser = ArgumentParser(description=m_description)

    parser.add_argument(
            "--mpid",
            default=None,
            type=str.lower,
            help="Materials Project id of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument(
            "--struct",
            default=None,
            type=str,
            help="Filename of the intermetallic structure.")

    parser.add_argument(
            "--mapi_key",
            default=None,
            help="Your Materials Project REST API key.\n" \
                 "For more info, please refer to " \
                 "www.materialsproject.org/open")

    parser.add_argument(
            "--cellmax",
            type=int,
            default=128,
            help="Maximum number of atoms in supercell.\n" \
                 "The default is 128\n" \
                 "Keep in mind the number of atoms in the supercell" \
                 "may vary from the provided number including the default.")

    args = parser.parse_args()
    vac_antisite_def_struct_gen(args.mpid, args.mapi_key, args.cellmax,
                                struct_file=args.struct)


def im_sol_sub_def_struct_gen():
    m_description = 'Command to generate solute substitution defect ' \
                    'structures for intermetallics.'

    parser = ArgumentParser(description=m_description)

    parser.add_argument(
            "--mpid",
            default=None,
            type=str.lower,
            help="Materials Project id of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument(
            "--struct",
            default=None,
            type=str,
            help="Filename of the intermetallic structure." \
                 "Supported file types include CIF, POSCAR/CONTCAR," \
                 "CHGCAR, LOCPOT, vasprun.xml, CSSR, Netcdf, and pymatgen JSONs.")

    parser.add_argument(
            "--solute",
            type=str,
            help="Solute Element")

    parser.add_argument(
            "--mapi_key",
            default=None,
            help="Your Materials Project REST API key.\n" \
                 "For more info, please refer to " \
                 "www.materialsproject.org/open")

    parser.add_argument(
            "--cellmax",
            type=int,
            default=128,
            help="Maximum number of atoms in supercell.\n" \
                 "The default is 128\n" \
                 "Keep in mind the number of atoms in the supercell" \
                 "may vary from the provided number including the default.")

    args = parser.parse_args()
    substitute_def_struct_gen(args.mpid, args.solute, args.mapi_key,
                              args.cellmax, struct_file=args.struct)


if __name__ == '__main__':
    # im_vac_antisite_def_struct_gen()
    im_sol_sub_def_struct_gen()
