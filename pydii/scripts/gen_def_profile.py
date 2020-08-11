#!/usr/bin/env python
"""
PyDII file for generating defect concentrations
"""
__author__ = "Bharat Medasani, Enze Chen"
__date__ = "Aug 10, 2020"

import os
from argparse import ArgumentParser

from monty.serialization import loadfn
from monty.json import MontyDecoder
from pydii.dilute_solution_model import compute_defect_density, solute_defect_density


def get_def_profile(mpid, T,  file_nm):
    raw_energy_dict = loadfn(file_nm, cls=MontyDecoder)

    e0 = raw_energy_dict['e0']
    struct = raw_energy_dict['structure']
    vacs = raw_energy_dict['vacancies']
    antisites = raw_energy_dict['antisites']
    vacs.sort(key=lambda entry: entry['site_index'])
    antisites.sort(key=lambda entry: entry['site_index'])
    for vac_def in vacs:
        if not vac_def:
            print('All vacancy defect energies not present.')
            continue
    for antisite_def in antisites:
        if not antisite_def:
            print('All antisite defect energies not present.')
            continue

    try:
        def_conc, def_en, mu = compute_defect_density(struct, e0, vacs,
                antisites, T, plot_style='gnuplot')
        return def_conc, def_en, mu
    except:
        raise


def get_solute_def_profile(mpid, solute, solute_conc, T, def_file, sol_file, 
                           trial_chem_pot):
    raw_energy_dict = loadfn(def_file, cls=MontyDecoder)
    sol_raw_energy_dict = loadfn(sol_file, cls=MontyDecoder)

    e0 = raw_energy_dict['e0']
    struct = raw_energy_dict['structure']
    vacs = raw_energy_dict['vacancies']
    antisites = raw_energy_dict['antisites']
    solutes = sol_raw_energy_dict['solutes']
    for vac_def in vacs:
        if not vac_def:
            print('All vacancy defect energies not present.')
            continue
    for antisite_def in antisites:
        if not antisite_def:
            print('All antisite defect energies not present.')
            continue
    for solute_def in solutes:
        if not solute_def:
            print('All solute defect energies not present.')
            continue

    try:
        def_conc = solute_defect_density(struct, e0, vacs,
                antisites, solutes, solute_concen=solute_conc, T=T,
                trial_chem_pot=trial_chem_pot, plot_style="gnuplot")
        return  def_conc
    except:
        raise


def im_vac_antisite_def_profile():
    m_description = 'Command to generate vacancy and antisite defect ' \
                    'concentration for intermetallics from the raw defect ' \
                    'energies.'

    parser = ArgumentParser(description=m_description)

    parser.add_argument(
            "--mpid",
            default=None,
            type=str.lower,
            help="Materials Project ID of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument(
            "--formula",
            default=None,
            type=str,
            help="Reduced formula of the intermetallic compound.")

    parser.add_argument('-T', "--temp", type=float, default=1000,
            help="Temperature in Kelvin")

    parser.add_argument("--file",
            default = None,
            help = "The default file is 'mpid'+'_raw_defect_energy.json'.\n" \
                   "If the file is named differently supply it.")

    args = parser.parse_args()

    if not args.mpid and not args.formula:
        print ('===========\nERROR: mpid is not given.\n===========')
        return

    if not args.file:
        if args.mpid:
            file_name = args.mpid + '_raw_defect_energy.json'
        elif args.formula:
            file_name = args.formula + '_raw_defect_energy.json'
    else:
        file_name = args.file

    conc_dat, en_dat, mu_dat = get_def_profile(args.mpid, args.temp, file_name)

    if args.mpid:
        prefix = args.mpid
    elif args.formula:
        prefix = args.formula

    if conc_dat:
        fl_nm = prefix + '_def_concentration.dat'
        with open(fl_nm, 'w') as fp:
            for row in conc_dat:
                print(row, file=fp)
        fl_nm = prefix + '_def_energy.dat'
        with open(fl_nm, 'w') as fp:
            for row in en_dat:
                print(row, file=fp)
        fl_nm = prefix + '_chem_pot.dat'
        with open(fl_nm, 'w') as fp:
            for row in mu_dat:
                print(row, file=fp)


def im_sol_sub_def_profile():
    m_description = 'Command to generate solute defect site preference ' \
                    'in an intermetallics from the raw defect energies.'

    parser = ArgumentParser(description=m_description)

    parser.add_argument(
            "--mpid",
            default=None,
            type=str.lower,
            help="Materials Project id of the intermetallic structure.\n" \
                 "For more info on Materials Project, please refer to " \
                 "www.materialsproject.org")

    parser.add_argument(
            "--formula",
            default=None,
            type=str,
            help="Reduced formula of the intermetallic compound.")

    parser.add_argument(
            "--solute",
            type=str,
            help="Solute Element")

    parser.add_argument(
            "--sol_conc",
            type=float,
            default=1.0,
            help="Solute Concentration in %. Default is 1%")

    parser.add_argument(
            "-T", "--temp",
            type=float,
            default=1000.0,
            help="Temperature in Kelvin")

    parser.add_argument("--trial_mu_file",
            default=None,
            help="Trial chemical potential in dict format stored in file.")

    args = parser.parse_args()

    if not args.mpid and not args.formula:
        print ('===========\nERROR: mpid is not given.\n===========')
        return
    if not args.solute:
        print ('===========\nERROR: Solute atom is not given.\n===========')
        return

    if args.mpid:
        prefix = args.mpid
    elif args.formula:
        prefix = args.formula

    def_file = prefix + '_raw_defect_energy.json'
    sol_file = prefix + '_solute-' + args.solute + '_raw_defect_energy.json'
    sol_conc = args.sol_conc / 100.0 # Convert from percentage
    if not os.path.exists(def_file):
        print ('===========\nERROR: Defect file not found.\n===========')
        return
    if not os.path.exists(sol_file):
        print ('===========\nERROR: Solute file not found.\n===========')
        return
    if args.trial_mu_file:
        trial_chem_pot = loadfn(args.trial_mu_file, cls=MontyDecoder)
    else:
        trial_chem_pot = None

    pt_def_conc = get_solute_def_profile(
            args.mpid, args.solute, sol_conc, args.temp, def_file,
            sol_file, trial_chem_pot=trial_chem_pot)

    if pt_def_conc:
        fl_nm = prefix + '_solute-' + args.solute + '_def_concentration.dat'
        with open(fl_nm, 'w') as fp:
            for row in pt_def_conc:
                print(row, file=fp)


if __name__ == '__main__':
    im_vac_antisite_def_profile()
    # im_sol_sub_def_profile()
