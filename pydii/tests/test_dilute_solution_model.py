# coding: utf-8
import unittest
import json
import os
import sympy

from monty.json import MontyDecoder
from pydii.dilute_solution_model import *

# Relative path to the directory containing this script
print('Dir: {}'.format(os.path.dirname(__file__)))
test_dir = os.path.join(os.path.dirname(__file__))

with open(os.path.join(test_dir, 'mp1048_defect_formation_energies.json')) as fp:
    formation_energy_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'mp1048_raw_defect_energies.json')) as fp:
    raw_energy_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'mp1487_raw_defect_energies.json')) as fp:
    mp1487_raw_energy_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'mp-2554_raw_defect_energies.json')) as fp:
    mp2554_raw_energy_dict = json.load(fp, cls=MontyDecoder)


@unittest.skipIf(not sympy, "sympy not present.")
class DiluteSolutionModelTest(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        self.e0 = raw_energy_dict['bulk_energy']
        self.asites = raw_energy_dict['antisites']
        self.vac = raw_energy_dict['vacancies']
        self.struct = raw_energy_dict['structure']
        self.T = 600

    def test_formation_energies_without_chem_pot(self):
        """
        Should generate formation energies without input chempot
        """
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_formation_energies_with_chem_pot(self):
        self.trial_mu = formation_energy_dict[str(self.T)]['chemical_potential']
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_plot_data_without_chem_pot(self):
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key,value in conc_data.items():
            self.assertIsNotNone(value)
        for key,value in mu_data.items():
            self.assertIsNotNone(value)
        for key,value in en_data.items():
            self.assertIsNotNone(value)

    def test_plot_data_with_chem_pot(self):
        self.trial_mu = formation_energy_dict[str(self.T)]['chemical_potential']
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key,value in conc_data.items():
            self.assertIsNotNone(value)
        for key,value in mu_data.items():
            self.assertIsNotNone(value)
        for key,value in en_data.items():
            self.assertIsNotNone(value)
        #print(plot_data['y'])


@unittest.skipIf(not sympy, "sympy not present.")
class DiluteSolutionModelSiteMultiplicityTest(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        print('Keys: {}'.format(mp2554_raw_energy_dict.keys()))
        self.e0 = mp2554_raw_energy_dict['e0']
        self.asites = mp2554_raw_energy_dict['antisites']
        self.vac = mp2554_raw_energy_dict['vacancies']
        self.struct = mp2554_raw_energy_dict['structure']
        self.T = 1000

    def test_formation_energies_without_chem_pot(self):
        """
        Should generate formation energies without input chempot
        """
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_plot_data_without_chem_pot(self):
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key,value in conc_data.items():
            self.assertIsNotNone(value)
        for key,value in mu_data.items():
            self.assertIsNotNone(value)
        for key,value in en_data.items():
            self.assertIsNotNone(value)


@unittest.skipIf(not sympy, "sympy not present.")
class DiluteSolutionModel1487Test(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        self.e0 = mp1487_raw_energy_dict['bulk_energy']
        self.asites = mp1487_raw_energy_dict['antisites']
        self.vac = mp1487_raw_energy_dict['vacancies']
        self.struct = mp1487_raw_energy_dict['structure']
        self.T = 1000

    def test_formation_energies_without_chem_pot(self):
        """
        Should generate formation energies without input chempot
        """
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_plot_data_without_chem_pot(self):
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T, generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key,value in conc_data.items():
            self.assertIsNotNone(value)
        for key,value in mu_data.items():
            self.assertIsNotNone(value)
        for key,value in en_data.items():
            self.assertIsNotNone(value)


@unittest.skipIf(not sympy, "sympy not present.")
class SoluteSiteFinderTest(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        self.e0 = mp1487_raw_energy_dict['bulk_energy']
        self.asites = mp1487_raw_energy_dict['antisites']
        self.vac = mp1487_raw_energy_dict['vacancies']
        self.solutes = mp1487_raw_energy_dict['solutes']
        self.struct = mp1487_raw_energy_dict['structure']
        self.T = 1000

    def test_plot_data_without_chem_pot(self):
        plot_data = solute_site_preference_finder(
            self.struct, self.e0, self.T, self.vac, self.asites, self.solutes,
            solute_concen=0.01)
        print(plot_data.keys())
        self.assertIsNotNone(plot_data)

    def still_wait_plot_data_with_chem_pot(self):
        plot_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='plot')
        self.assertIsNotNone(plot_data)
        for key,value in plot_data.items():
            self.assertIsNotNone(value)
        print(plot_data['y'])


if __name__ == "__main__":
    unittest.main()
