# coding: utf-8
from __future__ import unicode_literals
from __future__ import division

"""
Evaluate the defect concentration based on composition, temperature,
and defect energies using the "Dilute Solution Model"
Reference:
    C. Woodward, M. Asta, G. Kresse and J. Hafner.
    "Density of constitutional and thermal point defects in L12 Al3Sc,"
    Physical Review B, 63, 094103, 2001,
"""

__author__ = 'Bharat Medasani, Enze Chen'
__version__ = "0.3"
__maintainer__ = "Enze Chen"
__email__ = "chenze@berkeley.edu"
__status__ = "Beta"
__date__ = "2020/06/25"

import sys
import copy
import numpy as np
from six.moves import zip
from monty.fractions import gcd

try:
    import sympy
    from sympy import Symbol, Integer, Float, Matrix, exp, Eq
except:
    raise ImportError('Sympy module cannot be found.')

# physical constants
k_B = 8.6173324e-5  # eV/K

# Check the inputs
def check_input(def_list):
    """Check that all defect structures in list are not None.

    Args:
        def_list (list): List of defect structures

    Returns:
        True if all defect structures are present.
    """
    for defect in def_list:
        if not defect:
            return False
    return True


def dilute_solution_model(structure, e0, vac_defs, antisite_defs, T, trial_chem_pot = None, generate = 'plot'):
    """Compute the defect densities using the dilute solution model.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0 (float): The total energy of the undefected system.
            This is E0 from a VASP calculation.
        vac_defs (list): Vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from the
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (Fortran index).
        antisite_defs (list): Antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from the pymatgen.analysis.defects.point_defects.Vacancy class.
        T (float): Temperature in Kelvin
        trial_chem_pot (dict): Trial chemical potentials to speedup
            the plot generation. Format is {el1:mu1,...}. Optional.
        generate (string): Options are 'plot' or 'energy'
            Chemical potentials are also returned with the 'energy' option.
            If the 'energy' option is not chosen, a plot is generated.

    Returns:
        If generate='plot', the plot data is generated and returned in
        HighCharts format.
        If generate='energy', defect formation enthalpies and chemical
        potentials are returned.
    """

    if not check_input(vac_defs):
        raise ValueError('Vacancy energy is not defined')
    if not check_input(antisite_defs):
        raise ValueError('Antisite energy is not defined')

    formation_energies = {}
    formation_energies['vacancies'] = copy.deepcopy(vac_defs)
    formation_energies['antisites'] = copy.deepcopy(antisite_defs)
    for vac in formation_energies['vacancies']:
        del vac['energy']
    for asite in formation_energies['antisites']:
        del asite['energy']

    # Setup the system
    site_species = [vac_def['site_specie'] for vac_def in vac_defs]
    multiplicity = [vac_def['site_multiplicity'] for vac_def in vac_defs]
    m = len(set(site_species))  # distinct species
    n = len(vac_defs)           # inequivalent sites
    # print(f'Site species: {site_species}')
    # print(f'Multiplicity: {multiplicity}')
    # print(f'Total number of distinct species: {m}')
    # print(f'Total number of inequivalent sites: {n}')

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val / comm_div for val in multiplicity]
    e0 = e0 / comm_div
    T = Float(T)
    # print(f'Energy: {e0}')

    c0 = np.diag(np.ones(n))
    mu = [Symbol(f'mu{i}') for i in range(m)]
    # print(f'Mu list: {mu}')

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []       # Contains hash for site->mu map  e.g. [Al, Al, V]
    site_specie_set = set()     # e.g. {V, Al}
    for i in range(n):
        site_specie  = site_species[i]
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []            # e.g. [mu0, mu0, mu1] where mu0->Al, and mu1->V
    for i in range(n):
        site_specie  = site_species[i]
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []  # e.g. [(0, 2), (2, 3)] for Al and V
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m-1:
            high_ind = site_species.index(specie_order[i+1])
        else:
            high_ind = n
        specie_site_index_map.append((low_ind, high_ind))
    # print(f'Specie order: {specie_order}')
    # print(f'Site specie set: {site_specie_set}')
    # print(f'Site mu map: {site_mu_map}')
    # print(f'Specie site index map: {specie_site_index_map}')


    """
    dC: delta concentration matrix:
    dC[i,j,k]: Concentration change of atom i, due to presence of atom
               j on lattice site k.
    Special case is [i,i,i] which is considered as vacancy.
    Few cases: (1) dC[i,i,i] = -1 due to being vacancy special case
               (2) dC[k,k,i] = +1 due to increment in k at i lattice if i
                                  lattice type is of different element
               (3) dC[i,k,i] = -1 due to decrement of ith type atom due to
                                  presence of kth type atom on ith sublattice
                                  and kth type atom specie is different from
                                  ith sublattice atom specie
               (4) dC[i,k,k] =  0 due to no effect on ith type atom
               (5) dC[i,j,k] =  0 if i != j != k
    """
    dC = np.zeros((n, n, n), dtype=np.int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # case (2)
                if i == j and site_species[j] != site_species[k]:
                    dC[i, j, k] = 1
                # case (1) and (3)
                if i == k:
                    dC[i, j, k] = -1
                # case (4)
                if i != j and site_species[j] == site_species[k]:
                    dC[i, j, k] = 0
    # print(f'Delta concentration matrix:\n{dC}')

    # dE matrix: Flip energies (or raw defect energies)
    els = site_species.copy()
    dE = [[0 for _ in range(n)] for _ in range(n)]
    for j in range(n):
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]['energy']
            else:
                sub_specie = vac_defs[i]['site_specie']
                site_specie = vac_defs[j]['site_specie']
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def['site_index']) == j+1 and \
                           sub_specie == as_def['substitution_specie']:
                            dE[i][j] = as_def['energy']
                            break
    dE = np.array(dE)
    # print(f'Delta energies matrix:\n{dE}')

    # Initialization for concentrations
    # c(i, p) == presence of ith type atom on pth type site
    c = sympy.zeros(n, n)
    for i in range(n):
        for p in range(n):
            c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) \
                              for j in range(n)])
                flip = Integer(dC[i, epi, p]) * \
                       exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    c[i, p] += flip
    total_c = []
    for ind in specie_site_index_map:
        val = 0
        for i in range(*ind):
            val += sum([c[i, j] * multiplicity[j] for j in range(n)])
        total_c.append(val)
    c_ratio = [total_c[-1] / total_c[i] for i in range(m)]
    # print(f'Concentration matrix:\n{c}')
    # print(f'c_ratio:\n{c_ratio}')

    # Expression for Omega, the Grand Potential
    omega1 = e0 - sum([mu[site_mu_map[i]] * c0[i, i] * multiplicity[i] \
                       for i in range(n)])
    omega2 = 0
    used_dEs = []
    for p_r in range(n):
        for epi in range(n):
            sum_mu = sum([mu[site_mu_map[j]] * \
                          Float(dC[j, epi, p_r]) for j in range(n)])
            if p_r != epi and site_mu_map[p_r] == site_mu_map[epi]:
                continue
            if dE[epi, p_r] not in used_dEs:
                omega2 += k_B * T * multiplicity[p_r] * \
                          exp(-(dE[epi, p_r] - sum_mu) / (k_B * T))
                used_dEs.append(dE[epi, p_r])
    omega = omega1 - omega2
    # print(f'Omega: {omega}')

    # Compute composition range
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = sum(multiplicity[li:hi])/sum(multiplicity) * 100 - 1
    comp1_max = sum(multiplicity[li:hi])/sum(multiplicity) * 100 + 1
    delta = float(comp1_max - comp1_min) / 120.0    # dc = 1/60
    yvals = []
    for comp1 in np.arange(comp1_min, comp1_max+delta, delta):
        comp2 = 100 - comp1
        y = comp2 / comp1
        yvals.append(y)
    # print(f'Composition [min, max]: [{comp1_min}, {comp1_max}]')

    def reduce_mu():
        omega = [e0 - sum([mu[site_mu_map[i]] * sum(c0[i, :]) \
                           for i in range(n)])]
        x = sympy.solve(omega)
        return x

    def compute_mus_by_search():
        # Compute trial mu
        mu_red = reduce_mu()
        specie_concen = [sum(multiplicity[ind[0]:ind[1]]) \
                         for ind in specie_site_index_map]
        y_vect = [specie_concen[-1] / specie_concen[i] for i in range(m)]
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m-1)]
        vector_func.append(omega)
        min_diff = 1e10
        mu_vals = None
        c_val = None
        m1_min = -20.0
        if e0 > 0:
            m1_max = 10            # Search space needs to be modified
        else:
            m1_max = 0
        for m1 in np.arange(m1_min, m1_max, 0.01):
            m0 = mu_red[mu[0]].subs(mu[-1], m1)

            try:
                x = sympy.nsolve(vector_func, mu, [m0, m1])
            except:
                continue

            c_val = c.subs(dict(zip(mu,x)))
            specie_concen = []
            for ind in specie_site_index_map:
                specie_concen.append(sum([sum(c_val[i, :]) \
                                          for i in range(*ind)]))
            y_comp = [specie_concen[-1] / specie_concen[i] for i in range(m)]
            diff = np.sqrt(sum([pow(abs(y_comp[i] - y_vect[i]), 2) \
                                for i in range(m)]))
            if diff < min_diff:
                min_diff = diff
                mu_vals = x

        if mu_vals:
            mu_vals = [float(mu_val) for mu_val in mu_vals]
        else:
            raise ValueError()
        return mu_vals

    def compute_def_formation_energies():
        i = 0
        for vac_def in vac_defs:
            site_specie = vac_def['site_specie']
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def['energy']
            formation_energy = uncor_energy + mu_vals[ind]
            formation_energies['vacancies'][i]['formation_energy'] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1] - indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del-1:
                label = '$V_{' + site_specie + '}$'
            else:
                label = '$V_{' + site_specie + '_' + str(cur_ind) + '}$'
            formation_energies['vacancies'][i]['label'] = label
            i += 1
        i = 0
        for as_def in antisite_defs:
            site_specie = as_def['site_specie']
            sub_specie = as_def['substitution_specie']
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def['energy']
            formation_energy = uncor_energy + mu_vals[ind1] - mu_vals[ind2]
            formation_energies['antisites'][i]['formation_energy'] = formation_energy
            specie_ind = site_mu_map[i]
            indices = specie_site_index_map[specie_ind]
            specie_ind_del = indices[1] - indices[0]
            cur_ind = i - indices[0] + 1
            if not specie_ind_del-1:
                label = '$' + sub_specie + '_{' + site_specie + '}$'
            else:
                label = '$' + sub_specie + '_{' + site_specie + '_' + str(cur_ind) + '}$'
            formation_energies['antisites'][i]['label'] = label
            i += 1
        return formation_energies

    # If generate option is energy compute effective formation energies
    # at ideal stoichiometry and return the formation energies and chem pot.
    if generate == 'energy':
        print('Using energy generate method.')
        if not trial_chem_pot:
            mu_vals = compute_mus_by_search()
        else:
            try:
                mu_vals = [trial_chem_pot[element] for element in specie_order]
            except:
                mu_vals = compute_mus()

        formation_energies = compute_def_formation_energies()
        mu_dict = dict(zip(specie_order, mu_vals))
        return formation_energies, mu_dict

    if not trial_chem_pot:
        print('No trial chemical potential is given.')
        # Try computing mus by assuming one of the defects is dominant at 0.01
        # concen.  First vacancy is tried and then antisite

        # Generate trial mus assuming vacancy as dominant defect
        # for specie-0 at lower yval
        li0 = specie_site_index_map[0][0]
        hi0 = specie_site_index_map[0][1]
        li1 = specie_site_index_map[1][0]
        hi1 = specie_site_index_map[1][1]
        spec_mult = [sum(multiplicity[li0:hi0]), sum(multiplicity[li1:hi1])]
        ln_def_conc = 4.60517       # -ln(0.01)
        for i in range(li0, hi0):
            vac_flip_en =  vac_defs[i]['energy']
            mu_vals = [ln_def_conc * k_B * T - vac_flip_en]
            mu_vals.append((e0 - spec_mult[0] * mu_vals[0]) / spec_mult[1])
            comp_ratio = yvals[0]

            # Test if the trial mus are good
            vector_func = [comp_ratio - c_ratio[0]]
            vector_func.append(omega)
            try:
                mu_vals = sympy.nsolve(vector_func, mu, mu_vals)
                if mu_vals:
                    mu_vals = [float(mu_val) for mu_val in mu_vals]
                break
            except:     # Go for antisite as dominant defect
                mu_gs = [Symbol(f'mu_gs{j}') for j in range(m)]
                eqs = [mu_gs[0] - mu_gs[1] - (ln_def_conc * k_B * T - \
                       antisite_defs[i]['energy'])]
                eqs.append(spec_mult[0] * mu_gs[0] + spec_mult[1] * mu_gs[1] - e0)
                x = sympy.solve(eqs, mu_gs)
                mu_vals = [x[key] for key in sorted(x.keys(), key=lambda inp: inp.name)]

                try:
                    mu_vals = sympy.nsolve(vector_func, mu, mu_vals)
                    if mu_vals:
                        mu_vals = [float(mu_val) for mu_val in mu_vals]
                    break
                except:     # Go to the default option (search the space)
                    pass
        else:
            mu_vals = compute_mus_by_search()
    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except:
            mu_vals = compute_mus_by_search()
    print(f'Trial mu_vals: {mu_vals}')

    # COMPUTE MU's FOR ALL COMPOSITION RATIOS in the range
    # +/- 1% from the stoichiometry. MAIN SOLVER!!
    result = {}
    failed_y, failed_i = [], []
    for i, y in enumerate(yvals):
        vector_func = [y - c_ratio[0]]
        vector_func.append(omega)
        try:
            x = sympy.nsolve(vector_func, mu, mu_vals)
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except:
            failed_y.append(y)
            failed_i.append(i)
            continue
        result[y] = mu_vals
        x = None
    print(f'Failed {len(failed_i)} times.')

    def get_next_mu_val(i):
        if i >= len(yvals):
            return None

        y =  yvals[i+1]
        x = result.get(y, None)
        if x:
            mu_vals = [float(mu_val) for mu_val in x]
            return mu_vals
        else:
            return get_next_mu_val(i+1)

    def get_prev_mu_val(i):
        if i <= 0:
            return None

        y =  yvals[i-1]
        x = result.get(y, None)
        if x:
            mu_vals = [float(mu_val) for mu_val in x]
            return mu_vals
        else:
            return get_next_mu_val(i-1)

    # Try to get better trial mus for failed cases
    for j in range(len(failed_y)):
        i = failed_i[j]

        prev_mu_val = get_prev_mu_val(i)
        if not prev_mu_val:
            continue
        next_mu_val = get_next_mu_val(i)
        if not next_mu_val:
            continue

        y = failed_y[j]
        vector_func = [y - c_ratio[0]]
        vector_func.append(omega)
        trial_mu = list(map(lambda x: float(sum(x)) / len(x), \
                zip(prev_mu_val, next_mu_val)))
        try:
            x = sympy.nsolve(vector_func, mu, trial_mu)
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except:
            continue
        result[y] = mu_vals
        x = None

    if len(result.keys()) < len(yvals)/2:
        raise ValueError('Not sufficient data.')

    res = []
    new_mu_dict = {}
    # COMPUTE THE CONCENTRATION for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu, mu_val))) \
                       for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu, mu_val)))
        res1 = [float(total_c_val[0] / sum(total_c_val))]
        # Concentration of first element/over total concentration
        new_mu_dict[res1[0]] = mu_val
        # sum_c0 = sum([c0[i, i] for i in range(n)])
        for i in range(n):
            for j in range(n):
                if i == j:              # Vacancy
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]] + dE[i, i])/(k_B * T)))
                    res1.append(vac_conc)
                    if vac_conc > 1.0:
                        print('Warning! Vacancy concentration is > 1.')
                else:                   # Antisite
                    anti_conc = float(c_val[i, j] / c0[j, j])
                    res1.append(anti_conc)
                    if anti_conc > 1.0:
                        print('Warning! Antisite concentration is > 1.')
        res.append(res1)

    res = np.array(res)
    dtype = [(str('x'), np.float64)] + [(f'y{i}{j}', np.float64) \
                                        for i in range(n) for j in range(n)]
    res1 = np.sort(res.view(dtype), order=[str('x')], axis=0)

    conc_data = {}
    # Because all the plots have identical x-points, we store them in a single array
    conc_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    conc_data['x_label'] = els[0] + " mole fraction"
    conc_data['y_label'] = "Point defect concentration"
    conc = [[[] for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            # n x n array where each element is a list of concentrations
            conc[i][j] = [dat[0][i * n + j + 1] for dat in res1]

    # Get the vacancies concentration for each site
    y_data = []
    for i in range(n):
        data = conc[i][i]
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = f'{vac_string}{site_specie}}}$'
        else:
            label = f'{vac_string}{site_specie}_{cur_ind}}}$'
        y_data.append({'data':data, 'name':label})

    # Get the antisites concentration for each site
    for i in range(n):
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        for j in range(m):
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del-1:
                label = f'${sub_specie}_{{{site_specie}}}$'
            else:
                label = f'${sub_specie}_{{{site_specie}_{cur_ind}}}$'
            inds = specie_site_index_map[j]
            # Sum concentrations of multi-site species with same composition
            data = np.sum([conc[ind][i] for ind in range(*inds)], axis=0)
            data = data.tolist()
            y_data.append({'data':data, 'name':label})

    conc_data['y'] = y_data

    # Compute the vacancy formation energies
    def compute_vac_formation_energies(mu_vals):
        en = []
        for vac_def in vac_defs:
            site_specie = vac_def['site_specie']
            ind = specie_order.index(site_specie)
            uncor_energy = vac_def['energy']
            formation_energy = uncor_energy + mu_vals[ind]
            en.append(float(formation_energy))
        return en

    en_res = []
    for key in sorted(new_mu_dict.keys()):
        mu_val = new_mu_dict[key]
        en_res.append(compute_vac_formation_energies(mu_val))

    en_data = {'x_label':els[0] + ' mole fraction'}
    en_data['x'] = [dat[0][0] for dat in res1]         # x-axis data

    # Get vacancy formation energy for each site
    y_data = []
    for i, vac_def in enumerate(vac_defs):
        data = [data[i] for data in en_res]
        site_specie = vac_def['site_specie']
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = vac_string + site_specie + '}$'
        else:
            label = vac_string + site_specie + '_' + str(cur_ind) + '}$'
        y_data.append({'data':data, 'name':label})

    # Compute the antisite formation energies
    def compute_as_formation_energies(mu_vals):
        en = []
        for as_def in antisite_defs:
            site_specie = as_def['site_specie']
            sub_specie = as_def['substitution_specie']
            ind1 = specie_order.index(site_specie)
            ind2 = specie_order.index(sub_specie)
            uncor_energy = as_def['energy']
            formation_energy = uncor_energy + mu_vals[ind1] - mu_vals[ind2]
            en.append(formation_energy)
        return en

    en_res = []
    for key in sorted(new_mu_dict.keys()):
        mu_val = new_mu_dict[key]
        en_res.append(compute_as_formation_energies(mu_val))

    for i, as_def in enumerate(antisite_defs):
        data = [data[i] for data in en_res]
        site_specie = as_def['site_specie']
        sub_specie = as_def['substitution_specie']
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        if not specie_ind_del-1:
            label = '$' + sub_specie + '_{' + site_specie + '}$'
        else:
            label = '$' + sub_specie + '_{' + site_specie + '_' + \
                    str(cur_ind) + '}$'
        y_data.append({'data':data, 'name':label})

    en_data['y'] = y_data

    # Return chemical potential as well
    mu_data = {'x_label':els[0] + ' mole fraction', 'x':[]}
    mu_data['x'] = [dat[0][0] for dat in res1]         # x-axis data

    y_data = []
    for j in range(m):
        specie = specie_order[j]
        mus = [new_mu_dict[key][j] for key in sorted(new_mu_dict.keys())]
        y_data.append({'data':mus, 'name':specie})
    mu_data['y'] = y_data

    return conc_data, en_data, mu_data


def compute_defect_density(structure, e0, vac_defs, antisite_defs, T=800, \
                           trial_chem_pot = None, plot_style = "highcharts"):
    """
    Wrapper for the dilute_solution_model.
    The computed plot data is prepared based on plot_style.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        T: Temperature in Kelvin
        trial_chem_pot (optional): Trial chemical potentials to speedup
            the plot generation. Format is {el1:mu1, ...}
        plot_style (string): Allowed options are
            1) highcharts (default)
            2) gnuplot

    Returns:
        The plot data is generated and returned in asked format.
    """
    conc_data, en_data, mu_data = dilute_solution_model(
        structure, e0, vac_defs, antisite_defs, T, trial_chem_pot=trial_chem_pot)

    if plot_style == 'highcharts':
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data['xAxis'] = conc_data['x_label']
        hgh_chrt_data['yAxis'] = conc_data['y_label']

        series = []
        x = conc_data['x']
        for y_data in conc_data['y']:
            y = y_data['data']
            xy = zip(x,y)
            xy = [list(el) for el in xy]
            name = y_data['name'].strip('$')
            flds= name.split('_')
            def_string = flds[0]
            site_string = flds[1].strip('{}')
            name = def_string + "<sub>" + site_string + "</sub>"
            # series.append({'data':xy, 'name':y_data['name']})
            series.append({'data':xy, 'name':name})
        hgh_chrt_data['series'] = series
        return hgh_chrt_data
    elif plot_style == 'gnuplot':
        def data_to_rows(inp_data):
            rows = []
            labels = []
            labels.append(inp_data['x_label'])
            labels += [y['name'] for y in inp_data['y']]
            rows.append('#' + '\t'.join(labels))
            m = len(inp_data['x'])
            for i in range(m):
                data = []
                data.append(inp_data['x'][i])
                data += [y['data'][i] for y in inp_data['y']]
                data = [float(x) for x in data]
                rows.append('\t'.join(list(map(str, data))))
            return rows
        conc_rows = data_to_rows(conc_data)
        en_rows = data_to_rows(en_data)
        mu_rows = data_to_rows(mu_data)

        return conc_rows, en_rows, mu_rows
    else:
        raise ValueError(f'Plot style {plot_style} is not supported.')


# solute_site_preference_finder is based on dilute_solution_model and so most
# of the code is same. However, differences exist in setup and processing,
# hence the new function
def solute_site_preference_finder(structure, e0, vac_defs, antisite_defs, T, \
    solute_defs, solute_concen=0.01, trial_chem_pot = None):
    """
    Compute the solute defect densities using dilute solution model.
    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unitcell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from VASP calculation.
        T: Temperature in Kelvin
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite
            defect are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        solute_defs: List of solute defect parameters in the dictionary
            format. Similar to that of antisite defs, with solute specie
            specified in substitution_specie
        solute_concen: Solute concentration (in fractional value)
        trial_chem_pot: Trial chemical potentials to speed up the plot
            generation. Format is {el1:mu1,...}

    Returns:
        plot_data: The data for plotting the solute defect concentration.
    """

    if not check_input(vac_defs):
        raise ValueError('Vacancy energy is not defined')
    if not check_input(antisite_defs):
        raise ValueError('Antisite energy is not defined')

    # Setup the system
    site_species = [vac_def['site_specie'] for vac_def in vac_defs]
    solute_specie = solute_defs[0]['substitution_specie']
    site_species.append(solute_specie)
    multiplicity = [vac_def['site_multiplicity'] for vac_def in vac_defs]
    m = len(set(site_species))  # distinct species
    n = len(vac_defs)           # inequivalent sites
    # print(f'Site species: {site_species}')
    # print(f'Multiplicity: {multiplicity}')
    # print(f'Total number of distinct species: {m}')
    # print(f'Total number of inequivalent sites: {n}')

    # Reduce the system and associated parameters such that only distinctive
    # atoms are retained
    comm_div = gcd(*tuple(multiplicity))
    multiplicity = [val / comm_div for val in multiplicity]
    multiplicity.append(0)
    e0 = e0 / comm_div
    T = Float(T)
    # print(f'Energy: {e0}')

    c0 = np.diag(np.ones(n+1))
    c0[n, n] = 0
    mu = [Symbol(f'mu{i}') for i in range(m)]
    # print(f'Mu list: {mu}')

    # Generate maps for hashing
    # Generate specie->mu map and use it for site->mu map
    specie_order = []       # Contains hash for site->mu map    Eg: [Al, Ni, Fe]
    site_specie_set = set()             # e.g. {Ni, Al, Fe}
    for i in range(len(site_species)):
        site_specie = site_species[i]
        if site_specie not in site_specie_set:
            site_specie_set.add(site_specie)
            specie_order.append(site_specie)
    site_mu_map = []     # Eg: [mu0, mu1, mu2] where mu0->Al, mu1->Ni, mu2->Fe
    for i in range(len(site_species)):
        site_specie = site_species[i]
        j = specie_order.index(site_specie)
        site_mu_map.append(j)
    specie_site_index_map = []      # Eg: [(0,1), (1,2), (2,3)] for [Al, Ni ,Fe]
    for i in range(m):
        low_ind = site_species.index(specie_order[i])
        if i < m-1:
            high_ind = site_species.index(specie_order[i+1])
        else:
            high_ind = len(site_species)
        specie_site_index_map.append((low_ind, high_ind))
    # print(f'Specie order: {specie_order}')
    # print(f'Site specie set: {site_specie_set}')
    # print(f'Site mu map: {site_mu_map}')
    # print(f'Specie site index map: {specie_site_index_map}')

    """
    dC: delta concentration matrix:
    dC[i,j,k]: Concentration change of atom i, due to presence of atom
               j on lattice site k.
    Special case is [i,i,i] which is considered as vacancy.
    Few cases: (1) dC[i,i,i] = -1 due to being vacancy special case
               (2) dC[k,k,i] = +1 due to increment in k at i lattice if i
                                  lattice type is of different element
               (3) dC[i,k,i] = -1 due to decrement of ith type atom due to
                                  presence of kth type atom on ith sublattice
                                  and kth type atom specie is different from
                                  ith sublattice atom specie
               (4) dC[i,k,k] =  0 due to no effect on ith type atom
               (5) dC[i,j,k] =  0 if i != j != k
    """
    dC = np.zeros((n+1, n+1, n), dtype=np.int)
    for i in range(n):
        for k in range(n):
            for j in range(n+1):
                # case (2)
                if i == j and site_species[j] != site_species[k]:
                    dC[i, j, k] = 1
                # case (1) & (3)
                if i == k:
                    dC[i, j, k] = -1
                # case (4)
                if i != j and site_species[j] == site_species[k]:
                    dC[i, j, k] = 0
    # solute case
    dC[n, n, :] = 1
    # print(f'Delta concentration matrix:\n{dC}')

    # dE matrix: Flip energies (or raw defect energies)
    els = site_species.copy()
    dE = [[0 for _ in range(n)] for _ in range(n+1)]
    for j in range(n):
        site_specie = vac_defs[j]['site_specie']
        # Vacancy and antisite energies
        for i in range(n):
            if i == j:
                dE[i][j] = vac_defs[i]['energy']
            else:
                sub_specie = vac_defs[i]['site_specie']
                if site_specie == sub_specie:
                    dE[i][j] = 0
                else:
                    for as_def in antisite_defs:
                        if int(as_def['site_index']) == j+1 and \
                           sub_specie == as_def['substitution_specie']:
                            dE[i][j] = as_def['energy']
                            break
        # Solute energies
        for solute_def in solute_defs:
            def_site_ind = int(solute_def['site_index'])
            def_site_specie = solute_def['site_specie']
            if def_site_specie == site_specie and def_site_ind == j+1:
                dE[n][j] = solute_def['energy']
                break
    dE = np.array(dE)
    # print(f'Delta energies matrix:\n{dE}')

    # Initialization for concentrations
    # c(i, p) == presence of ith type atom on pth type site
    c = sympy.zeros(n+1, n)
    for i in range(n+1):
        for p in range(n):
            c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n+1):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) \
                              for j in range(n+1)])
                flip = Integer(dC[i, epi, p]) * \
                       exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    c[i, p] += flip
    total_c = []
    for ind in specie_site_index_map:
        val = 0
        for i in range(*ind):
            val += sum([c[i, j] * multiplicity[j] for j in range(n)])
        total_c.append(val)
    c_ratio = [total_c[i] / sum(total_c) for i in range(m)]  # different ratio from DSM!
    # print(f'Concentration matrix:\n{c}')
    # print(f'c_ratio:\n{c_ratio}')

    # Similar, just without the solute specie
    host_c = sympy.zeros(n, n)
    for i in range(n):
        for p in range(n):
            host_c[i, p] = Integer(c0[i, p])
            site_flip_contribs = []
            for epi in range(n):
                sum_mu = sum([mu[site_mu_map[j]] * Integer(dC[j, epi, p]) \
                              for j in range(n)])
                flip = Integer(dC[i, epi, p]) * \
                       exp(-(dE[epi, p] - sum_mu) / (k_B * T))
                if flip not in site_flip_contribs:
                    site_flip_contribs.append(flip)
                    host_c[i, p] += flip
    host_total_c = []
    for ind in specie_site_index_map[:-1]:
        val = 0
        for i in range(*ind):
            val += sum([host_c[i, j] * multiplicity[j] for j in range(n)])
        host_total_c.append(val)
    host_c_ratio = [host_total_c[i] / sum(host_total_c) for i in range(m-1)]

    # Expression for Omega, the Grand Potential
    omega1 = e0 - sum([mu[site_mu_map[i]] * c0[i, i] * multiplicity[i] \
                       for i in range(n)])
    omega = copy.deepcopy(omega1)
    used_dEs = []
    for p_r in range(n):
        for epi in range(n):
            sum_mu1 = sum([mu[site_mu_map[j]] * \
                          Float(dC[j, epi, p_r]) for j in range(n)])
            sum_mu = sum_mu1 - mu[site_mu_map[n]] * dC[n, epi, p_r]
            if p_r != epi and site_mu_map[p_r] == site_mu_map[epi]:
                continue
            if dE[epi, p_r] not in used_dEs:
                omega1 -= k_B * T * multiplicity[p_r] * \
                          exp(-(dE[epi, p_r] - sum_mu1) / (k_B * T))
                omega -=  k_B * T * multiplicity[p_r] * \
                          exp(-(dE[epi, p_r] - sum_mu) / (k_B * T))
                used_dEs.append(dE[epi, p_r])
    # print(f'Omega: {omega}')

    # Compute composition ranges
    max_host_specie_concen = 1 - solute_concen
    specie_concen = [sum(multiplicity[ind[0]:ind[1]]) \
                     for ind in specie_site_index_map]
    host_specie_concen_ratio = [specie_concen[i] / sum(specie_concen) * \
                                max_host_specie_concen for i in range(m)]
    host_specie_concen_ratio[-1] = solute_concen
    li = specie_site_index_map[0][0]
    hi = specie_site_index_map[0][1]
    comp1_min = sum(multiplicity[li:hi]) / sum(multiplicity) * \
                max_host_specie_concen - 0.01
    comp1_max = sum(multiplicity[li:hi]) / sum(multiplicity) * \
                max_host_specie_concen + 0.01
    delta = (comp1_max - comp1_min) / 50.0  # 0.04%
    # print(f'Composition range: [{comp1_min}, {comp1_max}]')
    # print(f'delta: {delta}')
    # print(f'max_host_specie_concen: {max_host_specie_concen}')
    # print(f'specie_concen: {specie_concen}')
    # print(f'host_specie_concen_ratio: {host_specie_concen_ratio}')

    def reduce_mu():
        host_concen = 1 - solute_concen
        new_c0 = c0.astype(float)
        for i in range(n):
            new_c0[i, i] = host_concen * c0[i, i]
        new_c0[n, n] = 2 * solute_concen  # double because of binary host
        omega = [e0 - sum([mu[site_mu_map[i]] * new_c0[i, i] for i in range(n+1)])]
        x = sympy.solve(omega)
        return x

    def compute_solute_mu_by_lin_search(host_mu_vals):
        mu_red = reduce_mu()
        # print(f'Reduced mu solution: {mu_red}')

        y_vect = host_specie_concen_ratio
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m-1)]
        vector_func.append(omega)
        mu_vals = None
        dm = 0.5
        m1_min = -15.0
        if e0 > 0:
            m1_max = 10            # Search space needs to be modified
        else:
            m1_max = 0
        for m1 in np.arange(m1_min, m1_max, dm):
            trial_mus = host_mu_vals + [m1]
            try:
                x = sympy.nsolve(vector_func, mu, trial_mus)
                if x:
                    mu_vals = [float(mu_val) for mu_val in x]
                break
            except:
                continue
        else:
            raise ValueError()
        return mu_vals

    def compute_mus():
        # Compute trial mu
        mu_red = reduce_mu()
        # print(f'Reduced mu solution: {mu_red}')

        y_vect = host_specie_concen_ratio
        vector_func = [y_vect[i] - c_ratio[i] for i in range(m)]
        vector_func.append(omega)
        mu_vals = None
        dm = 0.5
        m_min = -15.0
        if e0 > 0:
            m_max = 10            # Search space needs to be modified
        else:
            m_max = 0
        for m1 in np.arange(m_min, m_max, dm):
            for m2 in np.arange(m_min, m_max, dm):
                m0 = mu_red[mu[0]].subs([(mu[1], m1), (mu[2], m2)])
                try:
                    mu_vals = sympy.nsolve(vector_func, mu, [m0, m1, m2])
                    # Line needs to be modified to include all mus when n > 2
                except:
                    continue
                break
            if mu_vals:
                mu_vals = [float(mu_val) for mu_val in mu_vals]
                break
        else:
            raise ValueError("Couldn't find mus!")
        return mu_vals

    if not trial_chem_pot:
        print('No trial chemical potential is given.')
        # Try computing mus by assuming one of the defects is dominant at 0.01
        # concentration. First vacancy is tried and then antisite.

        # Generate trial mus assuming vacancy as dominant defect
        # for specie-0 at lower yval.
        li0 = specie_site_index_map[0][0]
        hi0 = specie_site_index_map[0][1]
        li1 = specie_site_index_map[1][0]
        hi1 = specie_site_index_map[1][1]
        spec_mult = [sum(multiplicity[li0:hi0]), sum(multiplicity[li1:hi1])]
        ln_def_conc = 4.60517
        # Try to initialize with host mus
        for i in range(li0, hi0):
            vac_flip_en =  vac_defs[i]['energy']
            mu_vals = [ln_def_conc * k_B * T - vac_flip_en]
            mu_vals.append((e0 - spec_mult[0] * mu_vals[0]) / spec_mult[1])
            comp_ratio = comp1_min

            # Test if the trial mus are good
            vector_func = [comp_ratio - host_c_ratio[0]]
            vector_func.append(omega1)
            try:
                # Get host mu values first, then solutes
                host_mu_vals = sympy.nsolve(vector_func, mu[:-1], mu_vals)
                if host_mu_vals:
                    host_mu_vals = [float(mu_val) for mu_val in host_mu_vals]
                compute_solute_mu_by_lin_search(host_mu_vals)
                break
            except:     # Go for antisite as dominant defect
                mu_gs = [Symbol(f'mu_gs{j}') for j in range(m-1)]
                eqs = [mu_gs[0] - mu_gs[1] - \
                       (ln_def_conc * k_B * T - antisite_defs[i]['energy'])]
                eqs.append(spec_mult[0] * mu_gs[0] + spec_mult[1] * mu_gs[1] - e0)
                x = sympy.solve(eqs, mu_gs)
                host_mu_vals = [x[key] for key in sorted(x.keys(), key=lambda inp: inp.name)]

                try:
                    host_mu_vals = sympy.nsolve(vector_func, mu[:-1], host_mu_vals)
                    if host_mu_vals:
                        host_mu_vals = [float(mu_val) for mu_val in host_mu_vals]
                    mu_vals = compute_solute_mu_by_lin_search(host_mu_vals)
                    break
                except:     # Go to the default option (search the space)
                    pass
        else:
            mu_vals = compute_mus()
    else:
        try:
            mu_vals = [trial_chem_pot[element] for element in specie_order]
        except:
            mu_vals = compute_mus()
    print(f'Trial mu_vals: {mu_vals}')

    # COMPUTE MU's FOR ALL COMPOSITION RATIOS in the range
    # +/- 1% from the stoichiometry
    result = {}
    for y in np.arange(comp1_min, comp1_max + delta, delta):
        y_vect = [y, max_host_specie_concen - y, solute_concen]
        vector_func = [y_vect[i] - c_ratio[i] for i in range(1, m)]
        vector_func.append(omega)

        try:
            x = sympy.nsolve(vector_func, mu, mu_vals)
            if x:
                mu_vals = [float(mu_val) for mu_val in x]
        except:
            continue
        result[y] = mu_vals

    res = []
    # COMPUTE THE CONCENTRATION for all the compositions
    for key in sorted(result.keys()):
        mu_val = result[key]
        total_c_val = [total_c[i].subs(dict(zip(mu, mu_val))) \
                       for i in range(len(total_c))]
        c_val = c.subs(dict(zip(mu, mu_val)))
        # Concentration of first element / total concen
        res1 = [float(total_c_val[0] / sum(total_c_val))]
        for i in range(n+1):
            for j in range(n):
                if i == j:              # Vacancy
                    vac_conc = float(exp(-(mu_val[site_mu_map[i]] + dE[i, i])/(k_B * T)))
                    res1.append(vac_conc)
                    if vac_conc > 1.0:
                        print('Warning! Vacancy concentration is > 1.')
                else:                   # Antisite
                    anti_conc = float(c_val[i, j] / c0[j, j])
                    res1.append(anti_conc)
                    if anti_conc > 1.0:
                        print('Warning! Antisite concentration is > 1.')
        res.append(res1)

    res = np.array(res)
    dtype = [(str('x'), np.float64)] + [(f'y{i}{j}', np.float64) \
                                        for i in range(n+1) for j in range(n)]
    res1 = np.sort(res.view(dtype), order=[str('x')], axis=0)

    conc_data = {}
    # Because all the plots have identical x-points, we store them in a single array
    conc_data['x'] = [dat[0][0] for dat in res1]         # x-axis data
    # Element whose composition is varied. For x-label
    conc_data['x_label'] = els[0] + " mole fraction"
    conc_data['y_label'] = "Point defect concentration"
    conc = [[[] for _ in range(n)] for _ in range(n+1)]
    for i in range(n+1): # Append vacancies
        for j in range(n):
            # n+1 x n array where each element is a list of concentrations
            conc[i][j] = [dat[0][i * n + j + 1] for dat in res1]


    y_data = []
    # Get the vacancies concentration for each site
    for i in range(n):
        data = conc[i][i]
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        vac_string = "$Vac_{"
        if not specie_ind_del-1:
            label = f'{vac_string}{site_specie}}}$'
        else:
            label = f'{vac_string}{site_specie}_{cur_ind}}}$'
        y_data.append({'data':data, 'name':label})

    # Get the antisites and solute concentration for each site
    for i in range(n):
        site_specie = els[i]
        specie_ind = site_mu_map[i]
        indices = specie_site_index_map[specie_ind]
        specie_ind_del = indices[1] - indices[0]
        cur_ind = i - indices[0] + 1
        for j in range(m):
            sub_specie = specie_order[j]
            if sub_specie == site_specie:
                continue
            if not specie_ind_del-1:
                label = f'${sub_specie}_{{{site_specie}}}$'
            else:
                label = f'${sub_specie}_{{{site_specie}_{cur_ind}}}$'
            inds = specie_site_index_map[j]
            # Sum concentrations of multi-site species with same composition
            data = np.sum([conc[ind][i] for ind in range(*inds)], axis=0)
            data = data.tolist()
            y_data.append({'data':data, 'name':label})

    conc_data['y'] = y_data
    return conc_data


def solute_defect_density(structure, e0, vac_defs, antisite_defs, solute_defs,
    solute_concen=0.01, T=800, trial_chem_pot=None, plot_style="highcharts"):
    """
    Wrapper for the solute_site_preference_finder.
    The computed plot data is prepared based on plot_style.

    Args:
        structure: pymatgen.core.structure.Structure object representing the
            primitive or unit cell of the crystal.
        e0: The total energy of the undefected system.
            This is E0 from the VASP calculation.
        vac_defs: List of vacancy defect parameters in the dictionary format.
            The keys of the dict associated with each vacancy defect are
            1) site_index, 2) site_specie, 3) site_multiplicity, and
            4) energy. 1-3 can be obtained from
            pymatgen.analysis.defects.point_defects.Vacancy class.
            Site index is expected to start with 1 (fortran index).
        antisite_defs: List of antisite defect parameters in the dictionary
            format. The keys of the dict associated with each antisite defect
            are 1) site_index, 2) site_specie, 3) site_multiplicity,
            4) substitution_specie, and 5) energy. 1-3 can be obtained
            from pymatgen.analysis.defects.point_defects.Vacancy class.
        solute_defs: List of solute defect parameters in the dictionary
            format. Similary to that of antisite defs, wtih solute specie
            specified in substitution_specie
        solute_concen: Solute concentration (in fractional value)
        T: Temperature in Kelvin
        trial_chem_pot (optional): Trial chemical potentials to speed up
            the plot generation. Format is {el1:mu1,...}
        plot_style (string): Allowed options are
            1) highcharts (default)
            2) gnuplot

    Returns:
        The plot data is generated and returned in asked format.
    """
    def_conc_data = solute_site_preference_finder(
        structure, e0, vac_defs, antisite_defs, T, solute_defs,
        solute_concen=solute_concen, trial_chem_pot=trial_chem_pot)

    if plot_style == 'highcharts':
        "Energy data is ignored in this mode"
        hgh_chrt_data = {}
        hgh_chrt_data['xAxis'] = def_conc_data['x_label']
        hgh_chrt_data['yAxis'] = def_conc_data['y_label']

        series = []
        x = def_conc_data['x']
        for y_data in def_conc_data['y']:
            y = y_data['data']
            xy = zip(x,y)
            xy = [list(el) for el in xy]
            name = y_data['name'].strip('$')
            flds= name.split('_')
            def_string = flds[0]
            site_string = flds[1].strip('{}')
            name = def_string + "<sub>" + site_string + "</sub>"
            #series.append({'data':xy, 'name':y_data['name']})
            series.append({'data':xy, 'name':name})
        hgh_chrt_data['series'] = series
        return hgh_chrt_data
    elif plot_style == 'gnuplot':
        def data_to_rows(inp_data, y_lbl_flg):
            rows = []
            labels = []
            labels.append(inp_data['x_label'])
            if y_lbl_flg:
                labels.append(inp_data['y_label'])
            else:
                labels += [y['name'] for y in inp_data['y']]
            rows.append('#' + '\t'.join(labels))
            m = len(inp_data['x'])
            for i in range(m):
                data = []
                data.append(inp_data['x'][i])
                data += [y['data'][i] for y in inp_data['y']]
                data = [float(x) for x in data]
                rows.append('\t'.join(list(map(str, data))))
            return rows
        pt_def_conc_rows = data_to_rows(def_conc_data, False)
        return pt_def_conc_rows
    else:
        raise ValueError(f'Plot style {plot_style} is not supported.')
