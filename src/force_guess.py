import numpy as np


def get_guess(molecule, mass_balance, total_n, percent_condensed, elements_in_solid, temperature, number_densities):
    small_percentage_guess = False
    for element in elements_in_solid:
        if percent_condensed[element][-1]['percent'] < 95 and element != "O":
            small_percentage_guess = True
            break

    if small_percentage_guess:
        if molecule == 'Fe1_s' or molecule == 'Ni1_s' or molecule == 'Co1_s':
            if molecule == 'Fe1_s':
                return 8.e-14
            elif molecule == 'Ni1_s':
                return 8.e-14
            elif molecule == 'Co1_s':
                return 8.e-14
        else:
            if molecule != 'Mg1Si1O3_ortho_s':
                if molecule == 'Ca1Al4O7_s':
                    return 3.e-14
                elif molecule == 'Si1O2_s':
                    return 1.e-6 * total_n
                elif molecule == 'Ca2Mg1Si2O7_s':
                    return 1.e-13
                elif molecule == 'Ca2Al2Si1O7_s':
                    return 8.e-14
                else:
                    return 5.e-14
            else:
                return 2.e-14
    else:
        if temperature < 1300. and 'Ti' in elements_in_solid:
            if 'Fe' in elements_in_solid:
                return 1.e-5 * number_densities.get('Fe')
            else:
                return 1.e-2 * number_densities.get('Ti')
        elif temperature > 1300. and 'Ti' in elements_in_solid:
            return 3.e-1 * number_densities.get('Ti')
        else:
            return 9.e-8 * total_n


def kick_old_guesses(names, condensing_solids, condensing_liquids, number_densities):
    new_number_density_guess = list(number_densities)
    for solid in condensing_solids:
        if solid != condensing_solids[-1] and condensing_solids[-1] != 'Fe1_s' and condensing_solids[-1] != 'Co1_s' and \
                condensing_solids[-1] != 'Ni1_s':
            index = names.index(solid)
            new_number_density_guess[index] *= 0.88
    for liquid in condensing_liquids:
        if liquid != condensing_liquids[-1]:
            index = names.index(liquid)
            new_number_density_guess[index] *= 0.88
    return np.array(new_number_density_guess)
