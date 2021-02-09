import number_density
from math import log10
import sys


def mass_balance(guess_number_density, ordered_names, gas_element_appearances_in_molecules,
                 gas_molecule_library, K_dict, number_densities_dict,
                 temperature, condensing_solids, solid_molecules_library):
    """
    This function works in conjunction with a root finder to find roots.
    This function takes guess number densities as an input and attempts to force the outputs to 0
    in order to find the roots, i.e. the "real" number densities.  The root finder then returns a
    list of the "real" number densities by exercising this function.
    :param guess_number_density: The initial guess for the algorithm to proceed.
    :param ordered_names: A list of gas *elements* and condensed solid *molecules*.
    :param gas_element_appearances_in_molecules:
    :param gas_molecule_library:
    :param K_dict:
    :param number_densities_dict:  We force n_i to produce N_el to match equation 5 of Unterborn & Panero 2017.
    :param temperature:
    :param condensing_solids:
    :return:
    """

    mass_balance_zero = {}  # force zero to be 0
    number_density_dict_tmp = {}  # the guess_number_density list is adjusted by the algorithm per iteration, so need to rebuild a temperature dictionary during each internal loop
    for index, element in enumerate(ordered_names):
        number_density_dict_tmp.update({element: guess_number_density[index]})

    # calculate gas element number density through guesses
    gas_element_mass_balance = number_density.number_density_element_gas(
        element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
        K_dict=K_dict, guess_number_densities=number_density_dict_tmp, temperature=temperature)
    solid_element_mass_balance = number_density.number_density_element_solid(guess_number_densities=guess_number_density, condensing_solids=condensing_solids, solid_molecule_library=solid_molecules_library, ordered_names=ordered_names)

    for element in gas_element_mass_balance.keys():
        gas_mb = gas_element_mass_balance[element]
        zero = log10(number_densities_dict[element]) - log10(gas_mb)  # want the difference between the partial pressure proportionality and calculated number density to be 0, or the ratio to be 1, log10 helps speed calculation
        mass_balance_zero.update({element: zero})

            
    for element in solid_element_mass_balance:
        solid_mb = solid_element_mass_balance[element]
        zero = 1 - ((gas_element_mass_balance[element] + solid_mb) / number_densities_dict[element])
        mass_balance_zero.update({element: zero})

    # if there are condensed elements, we want to force its partial pressure (i.e. activity) to 1.
    for solid in condensing_solids:
        partial_pressure = number_density.solid_partial_pressure(molecule=solid,
                                                                 guess_number_densities=number_density_dict_tmp,
                                                                 solid_molecules_library=solid_molecules_library,
                                                                 temperature=temperature)
        solid_criterion_zero = 1 - (partial_pressure / -K_dict[solid])
        mass_balance_zero.update({solid: solid_criterion_zero})
        
    ordered_mass_balance_zero = [mass_balance_zero[name] for name in ordered_names]

    # the aim is to return 0's, which the solver will work to find
    return ordered_mass_balance_zero

