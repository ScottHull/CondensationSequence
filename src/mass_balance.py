from src import number_density
from math import log10
import sys


def mass_balance(guess_number_density, ordered_names, gas_element_appearances_in_molecules,
                 gas_molecule_library, K_dict, number_densities_dict,
                 temperature, condensing_solids, solid_molecules_library, condensing_liquids, liquid_molecules_library):
    """
    This function works in conjunction with a root finder to find roots.
    This function takes guess number densities as an input and attempts to force the outputs to 0
    in order to find the roots, i.e. the "real" number densities.  The root finder then returns a
    list of the "real" number densities by exercising this function.

    The gas number density solver works through application of the ideal gas law, and the 0 criterion enforces that the
    mass balance of the gas element calculated through the partial pressure proportionality
    (equation 5, Unterborn & Panero 2017) and the ideal gas law calculated done here are equivalent.

    The solid number density solver works by assuming that the mass balance of the solid molecule is a direct
    proportionality between the guessed solid molecule number density and its elemental stoichiometry.  We then adjust
    the guess by forcing the partial pressure of the molecule to be 1.
    """

    mass_balance_zero = {}  # the y-values of the function with number densities as the x-input, want to be all 0's so that the x-values are roots.
    number_density_dict_tmp = {}  # the guess_number_density list is adjusted by the algorithm per iteration, so need to rebuild a temperature dictionary during each internal loop
    for index, element in enumerate(ordered_names):
        number_density_dict_tmp.update({element: guess_number_density[index]})

    # calculate gas element number density through guesses
    gas_element_mass_balance = number_density.number_density_element_gas(
        element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
        K_dict=K_dict, guess_number_densities=number_density_dict_tmp,
        temperature=temperature)  # application of the ideal gas law
    liquid_element_mass_balance = number_density.number_density_element_liquid(number_densities=number_density_dict_tmp,
                                                                               condensing_liquids=condensing_liquids,
                                                                               liquid_molecule_library=liquid_molecules_library,
                                                                               ordered_names=ordered_names)
    solid_element_mass_balance = number_density.number_density_element_solid(number_densities=number_density_dict_tmp,
                                                                             condensing_solids=condensing_solids,
                                                                             solid_molecule_library=solid_molecules_library,
                                                                             ordered_names=ordered_names)  # mass balance by a direct relationship between the guessed number density and the elemental stoichiometry

    for element in gas_element_mass_balance.keys():
        gas_mb = gas_element_mass_balance[element]
        zero = log10(number_densities_dict[element]) - log10(
            gas_mb)  # want the difference between the partial pressure proportionality and calculated number density to be 0, or the ratio to be 1, log10 helps speed calculation
        mass_balance_zero.update({element: zero})

    if len(condensing_solids) > 0:  # only perform if there are any condensing solids
        for element in solid_element_mass_balance:
            solid_mb = solid_element_mass_balance[element]
            zero = 1.0 - (log10(gas_element_mass_balance[element] + solid_mb) / log10(number_densities_dict[element]))
            mass_balance_zero.update({element: zero})

    # TODO: the sum in this zero might be invalid if there is also a solid
    if len(condensing_liquids) > 0:  # only perform if there are any condensing liquids
        for element in liquid_element_mass_balance:
            liquid_mb = liquid_element_mass_balance[element]
            zero = 1.0 - (log10(gas_element_mass_balance[element] + liquid_mb) / log10(number_densities_dict[element]))
            mass_balance_zero.update({element: zero})

    # if there are condensed elements, we want to force its partial pressure (i.e. activity) to 1.
    for solid in condensing_solids:
        partial_pressure = number_density.solid_partial_pressure(molecule=solid,
                                                                 guess_number_densities=number_density_dict_tmp,
                                                                 solid_molecules_library=solid_molecules_library,
                                                                 temperature=temperature)  # partial pressure product of component elements
        solid_criterion_zero = 1.0 - (partial_pressure * K_dict[
            solid])  # P_mol = K * prod(P_componets), want P_mol = 1 (i.e. activity = 1) for condensation criterion
        mass_balance_zero.update({solid: solid_criterion_zero})

    # if there are condensed elements, we want to force its partial pressure (i.e. activity) to 1.
    for liquid in condensing_liquids:
        partial_pressure = number_density.liquid_partial_pressure(molecule=liquid,
                                                                  guess_number_densities=number_density_dict_tmp,
                                                                  liquid_molecules_library=liquid_molecules_library,
                                                                  temperature=temperature)  # partial pressure product of component elements
        liquid_criterion_zero = 1.0 - (partial_pressure * K_dict[
            liquid])  # P_mol = K * prod(P_componets), want P_mol = 1 (i.e. activity = 1) for condensation criterion
        mass_balance_zero.update({liquid: liquid_criterion_zero})

    ordered_mass_balance_zero = [mass_balance_zero[name] for name in ordered_names]

    # the aim is to return 0's, which the solver will work to find
    return ordered_mass_balance_zero
