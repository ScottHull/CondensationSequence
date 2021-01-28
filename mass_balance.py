import number_density
from math import log10
import sys


def mass_balance(guess_number_density, ordered_element_names, gas_element_appearances_in_molecules, gas_molecule_library, K_dict, partial_pressures_dict,
                 temperature, condensing_solids):

    number_density_dict_tmp = {}  # the partial_pressures list is adjusted by the algorithm per iteration, so need to rebuild a temperature dictionary during each internal loop
    for index, element in enumerate(ordered_element_names):
        number_density_dict_tmp.update({element: guess_number_density[index]})
    # gas_number_densities = number_density.calculate_molecule_number_density(
    #     element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
    #     K_dict=K_dict, partial_pressures=partial_pressures_dict, temperature=temperature)
    element_number_densities_dict = number_density.partial_pressure_element(element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
         K_dict=K_dict, guess_number_densities=number_density_dict_tmp, temperature=temperature)
    element_number_densities = [element_number_densities_dict[element] for element in ordered_element_names]

    # final_dict = {}
    # for element in number_density_dict_tmp.keys():
    #     # log10(P_el / N_el)
    #     final_dict.update({element: log10(partial_pressures_dict[element]) - log10(element_mass_balance_dict[element])})
    # finals = [final_dict[element] for element in ordered_element_names]
    # return finals

    return element_number_densities

