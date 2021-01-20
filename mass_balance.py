import number_density
from math import log10
import sys


def mass_balance(partial_pressures, ordered_element_names, gas_element_appearances_in_molecules, gas_molecule_library, K_dict, partial_pressures_dict,
                 temperature, condensing_solids):
    partial_pressures_dict_tmp = {}
    for index, element in enumerate(ordered_element_names):
        partial_pressures_dict_tmp.update({element: partial_pressures[index]})
    # gas_number_densities = number_density.calculate_molecule_number_density(
    #     element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
    #     K_dict=K_dict, partial_pressures=partial_pressures_dict, temperature=temperature)
    gas_number_densities = number_density.mass_balance_original(element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
         K_dict=K_dict, partial_pressures=partial_pressures_dict_tmp, temperature=temperature)
    gas_activities = [gas_number_densities[element] for element in ordered_element_names]

    final_dict = {}
    for element in partial_pressures_dict_tmp.keys():
        final_dict.update({element: log10(partial_pressures_dict[element]) - log10(gas_number_densities[element])})
    finals = [final_dict[element] for element in ordered_element_names]
    return finals

