from src import number_density


# def __get_O_mass_balance(number_densities, ordered_names, gas_element_appearances_in_molecules,
#                          gas_molecule_library, K_dict,
#                          temperature, condensing_solids, solid_molecules_library, condensing_liquids,
#                          liquid_molecules_library):
#     mass_balance_zero = {}  # the y-values of the function with number densities as the x-input, want to be all 0's so that the x-values are roots.
#     number_density_dict_tmp = {}  # the number_densities list is adjusted by the algorithm per iteration, so need to rebuild a temperature dictionary during each internal loop
#     for index, element in enumerate(ordered_names):
#         number_density_dict_tmp.update({element: number_densities[index]})
#
#     # calculate gas element number density through guesses
#     gas_element_mass_balance = number_density.number_density_element_gas(
#         element_appearances_in_molecules=gas_element_appearances_in_molecules, molecule_library=gas_molecule_library,
#         K_dict=K_dict, guess_number_densities=number_density_dict_tmp,
#         temperature=temperature)  # application of the ideal gas law
#     liquid_element_mass_balance = number_density.number_density_element_liquid(number_densities=number_density_dict_tmp,
#                                                                                condensing_liquids=condensing_liquids,
#                                                                                liquid_molecule_library=liquid_molecules_library,
#                                                                                ordered_names=ordered_names)
#     solid_element_mass_balance = number_density.number_density_element_solid(number_densities=number_density_dict_tmp,
#                                                                              condensing_solids=condensing_solids,
#                                                                              solid_molecule_library=solid_molecules_library,
#                                                                              ordered_names=ordered_names)  # mass balance by a direct relationship between the guessed number density and the elemental stoichiometry
#     print(solid_element_mass_balance)
#
#
# def calculate_refractory_oxygen_percent(number_densities, ordered_names, gas_element_appearances_in_molecules,
#                                         gas_molecule_library, K_dict,
#                                         temperature, condensing_solids, solid_molecules_library, condensing_liquids,
#                                         liquid_molecules_library):
#     mass_balance_oxygen = __get_O_mass_balance(
#         number_densities=number_densities,
#         ordered_names=ordered_names,
#         gas_element_appearances_in_molecules=gas_element_appearances_in_molecules,
#         gas_molecule_library=gas_molecule_library,
#         K_dict=K_dict,
#         temperature=temperature,
#         condensing_solids=condensing_solids,
#         solid_molecules_library=solid_molecules_library,
#         condensing_liquids=condensing_liquids,
#         liquid_molecules_library=liquid_molecules_library
#     )
#     oxygen_number_solid = 0
#     refractory_oxygen_in_phase = {}
#     for solid in number_densities.keys():
#         if "_s" in solid:  # just making sure that it is a solid:
#             stoich = solid_molecules_library[solid]
#             if "O" in stoich.keys():
#                 O_stoich = stoich["O"]
#                 O_in_phase = O_stoich * number_densities[solid]
#                 oxygen_number_solid += O_in_phase
#                 refractory_oxygen_in_phase.update({solid: O_in_phase / mass_balance_oxygen * 100.0})
#     return oxygen_number_solid / mass_balance_oxygen * 100.0, refractory_oxygen_in_phase

def calculate_refractory_oxygen_percent(number_densities, solid_molecule_library, mass_balance):
    mass_balance_O = mass_balance["O"]
    refractory_O_pct = {}
    total_refractory_O_pct = 0
    for solid in number_densities.keys():
        if "_s" in solid:
            stoich = solid_molecule_library[solid]
            if "O" in stoich.keys():
                O_stoich = stoich['O']
                O_moles = O_stoich * number_densities[solid]
                refractory_O_pct.update({solid: O_moles / mass_balance_O * 100.0})
                total_refractory_O_pct += O_moles
    print(total_refractory_O_pct / mass_balance_O * 100.0, refractory_O_pct)
    return total_refractory_O_pct / mass_balance_O * 100.0, refractory_O_pct


