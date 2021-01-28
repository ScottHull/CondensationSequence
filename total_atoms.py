def get_num(partial_pressures, element_stoichiometry, K, temperature, R=8.314):
    """
    Does this calculate number density???
    :param partial_pressures:
    :param element_stoichiometry:
    :param K:
    :param temperature:
    :param R:
    :return:
    """
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    total = 0
    for atom in element_stoichiometry:  # loop through elements present in molecule
        stoich = element_stoichiometry[atom]  # get the stoich coefficient of the atom in the element
        if atom in fugacities:  # sum [ P_i * R * T ]^nu, nu = stoich coeff
            total += (partial_pressures[atom] * (R * temperature)) ** (stoich / 2)
        else:
            total += (partial_pressures[atom] * (R * temperature)) ** stoich
    total += K / (R * temperature)
    return total


def calculate_total_atoms(gas_element_appearances_in_molecules, solid_element_appearances_in_molecules,
                          partial_pressures,
                          condensing_solids, gas_molecule_library, solid_molecule_library, K_dict, temperature,
                          R=8.314):
    """

    :param gas_element_appearances_in_molecules: a dictionary where the element is the key and the list of gas molecules it appears in is the value
    :param solid_element_appearances_in_molecules: a dictionary where the element is the key and the list of solid molecules it appears in is the value
    :param partial_pressures: gas activities calculated from the root finding function in main
    :param condensing_solids: a list of condensing solids
    :param K_dict: a dictionary of K values
    :param temperature: current temperature
    :param R: gas constant
    :return:
    """

    number_density_gas = {}
    number_density_solid = {}

    # gasses
    for atom in gas_element_appearances_in_molecules.keys():  # loop through all of the input elements for gasses
        number_density_gas.update({atom: 0})
        molecules = gas_element_appearances_in_molecules[atom]  # retrieve the molecules in which the element exists
        for molecule in molecules:  # loop through all molecules in which the element appears
            K = K_dict[molecule]  # get the equilibrium constant
            molecule_stoich = gas_molecule_library[molecule]  # retrieve the molecule stoichiometry
            for component_atom in molecule_stoich.keys():
                stoich = molecule_stoich[component_atom]
                num = stoich * get_num(partial_pressures=partial_pressures, element_stoichiometry=stoich, K=K,
                                       temperature=temperature,
                                       R=R)  # calcualte number density
                number_density_gas[atom] += num  # sum the number density

    # solids
    if len(condensing_solids) > 0:
        for atom in solid_element_appearances_in_molecules.keys():  # loop through all of the input elements for gasses
            number_density_solid.update({atom: 0})
            molecules = solid_element_appearances_in_molecules[
                atom]  # retrieve the molecules in which the element exists
            for molecule in molecules:  # loop through all molecules in which the element appears
                K = K_dict[molecule]  # get the equilibrium constant
                molecule_stoich = solid_molecule_library[molecule]  # retrieve the molecule stoichiometry
                for component_atom in molecule_stoich.keys():
                    stoich = molecule_stoich[component_atom]
                    num = stoich * partial_pressures[component_atom]  # calcualte number density
                    number_density_solid[atom] += num  # sum the number density

    total_pressure = sum(number_density_gas.values()) + sum(number_density_solid.values())
    return total_pressure
