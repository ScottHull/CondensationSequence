def get_num(gas_activities, element_stoichiometry, K, temperature, R=8.314):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    total = 0
    for atom in element_stoichiometry:  # loop through elements present in molecule
        stoich = element_stoichiometry[atom]  # get the stoich coefficient of the atom in the element
        if atom in fugacities:
            total += (gas_activities[atom] * (R * temperature)) ** (stoich / 2)
        else:
            total += (gas_activities[atom] * (R * temperature)) ** stoich
    total += K - (R  * temperature)



def calculate_total_atoms(gas_element_appearances_in_molecules, solid_element_appearances_in_molecules, gas_activities,
                          condensing_solids, gas_molecule_library, solid_molecule_library, K_dict, temperature,
                          R=8.314):
    """

    :param gas_element_appearances_in_molecules: a dictionary where the element is the key and the list of gas molecules it appears in is the value
    :param solid_element_appearances_in_molecules: a dictionary where the element is the key and the list of solid molecules it appears in is the value
    :param gas_activities: gas activities calculated from the root finding function in main
    :param condensing_solids: a list of condensing solids
    :param K_dict: a dictionary of K values
    :param temperature: current temperature
    :param R: gas constant
    :return:
    """

    for atom in gas_element_appearances_in_molecules.keys():  # loop through all of the input element
        molecules = gas_element_appearances_in_molecules[atom]  # retrieve the molecules in which the element exists
        s = 0
        for molecule in molecules:  # loop through all molecules in which the element appears
            stoich = gas_molecule_library[molecule]  # retrieve the molecule stoichiometry
            K = K_dict[molecule]  # get the equilibrium constant
            num = get_num(gas_activities=gas_activities, element_stoichiometry=stoich, K=K, temperature=temperature,
                          R=R)
