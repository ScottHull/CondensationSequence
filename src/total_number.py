import sys


def molecule_number_density(number_densities, molecule_stoich, K, temperature, R=8.314):
    """
    Calculates the number density of a given molecule.
    K = P_mol / prod(P_components)
    P_i = n_i RT -> K = P_mol / prod(n_i RT) -> P_mol = K prod(n_i RT)
    n_mol = P_mol / RT = (K prod(n_i RT)) / RT
    :param number_densities:
    :param element_stoichiometry:
    :param K:
    :param temperature:
    :param R:
    :return:
    """
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    total = K / (R * temperature)  # sum of all component element partial pressures
    for atom in molecule_stoich.keys():  # loop through elements present in molecule
        stoich = molecule_stoich[atom]  # get the stoich coefficient of the atom in the element
        if atom in fugacities:  # sum [ P_i * R * T ]^nu, nu = stoich coeff
            total *= (number_densities[atom] * (R * temperature)) ** (stoich / 2)
        else:
            total *= (number_densities[atom] * (R * temperature)) ** stoich
    return total


def calculate_total_N(gas_element_appearances_in_molecules, solid_element_appearances_in_molecules,
                      liquid_element_appearances_in_molecules,
                      element_number_densities,
                      condensing_solids, condensing_liquids, gas_molecule_library, solid_molecule_library,
                      liquid_molecule_library, K_dict, temperature,
                      R=8.314e-2):
    """
    Returns the total number density of all elements in the system based on their active gas and solid molecules.
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
    number_density_liquid = {}
    number_density_solid = {}

    # gasses
    for atom in gas_element_appearances_in_molecules.keys():  # loop through all of the input elements for gasses
        number_density_gas.update({atom: 0})
        molecules = gas_element_appearances_in_molecules[atom]  # retrieve the molecules in which the element exists
        for molecule in molecules:  # loop through all molecules in which the element appears
            K = K_dict[molecule]  # get the equilibrium constant
            molecule_stoich = gas_molecule_library[molecule]  # retrieve the molecule stoichiometry
            stoich = gas_molecule_library[molecule][atom]
            number_density = molecule_number_density(
                number_densities=element_number_densities,
                molecule_stoich=molecule_stoich,
                K=K,
                temperature=temperature,
                R=R
            )  # calculate number density
            num = stoich * number_density
            number_density_gas[atom] += num  # sum the number density

    # solids
    if len(condensing_solids) > 0:
        for atom in solid_element_appearances_in_molecules.keys():  # loop through all of the input elements for solids
            number_density_solid.update({atom: 0})
            molecules = solid_element_appearances_in_molecules[
                atom]  # retrieve the molecules in which the element exists
            for molecule in molecules:  # loop through all molecules in which the element appears
                if molecule in condensing_solids:
                    K = K_dict[molecule]  # get the equilibrium constant
                    molecule_stoich = solid_molecule_library[molecule]  # retrieve the molecule stoichiometry
                    for component_atom in molecule_stoich.keys():
                        stoich = molecule_stoich[component_atom]
                        num = stoich * element_number_densities[molecule]  # calculate number density
                        number_density_solid[atom] += num  # sum the number density

    # liquids
    if len(condensing_liquids) > 0:
        for atom in liquid_element_appearances_in_molecules.keys():  # loop through all of the input elements for liquids
            number_density_liquid.update({atom: 0})
            molecules = liquid_element_appearances_in_molecules[
                atom]  # retrieve the molecules in which the element exists
            for molecule in molecules:  # loop through all molecules in which the element appears
                if molecule in condensing_liquids:
                    K = K_dict[molecule]  # get the equilibrium constant
                    molecule_stoich = liquid_molecule_library[molecule]  # retrieve the molecule stoichiometry
                    for component_atom in molecule_stoich.keys():
                        stoich = molecule_stoich[component_atom]
                        num = stoich * element_number_densities[molecule]  # calculate number density
                        number_density_liquid[atom] += num  # sum the number density

    total_pressure = sum(number_density_gas.values()) + sum(number_density_solid.values()) + sum(
        number_density_liquid.values())
    return total_pressure
