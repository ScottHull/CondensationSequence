import numpy as np
from math import log10
from copy import copy
from collections import Counter
import sys
import re
import warnings


def number_density_element_gas(element_appearances_in_molecules, molecule_library, K_dict, guess_number_densities,
                               temperature, R=8.3144621e-2):
    """
    Calculate the mass balance equation for product molecules as a function of partial pressure of the reactants.
    The ideal gas law is PV = nRT -> n/V = P/RT -> n/V is number density n_i -> n_i = P_i / RT
    The number density "n" of a gas is taken from the ideal gas law.  For example, for H2O:
    The reaction is:
    H2 (g) + 0.5 O2 (g) <-> H2O (g)
    And the equilibrium constant K is:
    K_H2O = P_H2O / (P_H2 * (P_O2)^0.5)
    You can rewrite K in terms of the number density:
    n_H2O = P_H2O / RT
    n_H2 = P_H2 / RT
    n_O2 = P_O2 / RT
    K_H2O = P_H2O / (P_H2 * (P_O2)^0.5) = (P_H2O / RT) / [(P_H2 / RT) * (P_O2 / RT)^0.5]
    Therefore, you can rearrange this equation to solve for the number density n_H2O:
    n_H2O = n_H2 * (n_O2)^0.5 * K_H2O * (RT)^0.5
    Below I will write all number densities n_i in terms of the partial pressure: n_i = P_i / RT.
    This function returns elemental mass balance equations, which is the sum of number densities of molecules that contain the element.
    For example, the mass balance of O is:
    N_0 = n_H2O _ n_CO + 2n_CO2
    where the coefficent of 2 is because O appears twice.  Therefore, the general equation is number densities time the stoichiometry subscript of the element in the molecule.
    Note that diatomic molecules (i.e. O2, H2, Cl2, etc) are written monatomically (i.e. with coefficients stoich / 2).
    :param molecules:
    :param molecule_library:
    :param K_dict:
    :param partial_pressures:
    :param temperature:
    :return:
    """
    element_number_densities = {}
    diatomic_molecules = ['H', 'Cl', 'F', 'O', 'N']  # diatomic molecules
    for element in element_appearances_in_molecules.keys():  # element is the element (i.e. Mg, Si, Fe, etc)
        atom_number_density = 0  # e.g. N_O = n_H2O + n_CO + n_CO2 + ...
        molecules = element_appearances_in_molecules[element]  # get the list of molecules in which the element appears
        for molecule in molecules:  # e.g. H2O
            K = K_dict[molecule]  # K value of the molecule
            molecule_number_density = K / (R * temperature)  # the full number density is calculated down below, i.e. P_H2O = (K_H2O * prod(P_reactant_elements)) / RT
            molecule_stoich = molecule_library[molecule]  # get the stoichiometry of the molecule
            for component_element in molecule_stoich:  # e.g. {H: 2, O: 1}
                if guess_number_densities[component_element] > 0:  # negative guesses break the simulation
                    component_stoich = molecule_stoich[component_element]  # e.g. 2 for H
                    power = copy(component_stoich)
                    if component_element in diatomic_molecules:
                        power /= 2.0  # e.g. for H: 2/2 = 1
                    if component_element != element:  # only multiply by stoich coefficient if this is the element in focus
                        component_stoich = 1
                    # i.e. n_H2O = K_H2O * (RT)^-1 * (n_H2 RT) * (n_O2 RT)^0.5
                    # for example, this function would multiply on the (n_O2 RT)^0.5 to the equation above
                    # component stioch adds on a coefficient, so it would add 2n_CO2 calculating N_O like
                    molecule_number_density *= component_stoich * (
                            guess_number_densities[component_element] * R * temperature) ** power
                else:  # force the solver to re-guess
                    molecule_number_density = 1e999
                    break
            atom_number_density += molecule_number_density  # i.e. N_O += 2n_CO2
        element_number_densities.update({element: atom_number_density})
    return element_number_densities


def get_all_condensing_elements_and_molecules(condensing_molecules):
    """
    Returns a dictionary of all condensing elements and a list of condensing solids in which the element appears.
    :param condensing_solids:
    :return:
    """
    elements_and_molecules = {}
    for molecule in condensing_molecules:
        stoich = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)
        for element in stoich:
            element = element[0]  # get the element name from the element-stoich tuple returned by the re function
            if element not in elements_and_molecules.keys():
                elements_and_molecules.update({element: []})
            elements_and_molecules[element].append(molecule)
    return elements_and_molecules


def number_density_element_solid(number_densities, condensing_solids, solid_molecule_library, ordered_names):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    solid_number_density = {}
    solids = get_all_condensing_elements_and_molecules(condensing_molecules=condensing_solids)
    for element in solids.keys():
        mass_balance = 0
        molecules = solids[element]
        for solid_molecule in molecules:
            stoich = solid_molecule_library[solid_molecule]
            if number_densities[solid_molecule] > 0:  # negative guesses break the solver
                coeff = stoich[element]
                if element in fugacities:
                    coeff /= 2.0
                mass_balance += coeff * number_densities[solid_molecule]
            else:  # force the solver out to re-sample
                mass_balance = 1e999
        solid_number_density.update({element: mass_balance})
    return solid_number_density


def number_density_element_liquid(number_densities, condensing_liquids, liquid_molecule_library, ordered_names):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    liquid_number_density = {}
    liquids = get_all_condensing_elements_and_molecules(condensing_molecules=condensing_liquids)
    for element in liquids.keys():
        mass_balance = 0
        molecules = liquids[element]
        for liquid_molecule in molecules:
            stoich = liquid_molecule_library[liquid_molecule]
            if number_densities[liquid_molecule] > 0:  # negative guesses break the solver
                coeff = stoich[element]
                if element in fugacities:
                    coeff /= 2.0
                mass_balance += coeff * number_densities[liquid_molecule]
            else:  # force the solver out to re-sample
                mass_balance = 1e999
        liquid_number_density.update({element: mass_balance})
    return liquid_number_density


def solid_partial_pressure(molecule, guess_number_densities, solid_molecules_library, temperature, R=8.3144621e-2):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    molecule_partial_pressure = 1
    stoich = solid_molecules_library[molecule]
    for element in stoich.keys():
        if guess_number_densities[element] < 0:  # kick solver if it is guessing negative values
            molecule_partial_pressure = 1e999
            break
        else:
            coeff = stoich[element]
            if element in fugacities:
                coeff /= 2.0
            molecule_partial_pressure *= (guess_number_densities[element] * R * temperature) ** coeff
    return molecule_partial_pressure


def liquid_partial_pressure(molecule, guess_number_densities, liquid_molecules_library, temperature, R=8.3144621e-2):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    molecule_partial_pressure = 1
    stoich = liquid_molecules_library[molecule]
    for element in stoich.keys():
        if guess_number_densities[element] < 0:  # kick solver if it is guessing negative values
            molecule_partial_pressure = 1e999
            break
        else:
            coeff = stoich[element]
            if element in fugacities:
                coeff /= 2.0
            molecule_partial_pressure *= (guess_number_densities[element] * R * temperature) ** coeff
    return molecule_partial_pressure
