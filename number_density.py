import numpy as np
from math import log10
import sys


def calculate_molecule_number_density(element_appearances_in_molecules, molecule_library, K_dict, partial_pressures,
                                      temperature, R=8.314):
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
    Below I will write all number densities n_i in terms of the partial pressure: n_i = P_i / RT
    :param molecules:
    :param molecule_library:
    :param K_dict:
    :param partial_pressures:
    :param temperature:
    :return:
    """
    mass_balance = {}
    diatomic_molecules = ['H', 'Cl', 'F', 'O', 'N']  # will have coefficients of 2 in balanced reactions
    for atom in element_appearances_in_molecules.keys():  # atom is the element (i.e. Mg, Si, Fe, etc)
        molecule_number_densities = []
        for molecule in element_appearances_in_molecules[atom]:  # get the list of molecules in which the element appears and loop through it
            if atom in diatomic_molecules:
                pass
            else:
                molecule_stoich = molecule_library[molecule]  # get the stoichiometry of the molecule
                K = K_dict[molecule]  # get the molecule's equilibrium constant K
                reactant_partial_pressures_sum = np.prod([partial_pressures[i] ** (1.0 / molecule_stoich[i]) for i in
                                                          molecule_stoich.keys()])  # this is the denominator of the K equation
                partial_pressure_molecule = K * reactant_partial_pressures_sum  # prod(P_products) = K * prod(P_reactants)
                number_density_molecule = partial_pressure_molecule / (R * temperature)  # n_i = P_i / RT
                molecule_number_densities.append(molecule_stoich[atom] * number_density_molecule)  # multiply number density by stoich coefficient
        mass_balance_atom = sum(molecule_number_densities)  # the total mass balance of the element is the sum of number densities with stoich coefficients
        mass_balance.update({atom: mass_balance_atom})
        if atom == "Cr":
            print(mass_balance_atom)
    return mass_balance


def mass_balance_original(element_appearances_in_molecules, molecule_library, K_dict, partial_pressures,
                                      temperature, R=8.3144621e-2):
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    RT = R * temperature
    out_dict = {}
    # print guess_dict
    for i in element_appearances_in_molecules.keys():  # for element, gas appearance molecule in element dict
        j = element_appearances_in_molecules[i]
        out = 0.
        for k in j:  # for gas molecule in the list
            # print k
            coefs = molecule_library.get(k)  # get the stoich of the gas molecule
            entry = 0.
            for x in coefs.keys():  # for element, stoich in the dict
                y = coefs[x]
                if partial_pressures[x] <= 0.:
                    entry = 1.e999
                    break
                else:
                    if x == i:  # if the stoich element equals the input element
                        if x in fugacities:
                            # sum [  (nu_i / 2) * log(P_i * R * T + log(nu_i) ] -> raise to power of 10 -> sum [  nu_i (P_i * R * T)^(n_i / 2)  ]
                            entry = entry + ((y / 2.) * log10(partial_pressures.get(x) * RT) + log10(y))
                        else:  # sum  [ nu_i * log(P_I * R * T) + log(y) ] -> raise to power of 10 -> sum [  nu_i (P_i * R * T)^(n_i)  ]
                            entry = entry + (y * log10(partial_pressures.get(x) * RT) + log10(y))
                    else:  # if the stoich element does NOT equal the input element
                        if x in fugacities:  # sum [  (n_i / 2) * log(P_i * R * T)  ]  -> raise to power of 10 ->  sum [  (P_i * R * T)^(n_i / 2)  ]
                            entry = entry + ((y / 2.) * log10(partial_pressures.get(x) * RT))
                        else:  # sum [  (n_i) * log(P_i * R * T)  ]  -> raise to power of 10 ->  sum [  (P_i * R * T)^(n_i)  ]
                            entry = entry + (y * log10(partial_pressures.get(x) * RT))
            # sum((P_i R T)^nu_i) + K/RT
            entry = entry + (log10(K_dict.get(k))) - log10(RT)
            entry = pow(10., entry)
            out += entry

        out_dict.update({i: out})
    return out_dict
