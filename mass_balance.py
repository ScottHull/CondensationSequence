import os
import sys
import re
from collections import Counter
from math import exp, log, log10
import collect_data


def get_equations(condensing_solids):
    Names = condensing_solids
    X_dict = {}
    Total = []
    Element_dict = {}
    Molecule_dict = {}

    for i in Names:
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        X_dict.update(dict(Mol_name))
        No_numbers = []
        for j in Mol_name:
            No_numbers.append(j[0])
            Total.append(j[0])
        Molecule_dict.update({i: Mol_name})
    X = []

    for i, j in dict.iteritems(X_dict):
        X.append(i)
        Element_dict.update({i: []})

    for i in Names:
        No_numbers = []
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)
        for j in Mol_name:
            No_numbers.append(j[0])
        for k in No_numbers:
            Z = Element_dict.get(k)
            Z.append(i)
            Element_dict[k] = Z
    Z = Counter(Total)
    for i, j in dict.iteritems(Molecule_dict):
        entry = []
        for k in j:
            entry.append((k[0], float(k[1])))
        Molecule_dict[i] = entry

    return Molecule_dict, X, Element_dict, Z


def mass_balance(element_appearances, K_dict, partial_pressures, gasses, temperature, condensing_solids):
    # note that element appearances should be a dict where the cation/anion is the key and a list of molecules that it appears in is the value
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    R = 8.3144621e-2  # units of L bar mol^-1 K^-1

    RT = R * temperature
    guess_dict = {}
    for key in partial_pressures.keys():
        guess_dict.update({key: 0})

    final_dict = {}
    finals = []
    out_dict = {}
    # print guess_dict
    for element in element_appearances.keys():
        out = 0
        for molecule in element_appearances[element]:
            # print k
            coefs = gasses[molecule]  # open an individual gas species
            entry = 0.0  # what is this?
            for atom in coefs:
                stoich_coeff = coefs[atom]
                if guess_dict[atom] <= 0:  # if the partial pressure becomes negative
                    entry = 1.e999
                    break
                else:
                    if atom == element:
                        if atom in fugacities:
                            entry = entry + ((stoich_coeff / 2.0) * log10(partial_pressures[atom] * RT) + log10(
                                stoich_coeff))
                        else:
                            entry = entry + (stoich_coeff * log10(partial_pressures[atom] * RT) + log10(stoich_coeff))
                    else:
                        if atom in fugacities:
                            entry = entry + ((stoich_coeff / 2.0) * log10(partial_pressures[atom] * RT))
                        else:
                            entry = entry + (stoich_coeff * log10(partial_pressures[atom] * RT))
            entry = 10 ** (entry + (K_dict[molecule]) - log10(RT))
            out += entry
        out_dict.update({element: out})

    for element in partial_pressures.keys():
        pp = partial_pressures[element]  # get the partial pressure value of the element
        final_dict.update({element: log10(pp) - log10(out_dict[element])})

    if len(condensing_solids) > 0:
        solid_molecule_dict, Y, solid_element_dict, Z = get_equations(condensing_solids)
        for molecule in solid_element_dict:  # open an element
            Sum = 0.
            for k in j:  # open molecule that contains that element
                if guess_dict[k] <= 0.:
                    Sum = 1.e999
                    break
                else:
                    for x in solid_molecule_dict.get(k):
                        if x[0] == i:
                            if x[0] in fugacities:
                                Sum += (float(x[1]) / 2.) * guess_dict[k]
                            else:
                                Sum += float(x[1]) * guess_dict[k]

            final_dict.update({i: 1. - ((log10(out_dict[i] + Sum) / log10(partial_pressures[i])))})
