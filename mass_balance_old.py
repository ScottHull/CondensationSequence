import os
import sys
import re
from collections import Counter
from math import exp, log, log10
import collect_data


def get_equations(condensing_solids):
    Names = condensing_solids  # molecule names of all condensing solid atoms
    X_dict = {}
    Total = []
    Element_dict = {}
    Molecule_dict = {}

    for i in Names:  # for atom in atoms
        Mol_name = re.findall(r'([A-Z][a-z]*)(\d*)', i)  # stoich of condensing solid molecules from atom names
        X_dict.update(dict(Mol_name))  # X_dict is a dictionary that contains the stoich
        No_numbers = []
        for j in Mol_name:
            No_numbers.append(j[0])  # append the atom name
            Total.append(j[0])  # append hte a
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



def Mass_balance_fun(guess_list, *args):
    # Element_dict is a dictionary of elements from given abundances and all of the gas molecules in which it appears
    # K is the dictionary of solid and gas K values
    # Par_Pres is a dictionary of partial pressures of given abundance elements
    # gasses is a dictionary of gasses and their corresponding stoichiometry
    # condensing_solids is a list of condensing solids
    # guess is the same as Par_Pres but in list form and corresponds to the list of inputted elements (note that these values may reference and older Par_pres dict?)
    fugacities = ['H', 'Cl', 'F', 'O', 'N']
    R = 8.3144621e-2
    # units of L bar mol^-1 K^-1
    Element_dict = args[0]
    K = args[1]
    Par_Pres = args[2]
    gasses = args[3]
    Name = args[4]
    RT = R * args[5]
    condensing_solids = args[6]
    guess_dict = {}

    test = {}
    counter = 0
    for i in Name:  # loop through given input elements list
        guess_dict[i] = guess_list[counter]  # build a guess dict with the element as the key and the (old) partial pressure as the value
        counter += 1

    final = []  # Sum of all activity terms for individual element
    outs = 1
    out = []
    final_dict = {}
    finals = []
    coefs = {}
    out_dict = {}

    for i, j in dict.iteritems(Element_dict):  # i is element, j is list of molecules in which it exists
        out = 0.
        for k in j:  # for molecule in the list of molecules that contains k
            # print k
            coefs = gasses.get(k)  # get the number of atoms in the gas species
            entry = 0.
            for x, y in dict.iteritems(coefs):  # x is the atom in the molecule, y is the atom's stoichiometry
                if guess_dict[x] <= 0.:  # if the partial pressure of the atom is < 0
                    entry = 1.e999
                    break
                else:
                    # these look like the mass balance equations for gasses given by Mg/Si paper
                    if x == i:  # if the atom from the molecule equals in the input atom
                        if x in fugacities:  # if the atom is in the fugacities list (2 from diatomic molecules?)
                            entry += ((y / 2.) * log10(guess_dict.get(x) * RT) + log10(y))
                        else:
                            # this is a summation: sum(nu * log(P_i * RT) * log(nu)) = sum(nu * log(nu * P_i + nu * RT)
                            entry += + (y * log10(guess_dict.get(x) * RT) + log10(y))
                    else:
                        if x in fugacities:
                            entry += + ((y / 2.) * log10(guess_dict.get(x) * RT))
                        else:
                            entry += entry + (y * log10(guess_dict.get(x) * RT))
            entry += entry + (K.get(k)) - log10(RT)
            entry = pow(10., entry)
            out += entry

        out_dict.update({i: out})
    for i, j in dict.iteritems(Par_Pres):
        final_dict.update({i: log10(j) - log10(out_dict.get(i))})


    # come back to this later
    if len(condensing_solids) > 0:

        solid_molecule_dict, Y, solid_element_dict, Z = get_data.get_equations(condensing_solids)
        for i, j in dict.iteritems(solid_element_dict):  # open an element
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

            final_dict.update({i: 1. - ((log10(out_dict[i] + Sum) / log10(Par_Pres.get(i))))})


    # end come back to this later

    counter = 0
    for i in Name:
        if i not in condensing_solids:  # if the molecule isn't condensed, append it to the final dict
            finals.append(final_dict.get(i))
    # condensing_solids = ['Al2O3']
    # solid_molecule_dict,Y,solid_element_dict,Z = get_data.get_equations(condensing_solids)

    # come back to this later
    if len(condensing_solids) > 0:

        for i, j in dict.iteritems(solid_molecule_dict):

            pressures = 0.
            for k in j:
                if guess_dict.get(k[0]) < 0.:
                    pressures = 1.e999
                    break
                else:
                    if k[0] in fugacities:
                        pressures += (float(k[1]) / 2.) * (log10(guess_dict.get(k[0]) * RT))
                    else:
                        pressures += (float(k[1])) * (log10(guess_dict.get(k[0]) * RT))

            if K.get(i) == 0.:
                finals.append(-K.get(i) - pressures)
            else:
                finals.append((1. - (pressures / -K.get(i))))

    # Check whether solid equilibria even matters...

    return finals
