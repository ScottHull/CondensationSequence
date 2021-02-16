import pandas as pd
import re
import bisect
import sys


def get_methods(path):
    """
    Returns a dictionary of the method for calculating K for each molecule.
    :param path:
    :return:
    """
    methods = {}
    df = pd.read_csv(path, delimiter="\t", header=None)
    for row in df.index:
        methods.update({df[0][row]: df[1][row]})
    return methods


def linear_interpol(x, x1, x2, y1, y2):
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    assert (x1 <= x)
    assert (x2 >= x)
    assert (x1 <= x2)

    alpha = (x - x1) / (x2 - x1)
    return (1.0 - alpha) * y1 + alpha * y2


def linear_interpol_out(x, *args):
    x1 = args[0]
    x2 = args[1]
    y1 = args[2]
    y2 = args[3]
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    # print "x1",x1
    # print "x2",x2
    # print "x",x
    assert (x1 <= x)
    assert (x2 >= x)
    assert (x1 <= x2)

    alpha = (x - x1) / (x2 - x1)
    return -10. - ((1. - alpha) * y1 + alpha * y2)


def lookup_and_interpolate(table_x, table_y, x_value):
    table_x = list(table_x)
    table_y = list(table_y)
    idx = bisect.bisect_left(list(table_x), x_value) - 1
    if idx < 0:
        return table_y[0]
    elif (idx < len(table_x) - 1):
        return linear_interpol(x_value, table_x[idx], table_x[idx + 1], table_y[idx], table_y[idx + 1])
    else:
        return table_y[idx]


def get_atoms_from_molecule(path, skiprows=0, solid=False, liquid=False):
    """
    Returns the stoichiometry of molecules as a dictionary.
    :param path:
    :param skiprows:
    :return:
    """
    molecules_dict = {}
    df = pd.read_csv(path, delimiter="\t", header=None, skiprows=skiprows)
    molecules = df[0]
    for m in molecules:
        molecule_name = m
        if solid:
            molecule_name = m + "_s"
        if liquid:
            molecule_name = m + "_l"
        molecules_dict.update({molecule_name: {}})
        atom_nums = re.findall(r'([A-Z][a-z]*)(\d*)', m)
        for i in atom_nums:
            name = i[0]
            molecules_dict[molecule_name].update({name: int(i[1])})
    return molecules_dict


def element_appearances_in_molecules(abundances, library):
    """
    Given the abundances dictionary, return a dictionary where the element is the key and the molecules it appears in
    as a list is the value.
    :param abundances:
    :param library:
    :return:
    """
    element_appearances = {}
    abundance_elements = abundances.keys()  # get a list of input elements
    for atom in abundance_elements:  # construct the empty return dictionary
        element_appearances.update({atom: []})
    for molecule in library.keys():
        atom_nums = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)  # returns something like [('Cl', '2'), ('Cr', '1')]
        for atom_stoich in atom_nums:  # get the atom:stoich tuple
            atom = atom_stoich[0]  # the atom name
            if atom in element_appearances:
                element_appearances[atom].append(molecule)
    return element_appearances


def linear_interpol_in(x, *args):
    x1 = args[0]
    x2 = args[1]
    y1 = args[2]
    y2 = args[3]
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    assert (x1 <= x)
    assert (x2 >= x)
    assert (x1 <= x2)

    alpha = (x - x1) / (x2 - x1)
    return (1.0 - alpha) * y1 + alpha * y2
