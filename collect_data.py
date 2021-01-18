import pandas as pd
import re
import bisect

def get_methods(path):
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


def lookup_and_interpolate(table_x, table_y, x_value):
    idx = bisect.bisect_left(table_x, x_value) - 1
    if idx < 0:
        return table_y[0]
    elif (idx < len(table_x) - 1):
        return linear_interpol(x_value, table_x[idx], table_x[idx + 1], table_y[idx], table_y[idx + 1])
    else:
        return table_y[idx]


def get_atoms_from_molecule(path, skiprows=0):
    molecules_dict = {}
    df = pd.read_csv(path, delimiter="\t", header=None, skiprows=skiprows)
    molecules = df[0]
    for m in molecules:
        molecules_dict.update({m: {}})
        atom_nums = re.findall(r'([A-Z][a-z]*)(\d*)', m)
        for i in atom_nums:
            molecules_dict[m].update({i[0]: i[1]})
    return molecules_dict
