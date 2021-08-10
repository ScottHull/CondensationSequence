from math import exp, log, log10
import re
import numpy as np
import pandas as pd
import sys

from src import collect_data


def read_janaf_file(path):
    temperature_ref = []
    delta_H_ref = []
    H_ref = []
    S_ref = []
    delta_G_ref = []
    K_ref = []
    df = open(path, "r")
    for aRow in df:
        if not aRow.startswith("#"):
            values = aRow.split('\t')
            values = [np.inf if x == 'INFINITE\n' else x for x in values]
            temperature_ref.append(float(values[0]))
            delta_H_ref.append(float(values[4]))
            H_ref.append(float(values[5]))
            S_ref.append(float(values[2]))
            delta_G_ref.append(float(values[6]))
            K_ref.append(float(values[7]))
    return temperature_ref, delta_H_ref, H_ref, S_ref, delta_G_ref, K_ref


def calc_solid_K(molecule, temperature, delta_H, S0, C0, C1, C2, C3):
    """
    Calculate the K equilibrium constant values of solid molecules based on the fitted input parameters.
    """
    K = 1
    atoms_in_molecule = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)
    reference_correction = 0.

    for atom in atoms_in_molecule:
        temperature_ref, delta_H_ref, H_ref, S_ref, delta_G_ref, K_ref = read_janaf_file(
            "data/Reference/" + atom[0] + "_ref.dat")

        H0 = collect_data.lookup_and_interpolate(temperature_ref, H_ref, 298.15) * 1000.
        change_H = collect_data.lookup_and_interpolate(temperature_ref, delta_H_ref, temperature) * 1000.
        change_S = collect_data.lookup_and_interpolate(temperature_ref, S_ref, temperature)
        fugacities = ['H', 'Cl', 'F', 'O', 'N']
        if atom[0] in fugacities:
            reference_correction += (float(atom[1]) / 2.) * (H0 + change_H - (temperature * change_S))
        else:
            reference_correction += (float(atom[1])) * (H0 + change_H - (temperature * change_S))

    del_G_formation_solid = delta_H - (temperature * S0) + C0 + C1 + C2 + C3
    del_G_formation_reaction = del_G_formation_solid - reference_correction

    # delta G = -RT ln(k)
    # ln(k) = -deltaG / RT
    # K = exp(-deltaG / RT)
    # Note: 10**((-deltaG / RT) * log10(exp(1))) = exp(-deltaG / RT)
    # Therefore, ln(K) = (-deltaG / RT) * log10(exp(1))
    K = np.exp((-del_G_formation_reaction / (8.3144621 * temperature)))  # K = exp(-deltaG/RT)

    # try:
    # K = exp((-del_G_formation_reaction / (8.3144621 * temperature)))  # K = exp(-deltaG/RT)
    # K = np.exp((-del_G_formation_reaction / (8.3144621 * temperature)))  # K = exp(-deltaG/RT)
    # except OverflowError:
    #     K = float('inf')

    return K


def solid_K_B8_B5_JB(molecule, temperature, k0, k1, k2, k3, S0, delta_H):
    C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
    C1 = 2. * k1 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
            temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
    C2 = -k2 * (((temperature ** -1.) - (298.15 ** -1.)) - (
            (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
    C3 = -k3 * ((((temperature ** -2.) - (298.15 ** -2.)) / 2.) - (
            (temperature / 3.) * (((temperature ** -3.) - (298.15 ** -3.)))))
    K = calc_solid_K(molecule=molecule, temperature=temperature, C0=C0, C1=C1, C2=C2, C3=C3, S0=S0, delta_H=delta_H)
    return K


def solid_K_B3(molecule, temperature, k0, k1, k2, k3, S0, delta_H):
    C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
    C1 = -k1 * (((temperature ** -1.) - (298.15 ** -1.)) - (
            (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
    C2 = 2. * k2 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
            temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
    C3 = k3 * ((log(temperature) - log(298.15)) + (temperature * ((temperature ** -1.) - (298.15 ** -1.))))
    K = calc_solid_K(molecule=molecule, temperature=temperature, C0=C0, C1=C1, C2=C2, C3=C3, S0=S0, delta_H=delta_H)
    return K


def solid_K_J(molecule, temperature, k0, k1, k2, k3, S0, delta_H):
    K = 1

    temperature_ref, delta_H_ref, H_ref, S_ref, delta_G_ref, K_ref = read_janaf_file("data/Janaf/" + molecule + ".dat")

    if temperature > float(temperature_ref[-1]):
        C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
        C1 = 2 * k1 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
                temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
        C2 = -k2 * (((temperature ** -1.) - (298.15 ** -1.)) - (
                (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
        C3 = -k3 * ((((temperature ** -2.) - (298.15 ** -2.)) / 2.) - (
                (temperature / 3.) * (((temperature ** -3.) - (298.15 ** -3.)))))
        K = calc_solid_K(molecule=molecule, temperature=temperature, C0=C0, C1=C1, C2=C2, C3=C3, S0=S0, delta_H=delta_H)
        return K
    else:
        # data in table is given as log10(K)
        try:
            K = 10 ** collect_data.lookup_and_interpolate(temperature_ref, K_ref, temperature)
        except OverflowError:
            K = float('inf')
        return K


def solid_K_R(molecule, temperature, k0, k1, k2, k3, S0, delta_H):
    C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
    C1 = k1 * ((0.5 * ((temperature ** 2) - (298.15 ** 2))) - (temperature * (temperature - 298.15)))
    C2 = 2. * k2 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
            temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
    C3 = -k3 * (((temperature ** -1.) - (298.15 ** -1.)) - (
            (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
    K = calc_solid_K(molecule=molecule, temperature=temperature, C0=C0, C1=C1, C2=C2, C3=C3, S0=S0, delta_H=delta_H)
    return K


def get_K_solids(temperature):
    """
    Initial function for calculating the K equilibrium constant values of solid molecules.
    """
    K_dict = {}
    molecules = pd.read_csv("data/Solids/Solids_Cp.dat", delimiter="\t", header=None, skiprows=1, index_col=False)

    for row in molecules.itertuples(index=False):
        molecule = row[0]
        method = row[1]
        delta_H = row[2]
        S0 = row[3]
        k0 = row[4]
        k1 = row[5]
        k2 = row[6]
        k3 = row[7]
        K = 0

        if method == 'B8' or method == 'B5' or method == 'JB':
            K = solid_K_B8_B5_JB(molecule=molecule, temperature=temperature, k0=k0, k1=k1, k2=k2, k3=k3, S0=S0,
                                 delta_H=delta_H)

        elif method == 'B3':
            K = solid_K_B3(molecule=molecule, temperature=temperature, k0=k0, k1=k1, k2=k2, k3=k3, S0=S0,
                           delta_H=delta_H)

        elif method == 'R':
            K = solid_K_R(molecule=molecule, temperature=temperature, k0=k0, k1=k1, k2=k2, k3=k3, S0=S0,
                          delta_H=delta_H)

        elif method == 'J':
            K = solid_K_J(molecule=molecule, temperature=temperature, k0=k0, k1=k1, k2=k2, k3=k3, S0=S0,
                          delta_H=delta_H)

        else:
            K = 1.0

        K_dict.update({molecule + "_s": K})

    return K_dict


def gas_F(temperature, A, B):
    log10K = A + (B / temperature)
    return 10 ** log10K


def gas_J(temperature, molecule, R=8.3144621):
    """
    Calculates K for gasses with the JANAF regime.
    :param temperature:
    :param molecule:
    :param R:
    :return:
    """
    temperature_ref, delta_H_ref, H_ref, S_ref, delta_G_ref, K_ref = read_janaf_file("data/Gasses/" + molecule + ".dat")
    delta_G_gas = collect_data.lookup_and_interpolate(temperature_ref, delta_G_ref,
                                                      temperature) * 1000  # interpolate deltaG
    if temperature < temperature_ref[-1]:
        # delta G = -RT ln(k)
        # ln(k) = -deltaG / RT
        # K = exp(-deltaG / RT)
        # Note: 10**((-deltaG / RT) * log10(exp(1))) = exp(-deltaG / RT)
        # Therefore, ln(K) = (-deltaG / RT) * log10(exp(1))
        K = exp(-delta_G_gas / (R * temperature))
        return K
    else:
        return 0


def gas_K(temperature, molecule, R=8.3144621):
    path = pd.read_csv("data/Gasses.dat", delimiter="\t", header=None, index_col=0)
    data = path.loc[molecule]
    H0 = float(data[2])
    S0 = float(data[3])
    A = float(data[4])
    B = float(data[5]) * 1e-3
    C = float(data[6]) * 1e6

    delH = H0 + A * (temperature - 298.15) + ((B / 2.) * ((temperature ** 2) - (298.15) ** 2)) - C * (
            (1. / temperature) - (1. / 298.15))
    S = S0 + A * (log(temperature) - log(298.15)) + B * (temperature - 298.15) - (C / 2.) * (
        ((temperature ** -2) - (298.15) ** -2))

    delG = delH - temperature * S  # deltaG = deltaH - TdeltaS
    molecule_atom_stoich = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)
    reference_correction = 0.

    for i in molecule_atom_stoich:
        Temp_ref = []
        del_H_ref = []
        S_ref = []
        G_ref = []
        H_ref = []

        path = open("data/Reference/" + i[0] + "_ref.dat", "rU")
        for aRow in path:
            if not aRow.startswith("#"):
                values = aRow.split('\t')
                Temp_ref.append(float(values[0]))
                del_H_ref.append(float(values[4]))
                H_ref.append(float(values[5]))
                S_ref.append(float(values[2]))

        H0 = collect_data.lookup_and_interpolate(Temp_ref, H_ref, 298.15) * 1000.
        change_H = collect_data.lookup_and_interpolate(Temp_ref, del_H_ref, temperature) * 1000.
        change_S = collect_data.lookup_and_interpolate(Temp_ref, S_ref, temperature)
        fugacities = ['H', 'Cl', 'F', 'O', 'N']
        if i[0] in fugacities:
            reference_correction += (float(i[1]) / 2.) * (H0 + change_H - (temperature * change_S))
        else:
            reference_correction += (float(i[1])) * (H0 + change_H - (temperature * change_S))
    del_G_corrected = delG - reference_correction
    K = exp(-del_G_corrected / (R * temperature))
    return K


def gas_P(temperature, molecule, R=8.3144621):
    path_molecule = open("data/Gasses/" + molecule + ".dat", "r")
    Temp_ref = []
    del_G_ref = []

    for aRow in path_molecule:
        values = aRow.split('\t')
        if not aRow.startswith('#'):
            Temp_ref.append(float(values[0]))
            del_G_ref.append(float(values[1]))
    G0 = collect_data.lookup_and_interpolate(Temp_ref, del_G_ref, temperature)

    del_G_gas = collect_data.lookup_and_interpolate(Temp_ref, del_G_ref, temperature)

    molecule_atom_stoich = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)
    reference_correction = 0.

    for i in molecule_atom_stoich:
        Temp_ref = []
        del_H_ref = []
        S_ref = []
        G_ref = []
        H_ref = []

        path = open("data/Reference/" + i[0] + "_ref.dat", "r")
        for aRow in path:
            if not aRow.startswith("#"):
                values = aRow.split('\t')
                Temp_ref.append(float(values[0]))
                del_H_ref.append(float(values[4]))
                H_ref.append(float(values[5]))
                S_ref.append(float(values[2]))

        H0 = collect_data.lookup_and_interpolate(Temp_ref, H_ref, 298.15) * 1000.
        change_H = collect_data.lookup_and_interpolate(Temp_ref, del_H_ref, temperature) * 1000.
        change_S = collect_data.lookup_and_interpolate(Temp_ref, S_ref, temperature)
        fugacities = ['H', 'Cl', 'F', 'O', 'N']
        if i[0] in fugacities:
            reference_correction += (float(i[1]) / 2.) * (H0 + change_H - (temperature * (change_S)))
        else:
            reference_correction += (float(i[1])) * (H0 + change_H - (temperature * (change_S)))
    del_G_corrected = del_G_gas - reference_correction
    K = exp(-del_G_corrected / (R * temperature))

    return K


def get_K_gas(temperature):
    """
    The initial function for calculating the K equilibrium constants for gasses.
    """
    gas_file = open("data/Gasses.dat", 'r')
    K_gas = {}
    for row in gas_file:
        if not row.startswith("#"):
            values = row.split('\t')
            molecule = values[0]
            method = values[1]  # define method for determining K value (i.e. J = JANAF)

            if method == "J":  # JANAF
                K = gas_J(temperature=temperature, molecule=molecule)
                K_gas.update({molecule: K})

            elif method == "K":
                K = gas_K(temperature=temperature, molecule=molecule)
                K_gas.update({molecule: K})

            elif method == "P":
                K = gas_P(temperature=temperature, molecule=molecule)
                K_gas.update({molecule: K})

            elif method == "F":
                A = float(values[2])
                B = float(values[3])
                K = gas_F(temperature=temperature, A=A, B=B)
                K_gas.update({molecule: K})

    return K_gas  # returns a dict of all K values for gasses


def liquid_J(temperature, molecule, R=8.3144621):
    """
    Calculates K for gasses with the JANAF regime.
    :param temperature:
    :param molecule:
    :param R:
    :return:
    """
    path_molecule = pd.read_csv("data/Liquids/" + molecule + ".dat", delimiter="\t", skiprows=2, header=None,
                                index_col=False)
    reference_temperature = path_molecule[0]
    delta_G_ref = path_molecule[6]
    delta_G_gas = collect_data.lookup_and_interpolate(reference_temperature, delta_G_ref,
                                                      temperature) * 1000  # interpolate deltaG
    if temperature < reference_temperature.iloc[-1]:
        # delta G = -RT ln(k)
        # ln(k) = -deltaG / RT
        # K = exp(-deltaG / RT)
        # Note: 10**((-deltaG / RT) * log10(exp(1))) = exp(-deltaG / RT)
        # Therefore, ln(K) = (-deltaG / RT) * log10(exp(1))
        K = exp(-delta_G_gas / (R * temperature))
        return K
    else:
        return 1


def liquid_F(temperature, A, B):
    return 10 ** (A + (B / temperature))


def get_K_liquid(temperature):
    """
    Initial function for calculating the K equilibrium constant values of liquid molecules.
    """
    K_dict = {}
    # molecules = pd.read_csv("data/Liquids.dat", delimiter="\t", header=None, skiprows=0, index_col=False)
    molecules = pd.read_csv("data/liquid_pseudospecies_fegley_1986.dat", delimiter="\t", header=None, skiprows=1,
                            index_col=False)

    for row in molecules.itertuples(index=False):
        molecule = row[0]
        method = row[1]

        if method == 'J':
            K = liquid_J(molecule=molecule, temperature=temperature)
        elif method == "F":
            K = liquid_F(temperature=temperature, A=row[2], B=row[3])
        else:
            K = 1

        K_dict.update({molecule + "_l": K})

    return K_dict


def get_K(gas_molecules, gas_methods, solid_molecules, liquid_molecules, temperature, solid, liquid, gas):
    """
    Returns a dictionary of all K equilibrium constants for both solid molecules (ending with _s) and gas molecules.
    """
    K_dict = {}
    if gas:
        K_gas = get_K_gas(temperature=temperature)
        K_dict.update(K_gas)
    if liquid:
        K_liquid = get_K_liquid(temperature=temperature)
        K_dict.update(K_liquid)
    if solid:
        K_solids = get_K_solids(temperature=temperature)
        K_dict.update(K_solids)
    return K_dict


def __vf_K(A, B, T):
    return 10 ** (A + (B / T))


def __get_constants_vf():
    base_liq_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 1", index_col="Reactant")
    gas_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 2", index_col="Product")
    liq_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 3", index_col="Product")
    constants = {}  # A, B in log10(K) = A + B/T
    for species in base_liq_df.index:
        A, B = base_liq_df["A"][species], base_liq_df["B"][species]
        constants.update({species + "_l": {"A": A, "B": B}})
    for species in gas_df.index:
        A, B = gas_df["A"][species], gas_df["B"][species]
        constants.update({species: {"A": A, "B": B}})
    for species in liq_df.index:
        A, B = liq_df["A"][species], liq_df["B"][species]
        constants.update({species + "_l": {"A": A, "B": B}})

    return constants


def get_K_vf(self, temperature):
    """
    K constants provided by methods of Vischer & Fegley 1987.
    :return:
    """
    base_oxides = ["SiO2_l", 'MgO', "FeO_l", "CaO_l", "Al2O3_l", "TIO2_l", "Na2O_l", "K2O_l", "ThO2_l", "UO2_l",
                   "PuO2_l"]
    self.constants = __get_constants_vf()
    K_dict = {}
    for species in self.constants.keys():
        K_dict.update(({
            species: __vf_K(A=self.constants[species]["A"], B=self.constants[species]["B"], T=temperature)
        }))
    print(K_dict)
