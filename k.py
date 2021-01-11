from math import exp, log, log10
import pandas as pd
import collect_data


def get_K_solids(temperature):
    molecules = pd.read_csv("Data/Solids/Solids_Cp.dat", header=None, skiprows=1)
    molecule = molecules[0]
    method = molecules[1]
    deltaH = molecule[2]
    S0 = molecules[3]
    k0 = molecules[4]
    k1 = molecules[5]
    k2 = molecules[6]
    k3 = molecules[7]
    C0 = 0
    C1 = 0
    C2 = 0
    C3 = 0

    for index, m in enumerate(molecule):
        if method[index] == 'B8' or method[index] == 'B5' or method[index] == 'JB':
            C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
            C1 = 2. * k1 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
                        temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
            C2 = -k2 * (((temperature ** -1.) - (298.15 ** -1.)) - (
                        (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
            C3 = -k3 * ((((temperature ** -2.) - (298.15 ** -2.)) / 2.) - (
                        (temperature / 3.) * (((temperature ** -3.) - (298.15 ** -3.)))))

        if method[index] == 'B3':
            C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
            C1 = -k1 * (((temperature ** -1.) - (298.15 ** -1.)) - (
                        (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))
            C2 = 2. * k2 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
                        temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
            C3 = k3 * ((log(temperature) - log(298.15)) + (temperature * ((temperature ** -1.) - (298.15 ** -1.))))

        if method[index] == 'R':
            C0 = k0 * ((temperature - 298.15) - (temperature * (log(temperature) - log(298.15))))
            C1 = k1 * ((0.5 * ((temperature ** 2) - (298.15 ** 2))) - (temperature * (temperature - 298.15)))
            C2 = 2. * k2 * (((temperature ** 0.5) - (298.15 ** 0.5)) + (
                        temperature * (((temperature ** -0.5)) - ((298.15 ** -0.5)))))
            C3 = -k3 * (((temperature ** -1.) - (298.15 ** -1.)) - (
                        (temperature / 2.) * (((temperature ** -2.) - (298.15 ** -2.)))))

        if method[index] == 'J':
            molecule_df = pd.read_csv("Data/Janaf/" + molecule[index] + ".dat", header=None, skiprows=2)
            temperature_ref = molecule_df[0]
            delta_H_ref = molecule_df[4]
            H_ref = molecule_df[5]
            S_ref = molecule_df[2]
            delta_G_ref = molecule_df[6]
            K_ref = molecule_df[7]

            if temperature > temperature_ref[-1]:
                del_H = float(molecule[2])
                S0 = float(molecule[3])
                k0 = float(molecule[4])
                k1 = float(molecule[5])
                k2 = float(molecule[6])
                k3 = float(molecule[7])
                C0 = k0 * ((T - 298.15) - (T * (log(T) - log(298.15))))
                C1 = 2. * k1 * (((T ** 0.5) - (298.15 ** 0.5)) + (T * (((T ** -0.5)) - ((298.15 ** -0.5)))))
                C2 = -k2 * (((T ** -1.) - (298.15 ** -1.)) - ((T / 2.) * (((T ** -2.) - (298.15 ** -2.)))))
                C3 = -k3 * ((((T ** -2.) - (298.15 ** -2.)) / 2.) - ((T / 3.) * (((T ** -3.) - (298.15 ** -3.)))))
            else:
                K = lookup_and_interpolate(Temp_ref, K_ref, T)
                return K

        Mol_name_split = re.findall(r'([A-Z][a-z]*)(\d*)', Name)
        reference_correction = 0.

        for i in Mol_name_split:
            Temp_ref = []
            del_H_ref = []
            S_ref = []
            G_ref = []
            H_ref = []

            path = open("Data/Reference/" + i[0] + "_ref.dat", "rU")
            for aRow in path:
                if not aRow.startswith("#"):
                    values = aRow.split('\t')
                    Temp_ref.append(float(values[0]))
                    del_H_ref.append(float(values[4]))
                    H_ref.append(float(values[5]))
                    S_ref.append(float(values[2]))

            H0 = lookup_and_interpolate(Temp_ref, H_ref, 298.15) * 1000.
            change_H = lookup_and_interpolate(Temp_ref, del_H_ref, T) * 1000.
            change_S = lookup_and_interpolate(Temp_ref, S_ref, T)
            # print "change_S",i[0],i[1],change_S
            if i[0] in fugacities:
                reference_correction += (float(i[1]) / 2.) * (H0 + change_H - (T * (change_S)))
            else:
                reference_correction += (float(i[1])) * (H0 + change_H - (T * (change_S)))

        del_G_formation_solid = del_H - (T * S0) + C0 + C1 + C2 + C3

        # print "ref",reference_correction/1e3
        # print "solid",del_G_formation_solid/1e3
        # print del_G_formation_solid/1e3
        # print reference_correction/1e3

        del_G_formation_reaction = del_G_formation_solid - reference_correction

        K = log10(exp(1)) * ((-del_G_formation_reaction / (8.3144621 * T)))

        return K


def calculate_K_gas(molecules, temperature):
    K_gas = {}
    for m in molecules:
        path_molecule = pd.read_csv("Data/Gasses/" + m + ".dat", delimiter="\t", skiprows=1, header=None)
        reference_temperature = path_molecule[0]
        delta_H = path_molecule[1]
        H_ref = path_molecule[4]
        S_ref = path_molecule[2]
        delta_G_ref = path_molecule[6]
        delta_G_gas = collect_data.lookup_and_interpolate(reference_temperature, delta_G_ref,
                                                          temperature) * 1000  # interpolate deltaG
        K = 0
        if temperature < reference_temperature[-1]:
            # delta G = -RT ln(k)
            # ln(k) = -deltaG / RT
            # K = exp(-deltaG / RT)
            K = exp(-delta_G_gas / (8.3144621 * temperature))
            return K
        K_gas.update({m: K})
    return K_gas


def get_K(gas_molecules, solid_molecules, temperature):
    K_dict = {}
    K_gas = calculate_K_gas(molecules=gas_molecules.keys(), temperature=temperature)
    K_solids = get_K_data(i, T)
    K_dict = K_gas.copy()
    K_dict.update(K_solids)
    return K_dict
