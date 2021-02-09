from math import log10
from scipy.optimize import brentq
import collect_data
import sys


def check_in(solids, number_densities, temperature, K_dict, condensing_solids, temperature_old, K_dict_old,
             number_densities_old, removed_solids_old, removed_solids, R=8.3144621e-2):
    """

    :param solids: A dictionary of solid molecules and their stoichiometry
    :param n_x: A dictionary of elements and their number densities
    :param T: Current temperature
    :param K_dict: A dictionary of all K values
    :param condensing_solids: A list of all currently condensed solids in the system
    :param T_old:
    :param K_dict_old:
    :param n_x_old:
    :param removed_solids_old:
    :param removed_solids:
    :return:
    """

    fugacities = ['H', 'Cl', 'F', 'O', 'N']  # diatomic molecules

    new_solids = []
    for solid_molecule in solids.keys():
        pressure_product = 1  # the product of all reactant partial pressures
        pressure_product_old = 1
        solid_molecule_stoich = solids[solid_molecule]
        for element in solid_molecule_stoich.keys():
            stoich = solid_molecule_stoich[element]
            if element in fugacities:
                pressure_product *= (number_densities[element] * R * temperature) ** (stoich / 2.0)
                if element in number_densities_old.keys():
                    pressure_product_old *= (number_densities_old[element] * R * temperature_old) ** (stoich / 2.0)
            else:
                pressure_product *= (number_densities[element] * R * temperature) ** (stoich)
                if element in number_densities_old.keys():
                    pressure_product_old *= (number_densities_old[element] * R * temperature_old) ** (stoich)

        # K * prod(P_i) = P_mol.  If P_mol >= 1, then the solid is stable (i.e. log10(1 / KP) < 0)
        # - log10(K) - log10(prod(P_i)) = log(1 / K prod(P_i)) -> raise to power of 10 -> 1 / K prod(P_I)
        condensation_criterion = 1.0 / (K_dict[solid_molecule] * pressure_product)  # i.e. "the chemical activity exceeds 1"
        if condensation_criterion < 1 and temperature != temperature_old:
            # if the partial pressure exceeds the equilibrium constant
            if solid_molecule not in condensing_solids:  # if the condensation criterion has previously been met
                condensation_criterion_old = 1.0 / (
                            K_dict_old[solid_molecule] * pressure_product_old)  # i.e. "the chemical activity exceeds 1"
                if condensation_criterion_old < 1:
                    new_solids.append([solid_molecule, temperature_old])
                else:  # if the condensation criterion has just been met
                    # interpolate the exact condensation temperature
                    # f(a) and f(b) must have different signs, so we employ log10 for a change between <1 and >1 (i.e. interpolate around log10(1) = 0)
                    test_temperature = brentq(collect_data.linear_interpol_in, temperature, temperature_old, args=(
                        temperature, temperature_old, log10(condensation_criterion), log10(condensation_criterion_old)),
                                              xtol=0.0000000001)
                    new_solids.append([solid_molecule, test_temperature])

    if len(new_solids) > 0:  # if there are any new solids
        solid_molecule = False
        test_temperature = 0  # the maximum interpolated temperature from the loop below
        for solid in new_solids:  # loop through solids
            solid_molecule = solid[0]  # the name of the solid molecule
            solid_condensation_temperature = solid[1]  # the temperature of condensation
            if solid_condensation_temperature > test_temperature:
                if solid_condensation_temperature < temperature_old:
                    test_temperature = solid_condensation_temperature
        if test_temperature != 0:
            return solid_molecule, test_temperature
        else:
            return False, False
    else:
        return False, False


def check_out(condensing_solids, number_density_solids, number_density_solids_old, temperature, temperature_old):
    temperature_out = {}
    dT = abs(temperature_old - temperature)  # absolute temperature difference
    for solid in condensing_solids:  # loop through condensing solids list
        new_solid = number_density_solids[solid]  # current number density
        old_solid = number_density_solids_old[solid]  # previous number density
        if old_solid is None:
            old_solid = -1.

        # if the solid's number density falls below the 1^-10 within the iteration, remove it from the system
        if new_solid <= 1.e-10 and old_solid >= 1.e-10:  # if the new solid meets the mole fraction criterion for existance
            out_temp = brentq(collect_data.linear_interpol_out, temperature, temperature_old, args=(
                temperature, temperature_old, log10(number_density_solids[solid]),
                log10(number_density_solids_old[solid]), dT),
                              xtol=0.0000000001)
            temperature_out.update({solid: out_temp})

    for solid in temperature_out.keys():
        out_temp = temperature_out[solid]
        if out_temp > 0:
            return solid, out_temp
        else:
            return False, False
    return False, False
