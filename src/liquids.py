from math import log10
from scipy.optimize import brentq
from src import collect_data
import sys


def check_in(liquids, number_densities, temperature, K_dict, condensing_liquids, temperature_old, K_dict_old,
             number_densities_old, removed_liquids, is_liquid, R=8.3144621e-2):
    """

    :param liquids: A dictionary of liquid molecules and their stoichiometry
    :param n_x: A dictionary of elements and their number densities
    :param T: Current temperature
    :param K_dict: A dictionary of all K values
    :param condensing_liquids: A list of all currently condensed liquids in the system
    :param T_old:
    :param K_dict_old:
    :param n_x_old:
    :param removed_liquids_old:
    :param removed_liquids:
    :return:
    """

    fugacities = ['H', 'Cl', 'F', 'O', 'N']  # diatomic molecules
    if not is_liquid:
        return False, False

    new_liquids = []
    for liquid_molecule in liquids.keys():
        pressure_product = 1  # the product of all reactant partial pressures
        pressure_product_old = 1
        liquid_molecule_stoich = liquids[liquid_molecule]
        for element in liquid_molecule_stoich.keys():
            stoich = liquid_molecule_stoich[element]
            if element in fugacities:
                pressure_product *= (number_densities[element] * R * temperature) ** (stoich / 2.0)
                if element in number_densities_old.keys():
                    pressure_product_old *= (number_densities_old[element] * R * temperature_old) ** (stoich / 2.0)
            else:
                pressure_product *= (number_densities[element] * R * temperature) ** (stoich)
                if element in number_densities_old.keys():
                    pressure_product_old *= (number_densities_old[element] * R * temperature_old) ** (stoich)

        # K * prod(P_i) = P_mol.  If P_mol >= 1, then the liquid is stable (i.e. log10(1 / KP) < 0)
        # - log10(K) - log10(prod(P_i)) = log(1 / K prod(P_i)) -> raise to power of 10 -> 1 / K prod(P_I)
        condensation_criterion = 1.0 / (
                K_dict[liquid_molecule] * pressure_product)  # i.e. "the chemical activity exceeds 1"
        if condensation_criterion < 1 and temperature != temperature_old:
            # if the partial pressure exceeds the equilibrium constant
            if liquid_molecule not in condensing_liquids:  # if the condensation criterion has previously been met
                condensation_criterion_old = 1.0 / (
                        K_dict_old[liquid_molecule] * pressure_product_old)  # i.e. "the chemical activity exceeds 1"
                if condensation_criterion_old < 1:
                    new_liquids.append([liquid_molecule, temperature_old])
                else:  # if the condensation criterion has just been met
                    # interpolate the exact condensation temperature
                    # f(a) and f(b) must have different signs, so we employ log10 for a change between <1 and >1 (i.e. interpolate around log10(1) = 0)
                    test_temperature = brentq(collect_data.linear_interpol_in, temperature, temperature_old, args=(
                        temperature, temperature_old, log10(condensation_criterion), log10(condensation_criterion_old)),
                                              xtol=0.0000000001)
                    new_liquids.append([liquid_molecule, test_temperature])

    if len(new_liquids) > 0:  # if there are any new liquids
        liquid_molecule = False
        test_temperature = 0  # the maximum interpolated temperature from the loop below
        for liquid in new_liquids:  # loop through liquids
            liquid_molecule = liquid[0]  # the name of the liquid molecule
            liquid_condensation_temperature = liquid[1]  # the temperature of condensation
            if liquid_condensation_temperature > test_temperature:
                if liquid_condensation_temperature < temperature_old:
                    test_temperature = liquid_condensation_temperature
        if test_temperature != 0:
            return liquid_molecule, test_temperature
        else:
            return False, False
    else:
        return False, False


def check_out(condensing_liquids, number_density_liquids, number_density_liquids_old, temperature, temperature_old):
    temperature_out = {}
    dT = abs(temperature_old - temperature)  # absolute temperature difference
    for liquid in condensing_liquids:  # loop through condensing liquids list
        new_liquid = number_density_liquids[liquid]  # current number density
        old_liquid = number_density_liquids_old[liquid]  # previous number density
        if old_liquid is None:
            old_liquid = -1.

        # if the liquid's number density falls below the 1^-10 within the iteration, remove it from the system
        if new_liquid <= 1.e-10 and old_liquid > 1.e-10:  # if the new liquid meets the mole fraction criterion for existance
            out_temp = brentq(collect_data.linear_interpol_out, temperature, temperature_old, args=(
                temperature, temperature_old, log10(number_density_liquids[liquid]),
                log10(number_density_liquids_old[liquid]), dT),
                              xtol=0.0000000001)
            temperature_out.update({liquid: out_temp})

    for liquid in temperature_out.keys():
        out_temp = temperature_out[liquid]
        if out_temp > 0:
            return liquid, out_temp
        else:
            return False, False
    return False, False
