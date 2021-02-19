import copy
from math import sqrt
import sys

import numpy as np
from scipy.optimize import root

from src import collect_data
from src import k
from src import liquids
from src import mass_balance
from src import solids
from src import total_number


class Condensation:

    def __init__(self, start_temperature, end_temperature, abundances, total_pressure, dT=2,
                 gas=True, liquid=True, solid=True):
        self.IS_GAS = gas  # if we are to equilibrate solids
        self.IS_LIQUID = liquid  # if we are to equilibrate liquids
        self.IS_SOLID = solid  # if we are to equilibrate solids

        self.R = 8.3144621e-2
        self.K = {}
        self.start_temperature = start_temperature
        self.temperature = start_temperature
        self.end_temperature = end_temperature
        self.dT = dT
        self.initial = True
        self.previous_temperature = start_temperature
        self.total_pressure = total_pressure
        self.gas_methods = collect_data.get_methods(
            path="data/Gasses.dat")  # methods for calculating K for gasses (i.e. JANAF)
        self.gas_molecules_library = collect_data.get_atoms_from_molecule(path="data/Gasses.dat")
        self.liquid_molecules_library = collect_data.get_atoms_from_molecule(path="data/Liquids.dat", liquid=True)
        self.solid_molecules_library = collect_data.get_atoms_from_molecule(path="data/Solids_Cp.dat", skiprows=1,
                                                                            solid=True)
        self.abundances = abundances
        self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
        self.element_gas_appearances = collect_data.element_appearances_in_molecules(abundances=self.abundances,
                                                                                     library=self.gas_molecules_library)
        self.element_liquid_appearances = collect_data.element_appearances_in_molecules(abundances=self.abundances,
                                                                                        library=self.liquid_molecules_library)
        self.element_solid_appearances = collect_data.element_appearances_in_molecules(abundances=self.abundances,
                                                                                       library=self.solid_molecules_library)

        self.mass_balance = self.calculate_mass_balance()
        self.condensing_solids = []
        self.condensing_liquids = []
        self.removed_solids = []
        self.removed_liquids = []
        self.previous_removed_solids = []
        self.previous_removed_liquids = []
        self.K = {}  # a dictionary of all K equilibrium constant values
        self.previous_K = {}  # tracks K values from the previous temperature iteration
        self.number_densities = {}  # element number densities
        self.previous_number_densities = {}  # tracks number density solution from the previous temperature iteration
        self.number_densities_solids = {}  # number densities of all solids
        self.previous_number_densities_solids = {}  # number densities of all solids
        self.number_densities_liquids = {}  # number densities of all liquids
        self.previous_number_densities_liquids = {}  # number densities of all liquids
        self.elements_in_solid = []  # all elements currently occupying a condensing solid
        self.elements_in_liquid = []  # all elements currently occupying a condensing liquid
        self.total_elements_condensed = {}  # tracks the number density of all elements in the condensed phase
        self.percent_element_condensed = {}
        for element in self.abundances.keys():  # initial setup
            self.percent_element_condensed.update({element: []})
            self.total_elements_condensed.update({element: 0})
        self.total_number = 0

        self.any_in = True
        self.any_out = True
        self.errors = []
        self.error_threshold = 0

        self.tracked_solids = {}  # tracks in and out temperatures of the solid
        self.tracked_liquids = {}  # tracks in and out temperatures of the liquid

    def solve(self):
        number_densities_from_partial_pressure = self.calculate_mass_balance()
        # list of elements given in the input
        names = list(self.number_densities.keys())
        guess_number_density = np.array(
            [self.number_densities[p] for p in names])  # a list of abundances corresponding to the elements list
        if not self.initial:
            # a list of abundances corresponding to the elements list
            for s in self.condensing_solids:
                if s not in names:
                    names.append(s)
                    guess_number_density = np.append(guess_number_density,
                                                     1.e-13)  # append a guess for all existing solids
            for s in self.condensing_liquids:
                if s not in names:
                    names.append(s)
                    guess_number_density = np.append(guess_number_density,
                                                     1.e-13)  # append a guess for all existing solids
        self.initial = False
        args = (names, self.element_gas_appearances, self.gas_molecules_library, self.K,
                number_densities_from_partial_pressure,
                self.temperature, self.condensing_solids, self.solid_molecules_library, self.condensing_liquids,
                self.liquid_molecules_library)
        # finds the root of the mass balance equation where the root is the number density (initial guess is just the partial pressures)
        print("Solving system...")
        number_densities = root(mass_balance.mass_balance, guess_number_density, args=args, method='lm',
                                options={'maxiter': 100000000,
                                         'ftol': 1.e-15})  # calculate the activities of gas phases corresponding to the element dict
        self.number_densities = dict(zip(names,
                                         number_densities.x))  # make a dictionary where the element is the key and the activity is the value
        self.errors = mass_balance.mass_balance(number_densities.x, *args)
        self.error_threshold = sqrt(
            sum([i ** 2 for index, i in enumerate(self.errors) if names[index] not in self.condensing_solids]))
        self.errors = dict(zip([n for n in names if n not in self.condensing_solids and n
                                not in self.condensing_liquids], self.errors))
        print("Solved system!")

    def normalize_abundances(self, abundances):
        # norm = {}
        # total = sum(abundances.values())
        # for element in abundances.keys():
        #     norm.update({element: self.abundances[element] / total})
        # return norm

        H = self.abundances['H']
        He = self.abundances['He']
        Ne = self.abundances['Ne']
        Ar = self.abundances['Ar']
        X = 0
        if self.initial:
            K = k.get_K_gas(molecules=self.gas_molecules_library.keys(), methods=self.gas_methods,
                            temperature=self.temperature)
            X = K["H1"] ** 2
        else:
            X = self.K["H1"] ** 2
        H2_coef_a = 4. + X
        H2_coef_b = -((4. * H) + (X * H))
        H2_coef_c = H ** 2

        H2 = ((-H2_coef_b) - np.sqrt(pow(H2_coef_b, 2.) - (4. * H2_coef_a * H2_coef_c))) / (2. * H2_coef_a)
        monoH = H - 2. * H2
        total = monoH + H2 + He + Ne + Ar
        normalized_abundances = {}
        for element in self.abundances.keys():
            normalized_abundances.update({element: self.abundances[element] / total})
        return normalized_abundances

    def calculate_mass_balance(self):
        """
        This is equation 5 from Unterborn and Panero 2017, where N_i is proportional to the partial pressure of element i.
        N_i = (a_i / a_total) * (P_total / RT)
        :return:
        """
        mass_balance_from_pp = {}
        for element in self.normalized_abundances.keys():
            N_i = self.normalized_abundances[element] * (self.total_pressure / (self.R * self.temperature))
            mass_balance_from_pp.update({element: N_i})
        return mass_balance_from_pp

    def equilibrate_solids(self):
        """
        solids.check_in: returns the incoming solid and the interpolated temperature at which it appears.  If no new
            solid appears, it will return False and self.any_in will be set to True, thereby breaking the while loop
            outside of this function.
        solids.check_out: Will trigger the removal of a solid from the system if it falls below the threshold number
            density, thereby also returning its removal temperature.
        """

        if self.IS_SOLID:
            print("Equilibrating potential solids...")

        # in_solid is the solid that condenses
        # in_solid_temperature is the temperature at which the solid condenses
        in_solid, in_solid_temp = solids.check_in(
            solids=self.solid_molecules_library,
            number_densities=self.number_densities,
            temperature=self.temperature,
            K_dict=self.K,
            condensing_solids=self.condensing_solids,
            temperature_old=self.previous_temperature,
            K_dict_old=self.previous_K,
            number_densities_old=self.previous_number_densities,
            removed_solids=self.removed_solids,
            removed_solids_old=self.previous_removed_solids,
            is_solid=self.IS_SOLID
        )

        if self.error_threshold > 0.05:  # if the error is smaller than the threshold
            in_solid = False
            in_solid_temp = 0
        if not in_solid:  # if there is no new solid entering, break the loop
            self.any_in = False  # break out of check-in loop
            in_solid_temp = 0

        out_solid, out_solid_temp = solids.check_out(condensing_solids=self.condensing_solids,
                                                     number_density_solids=self.number_densities_solids,
                                                     number_density_solids_old=self.previous_number_densities_solids,
                                                     temperature=self.temperature,
                                                     temperature_old=self.previous_temperature)

        if not out_solid:
            self.any_out = False
            out_solid_temp = 0

        return in_solid, in_solid_temp, out_solid, out_solid_temp

    def equilibrate_liquids(self):
        """
        liquids.check_in: returns the incoming liquid and the interpolated temperature at which it appears.  If no new
            liquid appears, it will return False and self.any_in will be set to True, thereby breaking the while loop
            outside of this function.
        liquids.check_out: Will trigger the removal of a liquid from the system if it falls below the threshold number
            density, thereby also returning its removal temperature.
        """
        if self.IS_LIQUID:
            print("Equilibrating potential liquids...")

        # in_liquid is the liquid that condenses
        # in_solid_temperature is the temperature at which the liquid condenses
        in_liquid, in_liquid_temp = liquids.check_in(
            liquids=self.liquid_molecules_library,
            number_densities=self.number_densities,
            temperature=self.temperature,
            K_dict=self.K,
            condensing_liquids=self.condensing_liquids,
            temperature_old=self.previous_temperature,
            K_dict_old=self.previous_K,
            number_densities_old=self.previous_number_densities,
            removed_liquids=self.removed_liquids,
            removed_liquids_old=self.previous_removed_liquids,
            is_liquid=self.IS_LIQUID
        )

        if self.error_threshold > 1 * 10 ** -13:  # if the error is smaller than the threshold
            in_liquid = False
            in_liquid_temp = 0
        if not in_liquid:  # if there is no new liquid entering, break the loop
            self.any_in = False  # break out of check-in loop
            in_liquid_temp = 0

        out_liquid, out_liquid_temp = liquids.check_out(condensing_liquids=self.condensing_liquids,
                                                        number_density_liquids=self.number_densities_liquids,
                                                        number_density_liquids_old=self.previous_number_densities_liquids,
                                                        temperature=self.temperature,
                                                        temperature_old=self.previous_temperature)

        if not out_liquid:
            self.any_out = False
            out_liquid_temp = 0

        return in_liquid, in_liquid_temp, out_liquid, out_liquid_temp

    def calculate_total_number(self):
        # the sum of all number densities across gasses and solids
        self.total_number = total_number.calculate_total_N(
            gas_element_appearances_in_molecules=self.element_gas_appearances,  # all active gas molecules in system
            solid_element_appearances_in_molecules=self.element_solid_appearances,
            liquid_element_appearances_in_molecules=self.element_liquid_appearances,
            # all active solid molecules in system
            element_number_densities=self.number_densities,  # computed elemental number densities from root finder
            condensing_solids=self.condensing_solids,  # list of all condensing solids
            gas_molecule_library=self.gas_molecules_library,  # stoich library for all supported gas molecules
            solid_molecule_library=self.solid_molecules_library,  # stoich library for all supported solid molecules
            liquid_molecule_library=self.liquid_molecules_library,
            condensing_liquids=self.condensing_liquids,
            K_dict=self.K,  # all equilibrium constants
            temperature=self.temperature,
        )

        for i in self.condensing_solids:
            self.number_densities_solids.update({i: self.number_densities[i] / self.total_number})

        for i in self.condensing_liquids:
            self.number_densities_liquids.update({i: self.number_densities[i] / self.total_number})

    def calculate_percentage_condensed(self):
        # get the total number density of the condensed elements
        for s in self.condensing_solids:
            solid_stoich = self.solid_molecules_library[s]
            for element in solid_stoich:
                stoich = solid_stoich[element]
                element_number_density = self.number_densities[element]
                self.total_elements_condensed[element] += stoich * element_number_density
        for s in self.condensing_liquids:
            liquid_stoich = self.liquid_molecules_library[s]
            for element in liquid_stoich:
                stoich = liquid_stoich[element]
                element_number_density = self.number_densities[element]
                self.total_elements_condensed[element] += stoich * element_number_density

        for element in self.total_elements_condensed:
            element_number_density = self.total_elements_condensed[element]
            percent_condensed = element_number_density / self.total_number * 100.0
            self.percent_element_condensed[element].append({
                "temperature": self.temperature,
                "percent": percent_condensed
            })

        if self.percent_element_condensed["Si"][-1]["percent"] > 99.99 and self.percent_element_condensed["Mg"][-1]["percent"] > 99.99 and self.percent_element_condensed["Fe"][-1]["percent"] > 99.99:
            print("[*] Refractory elements exhausted!")
            print(self.percent_element_condensed["Fe"])
            sys.exit()

    def sequence(self):
        """
        Partial pressure: (n_i / n_total) * (P_total / RT)

        K : exp(-deltaG_0 / RT) (depends on given method, this is for J)

        Get Equations:
            1. Take a list of condensing solid molecules
            2. For molecule in (1), calculate the stoichiometry (i.e. H20: 2H, 1O)
            3.

        Mass Balance: If the element is NOT in list fugacities,
        :return:
        """
        print("Initial setup...")
        self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
        self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
                         temperature=self.temperature, gas_methods=self.gas_methods,
                         liquid_molecules=self.liquid_molecules_library, gas=self.IS_GAS, liquid=self.IS_LIQUID,
                         solid=self.IS_SOLID)
        self.mass_balance = self.calculate_mass_balance()

        # a list of abundances corresponding to the elements list
        self.number_densities = self.mass_balance

        while self.temperature > self.end_temperature:
            print("AT TEMPERATURE: {}".format(self.temperature))

            self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
            self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
                             temperature=self.temperature, gas_methods=self.gas_methods,
                             liquid_molecules=self.liquid_molecules_library, gas=self.IS_GAS, liquid=self.IS_LIQUID,
                             solid=self.IS_SOLID)

            self.solve()  # run the root solver

            # calculate solid number densities
            self.calculate_total_number()

            # check in stable solid molecules into the system
            while self.any_in and self.any_out:  # enter a while loop

                in_solid, in_solid_temp, out_solid, out_solid_temp = self.equilibrate_solids()
                in_liquid, in_liquid_temp, out_liquid, out_liquid_temp = self.equilibrate_liquids()

                # out_solid_temp and in_solid_temp are 0 if an actual out_solid_temp or in_solid_temp are calculated.
                if in_solid_temp > out_solid_temp:  # if the appearance temperature is greater than the disappearance temperature
                    print("IN SOLID: {} ({} K)".format(in_solid, in_solid_temp))
                    if in_solid not in self.tracked_solids.keys():
                        self.tracked_solids.update(
                            {in_solid: {"In Temperature": in_solid_temp, "Out Temperature": None}})
                    self.temperature = in_solid_temp
                    self.condensing_solids.append(in_solid)
                    self.any_out = False  # a solid has not dropped out
                    self.any_in = True  # there is a new solid in
                    # track element in solid phase if it not already is being done
                    for s in self.solid_molecules_library[in_solid]:
                        if s not in self.elements_in_solid:
                            self.elements_in_solid.append(s)
                elif out_solid_temp > in_solid_temp:  # if the disappearance temperature is greater than the appearance temperature
                    print("OUT SOLID: {} ({} K)".format(out_solid, out_solid_temp))
                    self.tracked_solids[out_solid]["Out Temperature"] = out_solid_temp
                    self.temperature = out_solid_temp
                    self.any_out = True  # there are exiting solids
                    self.any_in = False  # there are no new solids
                    self.temperature = out_solid_temp  # set the temperature to the solid-out temperature
                    self.condensing_solids.remove(
                        out_solid)  # remove the outgoing solid from the condensing solids list
                    del self.number_densities_solids[
                        out_solid]  # remove the outgoing solid from the solids number density dict
                    del self.number_densities[out_solid]

                    self.removed_solids.append(out_solid)  # track the exit of the solid

                # out_liquid_temp and in_liquid_temp are 0 if an actual out_liquid_temp or in_liquid_temp are calculated.
                if in_liquid_temp > out_liquid_temp:  # if the appearance temperature is greater than the disappearance temperature
                    print("IN liquid: {} ({} K)".format(in_liquid, in_liquid_temp))
                    if in_liquid not in self.tracked_liquids.keys():
                        self.tracked_liquids.update(
                            {in_liquid: {"In Temperature": in_liquid_temp, "Out Temperature": None}})
                    self.temperature = in_liquid_temp
                    self.condensing_liquids.append(in_liquid)
                    self.any_out = False  # a liquid has not dropped out
                    self.any_in = True  # there is a new liquid in
                    # track element in liquid phase if it not already is being done
                    for s in self.liquid_molecules_library[in_liquid]:
                        if s not in self.elements_in_liquid:
                            self.elements_in_liquid.append(s)
                elif out_liquid_temp > in_liquid_temp:  # if the disappearance temperature is greater than the appearance temperature
                    print("OUT liquid: {} ({} K)".format(out_liquid, out_liquid_temp))
                    self.tracked_liquids[out_liquid]["Out Temperature"] = out_liquid_temp
                    self.temperature = out_liquid_temp
                    self.any_out = True  # there are exiting liquids
                    self.any_in = False  # there are no new liquids
                    self.temperature = out_liquid_temp  # set the temperature to the liquid-out temperature
                    self.condensing_liquids.remove(
                        out_liquid)  # remove the outgoing liquid from the condensing liquids list
                    del self.number_densities_liquids[
                        out_liquid]  # remove the outgoing liquid from the liquids number density dict
                    del self.number_densities[out_liquid]

                    self.removed_liquids.append(out_liquid)  # track the exit of the liquid

            if self.any_out or self.any_in:  # if there are any solids introduced or dropped by the above loop, then recalculate the number densities and abundances
                self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
                self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
                                 temperature=self.temperature, gas_methods=self.gas_methods,
                                 liquid_molecules=self.liquid_molecules_library, gas=self.IS_GAS, liquid=self.IS_LIQUID,
                                 solid=self.IS_SOLID)

                # resolve the system for new solids
                self.mass_balance = self.calculate_mass_balance()
                self.solve()

                # calculate solid number densities
                self.calculate_total_number()

                # calculate the percentage of each element condensed
                self.calculate_percentage_condensed()

                print("Finished equilibrating potential solids!")

            self.previous_number_densities = copy.copy(self.number_densities)
            self.previous_number_densities_solids = copy.copy(self.number_densities_solids)
            self.previous_number_densities_liquids = copy.copy(self.number_densities_liquids)
            self.previous_K = copy.copy(self.K)
            self.previous_removed_solids = copy.copy(self.removed_solids)
            self.previous_removed_liquids = copy.copy(self.removed_liquids)
            self.previous_temperature = self.temperature
            self.temperature -= self.dT
            if self.error_threshold > 1 * 10**-13:
                print("[!] Substantial solution error: {}".format(self.error_threshold))
            if self.error_threshold < 0.05:
                self.any_in = True
                self.any_out = True
            else:
                print("[!] Substantial solution error. Aborting... ({})".format(self.error_threshold))
                print("Errors: {}".format(self.errors))
                print("Percent condensed: {}".format(self.percent_element_condensed))
                sys.exit()
            if self.IS_SOLID:
                if len(self.condensing_solids) > 0:
                    print("Stable solids: {}".format(self.condensing_solids))
            if self.IS_LIQUID:
                if len(self.condensing_liquids) > 0:
                    print("Stable liquids: {}".format(self.condensing_liquids))