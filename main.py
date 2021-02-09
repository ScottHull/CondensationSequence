import collect_data
import k
import mass_balance
import solids
import total_atoms
import numpy as np
from scipy.optimize import root
from math import sqrt
import sys


class Condensation:

    def __init__(self, start_temperature, end_temperature, abundances, total_pressure, dT=2):
        self.R = 8.3144621e-2
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
        self.solid_molecules_library = collect_data.get_atoms_from_molecule(path="data/Solids_Cp.dat", skiprows=1,
                                                                            solid=True)
        self.abundances = abundances
        self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
        self.element_gas_appearances = collect_data.element_appearances_in_molecules(abundances=self.abundances,
                                                                                     library=self.gas_molecules_library)
        self.element_solid_appearances = collect_data.element_appearances_in_molecules(abundances=self.abundances,
                                                                                       library=self.solid_molecules_library)

        self.mass_balance = self.calculate_mass_balance()
        self.condensing_solids = []
        self.removed_solids = []
        self.previous_removed_solids = []
        self.K = {}  # a dictionary of all K equilibrium constant values
        self.previous_K = {}  # tracks K values from the previous temperature iteration
        self.number_densities = {}  # element number densities
        self.previous_number_densities = {}  # tracks number density solution from the previous temperature iteration
        self.number_densities_solids = {}  # number densities of all solids
        self.elements_in_solid = []  # all elements currently occupying a condensing solid
        self.total_elements_condensed = {}  # tracks the number density of all elements in the condensed phase
        self.percent_element_condensed = {}
        for element in self.abundances.keys():  # initial setup
            self.percent_element_condensed.update({element: 0})
            self.total_elements_condensed.update({element: 0})
        
        self.any_in = True
        self.any_out = True
        self.errors = []

    def solve(self):
        number_densities_from_partial_pressure = self.calculate_mass_balance()
        # list of elements given in the input
        names = list(self.normalized_abundances.keys())
        guess_number_density = np.array([self.number_densities[p] for p in names])  # a list of abundances corresponding to the elements list
        if not self.initial:
            # a list of abundances corresponding to the elements list
            for s in self.condensing_solids:
                names.append(s)
                guess_number_density = np.append(guess_number_density, 1.e-13)  # append a guess for all existing solids
        self.initial = False
        args = (names, self.element_gas_appearances, self.gas_molecules_library, self.K,
                number_densities_from_partial_pressure,
                self.temperature, self.condensing_solids, self.solid_molecules_library)
        # finds the root of the mass balance equation where the root is the number density (initial guess is just the partial pressures)
        print("Solving system...")
        number_densities = root(mass_balance.mass_balance, guess_number_density, args=args, method='lm',
                                options={'maxiter': 100000000,
                                         'ftol': 1.e-15})  # calculate the activities of gas phases corresponding to the element dict
        self.number_densities = dict(zip(names,
                                         number_densities.x))  # make a dictionary where the element is the key and the activity is the value
        self.errors = mass_balance.mass_balance(number_densities.x, *args)
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
        K = k.get_K_gas(molecules=self.gas_molecules_library.keys(), methods=self.gas_methods,
                        temperature=self.temperature)
        X = K["H1"] ** 2
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


    def equilibrate_solids(self, error_threshold):
        # in_solid is the solid that condenses
        # in_temperature is the temperature at which the solid condenses
        in_solid, in_temp = solids.check_in(
            solids=self.solid_molecules_library,
            number_densities=self.number_densities,
            temperature=self.temperature,
            K_dict=self.K,
            condensing_solids=self.condensing_solids,
            temperature_old=self.previous_temperature,
            K_dict_old=self.previous_K,
            number_densities_old=self.previous_number_densities,
            removed_solids=self.removed_solids,
            removed_solids_old=self.previous_removed_solids
        )
        if error_threshold > 1 * 10 ** -13:  # if the error is smaller than the threshold
            in_solid = False
            in_temp = 0
        if in_solid == False:  # ???
            self.any_in = False  # break out of check-in loop
            in_temp = 0

        out_solid, out_temp = solids.check_out(condensing_solids=self.condensing_solids,
                                               number_density_solids=self.number_densities,
                                               number_density_solids_old=self.previous_number_densities,
                                               temperature=self.temperature,
                                               temperature_old=self.previous_temperature)

        return in_solid, in_temp, out_solid, out_temp


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
                         temperature=self.temperature, gas_methods=self.gas_methods)
        self.mass_balance = self.calculate_mass_balance()


        # a list of abundances corresponding to the elements list
        self.number_densities = self.mass_balance

        while self.temperature > self.end_temperature:
            print("Entering solution loop...")
            print("AT TEMPERATURE: {}".format(self.temperature))
            self.any_in = True
            self.any_out = True
            self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
            self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
                             temperature=self.temperature, gas_methods=self.gas_methods)

            self.solve()  # run the root solver

            # the sum of all number densities across gasses and solids
            total_N = total_atoms.calculate_total_N(
                gas_element_appearances_in_molecules=self.element_gas_appearances,  # all active gas molecules in system
                solid_element_appearances_in_molecules=self.element_solid_appearances,  # all active solid molecules in system
                element_number_densities=self.number_densities,  # computed elemental number densities from root finder
                condensing_solids=self.condensing_solids,  # list of all condensing solids
                gas_molecule_library=self.gas_molecules_library,  # stoich library for all supported gas molecules
                solid_molecule_library=self.solid_molecules_library,  # stoich library for all supported solid molecules
                K_dict=self.K,  # all equilibrium constants
                temperature=self.temperature,
            )

            for i in self.condensing_solids:
                # the solid number density is equal to the fractional number density of the molecule?
                self.number_densities_solids.update({i: self.number_densities[i] / total_N})

            # make sure that the guess number densities return 0's, or very very close to 0's.
            error_threshold = sqrt(sum([i ** 2 for i in self.errors]))

            # check in stable solid molecules into the system
            while self.any_in and self.any_out:  # enter a while loop
                print("Equilibrating potential solids...")
                """
                Within this while loop, the following occurs
                1. Check in solids
                2. Check out solids
                3. If check in temp > check out temp: (check in solid)
                    a. Get total N (final)
                    b. Append solid to condensing solids list
                    c. Set self.any_out to False and self.any_in to True (there are new solids)
                    d. Append solid elements to "elements_in_solids" list if it doesn't already exist
                    e. Adjust guesses based on solid in order to keep the solver moving
                    f. Append solid to Name and guess
                    g. Reset parameters so that the next loop has the updated information
                    h. Recalculate abundance and guess number density
                4. If check out temp > check in temp: (check out solid)
                    a. Set self.any_out to True and self.any_in to False
                    b. Remove solid from condensing_solids list, num_solids dict, guess_dict, append to removed solids list
                    c. Recalculate abundance and guess number density
                    d. Reset parameters so that the next loop has the updated information
                5. If self.any_in is True or self.any_out is True (solids have entered or left the system)
                    a. Recalculate K
                    b. Recalculate n_i with root solver
                6. Recalculate n for solids
                7. Reset parameters
                8. The total elements condensed (tracked in the total_elements_condensed dict) is the previous plus (stoich * n_i) [using elements_in_solid dict] <- looks to just be a sum of the total number n in the condensed phase
                9. The percent elements condensed is the total_elements_condensed (number density) value divided by the total number density
                """

                in_solid, in_temp, out_solid, out_temp = self.equilibrate_solids(error_threshold=error_threshold)
                if in_solid == "Al2O3_s":
                    print(in_solid, in_temp, out_solid, out_temp)

                if out_solid == False:
                    self.any_out = False
                    out_temp = 0

                if in_temp > out_temp:  # if the appearance temperature is greater than the disappearance temperature
                    print("IN SOLID: {} ({} K)".format(in_solid, in_temp))
                    self.condensing_solids.append(in_solid)
                    self.any_out = False  # a solid has not dropped out
                    self.any_in = True  # there is a new solid in
                    # track element in solid phase if it not already is being done
                    for s in self.solid_molecules_library[in_solid]:
                        if s not in self.elements_in_solid:
                            self.elements_in_solid.append(s)
                elif out_temp > in_temp:  # if the disappearance temperature is greater than the appearance temperature
                    print("OUT SOLID: {} ({} K)".format(out_solid, out_temp))
                    self.any_out = True  # there are exiting solids
                    self.any_in = False  # there are no new solids
                    self.temperature = out_temp  # set the temperature to the solid-out temperature
                    self.condensing_solids.remove(out_solid)  # remove the outgoing solid from the condensing solids list
                    del self.number_densities_solids[out_solid]  # remove the outgoing solid from the solids number density dict

                    self.removed_solids.append(out_solid)  # track the exit of the solid

            print("Finished equilibrating potential solids!")
            # get the total number density of the condensed elements
            for s in self.condensing_solids:
                solid_stoich = self.solid_molecules_library[s]
                for element in solid_stoich:
                    stoich = solid_stoich[element]
                    element_number_density = self.number_densities[element]
                    self.total_elements_condensed[element] += stoich * element_number_density

            for element in self.total_elements_condensed:
                element_number_density = self.total_elements_condensed[element]
                self.percent_element_condensed[element] = element_number_density / self.mass_balance[element]

            self.previous_number_densities = self.number_densities
            self.previous_K = self.K
            self.previous_removed_solids = self.removed_solids
            self.previous_temperature = self.temperature
            self.temperature -= self.dT












a = {'Ni': 1.66e+16, 'C': 2.692e+18, 'F': 363100000000000.0, 'H': 1e+22, 'K': 1072000000000000.0,
     'Mn': 2692000000000000.0, 'Mg': 3.981e+17, 'O': 4.898e+18, 'Ne': 8.511e+17, 'P': 2570000000000000.0,
     'S': 1.318e+17, 'Ti': 891300000000000.0, 'N': 6.761e+17, 'Co': 977200000000000.0, 'Cl': 3162000000000000.0,
     'Ca': 2.188e+16, 'Si': 3.236e+17, 'Al': 2.818e+16, 'Ar': 2.512e+16, 'Fe': 3.162e+17, 'Na': 1.738e+16,
     'Cr': 4365000000000000.0, 'He': 8.511e+20}
c = Condensation(
    start_temperature=1732,
    end_temperature=200,
    abundances=a,
    total_pressure=1 * 10 ** -3
)
c.sequence()
