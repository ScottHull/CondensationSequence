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
        self.R = 8.314
        self.start_temperature = start_temperature
        self.temperature = start_temperature
        self.end_temperature = end_temperature
        self.dT = dT
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
        self.partial_pressures = self.mass_balance()
        self.previous_partial_pressures = {}
        self.condensing_solids = []
        self.removed_solids = []
        self.previous_removed_solids = []
        self.K = {}
        self.previous_K = {}
        self.number_densities = {}
        self.previous_number_densities = {}
        self.number_densities_solids = {}
        self.elements_in_solid = []
        self.total_elements_condensed = {}

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

    def mass_balance(self):
        """
        This is equation 5 from Unterborn and Panero 2017, where N_i is proportional to the partial pressure of element i.
        N_i = (a_i / a_total) * (P_total / RT)
        :return:
        """
        mass_balance_from_pp = {}
        for element in self.normalized_abundances.keys():  # what is the 10**-2 for?
            N_i = self.normalized_abundances[element] * (self.total_pressure / ((self.R * 10 ** -2) * self.temperature))
            mass_balance_from_pp.update({element: N_i})
        return mass_balance_from_pp


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

        # normalize the given abundances to 100 mole% (sum = 1)
        # self.normalized_abundances = self.normalize_abundances(abundances=self.normalized_abundances)
        # calculate the partial pressure of each element in the system
        # self.partial_pressures = self.mass_balance()
        # return a dictionary of K values
        self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
                         temperature=self.temperature, gas_methods=self.gas_methods)
        number_densities_from_partial_pressure = self.mass_balance()
        # list of elements given in the input
        elements = self.normalized_abundances.keys()

        # a list of abundances corresponding to the elements list
        guess_number_density = np.array([number_densities_from_partial_pressure[element] for element in elements])
        args = (elements, self.element_gas_appearances, self.gas_molecules_library, self.K, number_densities_from_partial_pressure,
                self.temperature, self.condensing_solids, self.solid_molecules_library)
        # finds the root of the mass balance equation where the root is the number density (initial guess is just the partial pressures)
        number_densities = root(mass_balance.mass_balance, guess_number_density, args=args, method='lm',
                                options={'maxiter': 100000000,
                                         'ftol': 1.e-15})  # calculate the activities of gas phases corresponding to the element dict
        self.number_densities = dict(zip(elements,
                                         number_densities.x))  # make a dictionary where the element is the key and the activity is the value

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

        errors = mass_balance.mass_balance(number_densities.x, *args)
        error_threshold = sqrt(sum([i ** 2 for i in errors]))

        # check in stable solid molecules into the system
        any_in = True
        any_out = True
        while any_in and any_out:  # enter a while loop
            # in_solid is the solid that condenses
            # in_temperature is the temperature at which the solid condenses
            in_solid, in_temp = solids.check_in(
                solids=self.solid_molecules_library,
                number_densities=self.number_densities,
                temperature=self.temperature,
                K_dict=self.K,
                condensing_solids=self.condensing_solids,
                temperature_old=self.temperature + self.dT,
                K_dict_old=self.previous_K,
                number_densities_old=self.previous_number_densities,
                removed_solids=self.removed_solids,
                removed_solids_old=self.previous_removed_solids
            )

            if error_threshold > 1 * 10 ** -13:  # if the error is smaller than the threshold
                in_solid = False
                in_temp = 0
            else:  # if the error is larger than the threshold
                any_in = False
                in_temp = 0

            new_solid = [in_solid]

            out_solid, out_temp = solids.check_out(condensing_solids=self.condensing_solids,
                                                   number_density_solids=self.number_densities,
                                                   number_density_solids_old=self.previous_number_densities,
                                                   temperature=self.temperature,
                                                   temperature_old=self.temperature + self.dT)

            if out_solid == False:
                any_out = False
                out_temp = 0

            if in_temp > out_temp:  # if the appearance temperature is greater than the disappearance temperature
                self.condensing_solids.append(in_solid)
                any_out = False  # a solid has not dropped out
                any_in = True  # there is a new solid in
                # track element in solid phase if it not already is being done
                for s in self.solid_molecules_library[in_solid]:
                    if s not in self.elements_in_solid:
                        self.elements_in_solid.append(s)
            else:  # if the disappearance temperature is greater than the appearance temperature
                any_out = True  # there are exiting solids
                any_in = False  # there are no new solids
                self.temperature = out_temp  # set the temperature to the solid-out temperature
                self.condensing_solids.remove(out_solid)  # remove the outgoing solid from the condensing solids list
                del self.number_densities_solids[out_solid]  # remove the outgoing solid from the solids number density dict

                self.removed_solids.append(out_solid)  # track the exit of the solid
                self.normalized_abundances = self.normalize_abundances(self.abundances)
                self.partial_pressures = self.mass_balance()

            # if any_out is True or any_in is True:  # must reset parameters if a condensation change is ocurring
            #     counter_out = 0
            #     self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solid_molecules_library,
            #              temperature=self.temperature, gas_methods=self.gas_methods)
            #     self.partial_pressures = self.mass_balance()
            #     args = (Element_dict, K_dict, Par_Pres_dict, gasses, Name, T, condensing_solids)
            #
            #     gas_activity = root(fun.Mass_balance_fun, guess, args=args, method='lm',options={'maxiter': 100000000,'ftol':1.e-15})
            #     guess = gas_activity.x
            #     guess_dict = dict(zip(Name, guess))
            #
            #     errors = []
            #     error_dict = {}
            #     lst_sqr = 0.
            #     errors = fun.Mass_balance_fun(gas_activity.x, *args)
            #     error_dict = dict(zip(Name, errors))
            #
            #     for i,j in dict.iteritems(error_dict):
            #         if i not in condensing_solids:
            #             lst_sqr += pow(j, 2.)
            #
            #     n_x = dict(zip(Name, gas_activity.x))
            #     Sum_Par_pres = get_data.get_total_atoms(Element_dict, n_x, gasses, condensing_solids, K_dict, RT,
            #                                             Element_solid_dict, solids)
            #
            #     guess = gas_activity.x
            #     guess_dict = dict(zip(Name, guess))










a = {'Ni': 1.66e+16, 'C': 2.692e+18, 'F': 363100000000000.0, 'H': 1e+22, 'K': 1072000000000000.0,
     'Mn': 2692000000000000.0, 'Mg': 3.981e+17, 'O': 4.898e+18, 'Ne': 8.511e+17, 'P': 2570000000000000.0,
     'S': 1.318e+17, 'Ti': 891300000000000.0, 'N': 6.761e+17, 'Co': 977200000000000.0, 'Cl': 3162000000000000.0,
     'Ca': 2.188e+16, 'Si': 3.236e+17, 'Al': 2.818e+16, 'Ar': 2.512e+16, 'Fe': 3.162e+17, 'Na': 1.738e+16,
     'Cr': 4365000000000000.0, 'He': 8.511e+20}
c = Condensation(
    start_temperature=2500,
    end_temperature=200,
    abundances=a,
    total_pressure=1 * 10 ** -3
)
c.sequence()
