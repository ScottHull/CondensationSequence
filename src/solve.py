import numpy as np
from scipy.optimize import root
from math import sqrt

from src import mass_balance


class Solve:
    """
    Should the number if iterations exceed the max number and the error still be too large, we need to kick the solver
    in a different direction.
    - We can track the error to see if the guess was incorrect.
    - We can adjust the initial guess based on this direction.
    """

    def __init__(self, condensing_solids, number_densities, element_gas_appearances, gas_molecules_library, K,
                 number_densities_from_partial_pressure, temperature, solid_molecules_library):
        self.number_densities = number_densities
        self.element_gas_appearances = element_gas_appearances
        self.gas_molecules_library = gas_molecules_library
        self.K = K
        self.number_densities_from_partial_pressure = number_densities_from_partial_pressure
        self.temperature = temperature
        self.solid_molecules_library = solid_molecules_library
        self.liquid_molecules_library = {}
        self.condensing_solids = condensing_solids
        self.condensing_liquids = []
        self.names = []

        self.initial_guess = {}
        self.current_guess = {}
        self.initial = True

        self.error = None
        self.error_count = 1
        self.error_threshold = 1 * 10 ** 10  # placeholder to keep while-loop running
        self.ERROR_THRESHOLD_SUCCESS = 1 * 10 ** -13
        self.MAX_ERROR_COUNT = 50
        self.error_history = []

    def __setup(self):
        names = list(self.number_densities.keys())
        guess_number_density = np.array(
            [self.number_densities[p] for p in names])  # a list of abundances corresponding to the elements list
        # a list of abundances corresponding to the elements list
        for s in self.condensing_solids:
            if s not in names:
                names.append(s)
                guess_number_density = np.append(guess_number_density,
                                                 1 * 10 ** -13)  # append a guess for all existing solids
        return names, guess_number_density

    def __solve(self, names, guess):
        args = (names, self.element_gas_appearances, self.gas_molecules_library, self.K,
                self.number_densities_from_partial_pressure,
                self.temperature, self.condensing_solids, self.solid_molecules_library, self.condensing_liquids,
                self.liquid_molecules_library)
        # finds the root of the mass balance equation where the root is the number density (initial guess is just the partial pressures)
        number_densities = root(mass_balance.mass_balance, guess, args=args, method='lm',
                                options={'maxiter': 100000000,
                                         'ftol': 1.e-15})  # calculate the activities of gas phases corresponding to the element dict
        self.number_densities = dict(zip(names,
                                         number_densities.x))  # make a dictionary where the element is the key and the activity is the value
        errors = mass_balance.mass_balance(number_densities.x, *args)
        error_threshold = sqrt(
            sum([i ** 2 for index, i in enumerate(errors) if names[index] not in self.condensing_solids]))
        return number_densities.x, error_threshold

    def __kick_guess(self):
        self.initial_guess[-1] *= 1 * 10 ** -1

    def solver(self):
        if self.initial:
            print(">>> INITIAL SOLVER SETUP")
            self.names, self.initial_guess = self.__setup()
            self.current_guess = self.initial_guess
            self.current_guess, self.error_threshold = self.__solve(names=self.names, guess=self.initial_guess)
            self.error_history.append(self.error_threshold)
            self.initial = False
        while self.error_count <= self.MAX_ERROR_COUNT and self.error_threshold > self.ERROR_THRESHOLD_SUCCESS:
            print(">>> IN MAIN ERROR LOOP ({}, count {})".format(self.error_threshold, self.error_count))
            number_densities, error_threshold = self.__solve(names=self.names, guess=self.initial_guess)
            self.error_history.append(error_threshold)
            # self.current_guess = number_densities
            self.error_threshold = error_threshold
            self.temperature -= 0.0001
            self.error_count += 1
        if self.error_count > self.MAX_ERROR_COUNT and self.error_threshold > self.ERROR_THRESHOLD_SUCCESS:
            print(">>> [!] Restarting solver... (trouble solid: {})...".format(self.condensing_solids[-1]))
            self.__kick_guess()
            self.error_count = 0
            self.current_guess = self.initial_guess
            self.solver()
        elif self.error_threshold <= self.ERROR_THRESHOLD_SUCCESS:
            print(">>> [*] SOLVER SUCCESS ({})".format(self.error_threshold))
            return self.number_densities, self.error_threshold

# def solve(self):
#     self.mass_balance = self.calculate_mass_balance()
#     number_densities_from_partial_pressure = self.calculate_mass_balance()
#     # list of elements given in the input
#     names = list(self.number_densities.keys())
#     guess_number_density = np.array(
#         [self.number_densities[p] for p in names])  # a list of abundances corresponding to the elements list
#     if not self.initial:
#         # a list of abundances corresponding to the elements list
#         for s in self.condensing_solids:
#             if s not in names:
#                 names.append(s)
#                 guess_number_density = np.append(guess_number_density,
#                                                  force_guess.get_guess(
#                                                      molecule=s,
#                                                      mass_balance=self.mass_balance,
#                                                      total_n=self.total_number,
#                                                      elements_in_solid=self.elements_in_solid,
#                                                      temperature=self.temperature,
#                                                      number_densities=self.number_densities,
#                                                      percent_condensed=self.percent_element_condensed
#                                                  ))  # append a guess for all existing solids
#         for s in self.condensing_liquids:
#             if s not in names:
#                 names.append(s)
#                 guess_number_density = np.append(guess_number_density,
#                                                  1.e-13)  # append a guess for all existing solids
#     self.initial = False
#     args = (names, self.element_gas_appearances, self.gas_molecules_library, self.K,
#             number_densities_from_partial_pressure,
#             self.temperature, self.condensing_solids, self.solid_molecules_library, self.condensing_liquids,
#             self.liquid_molecules_library)
#     # finds the root of the mass balance equation where the root is the number density (initial guess is just the partial pressures)
#     print("Solving system...")
#     number_densities = root(mass_balance.mass_balance, guess_number_density, args=args, method='lm',
#                             options={'maxiter': 100000000,
#                                      'ftol': 1.e-15})  # calculate the activities of gas phases corresponding to the element dict
#     self.number_densities = dict(zip(names,
#                                      number_densities.x))  # make a dictionary where the element is the key and the activity is the value
#     self.errors = mass_balance.mass_balance(number_densities.x, *args)
#     self.error_threshold = sqrt(
#         sum([i ** 2 for index, i in enumerate(self.errors) if names[index] not in self.condensing_solids]))
#     self.errors = dict(zip([n for n in names if n not in self.condensing_solids and n
#                             not in self.condensing_liquids], self.errors))
#
#     print("Solved system!")
