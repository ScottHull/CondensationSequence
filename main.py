import collect_data
import k
import mass_balance
import numpy as np
from scipy.optimize import root


class Condensation:

    def __init__(self, start_temperature, end_temperature, abundances, total_pressure, dt=2):
        self.R = 8.314
        self.start_temperature = start_temperature
        self.temperature = start_temperature
        self.end_temperature = end_temperature
        self.dt = dt
        self.total_pressure = total_pressure
        self.gas_methods = collect_data.get_methods(
            path="data/Gasses.dat")  # methods for calculating K for gasses (i.e. JANAF)
        self.gas_molecules_library = collect_data.get_atoms_from_molecule(path="data/Gasses.dat")
        self.solids_molecules_library = collect_data.get_atoms_from_molecule(path="data/Solids_Cp.dat", skiprows=1)
        self.abundances = abundances
        self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
        self.partial_pressures = self.partial_pressure()
        self.condensing_solids = []
        self.K = {}

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

    def partial_pressure(self):
        # P_i = (n_i / n_total) * (P_total / (RT))
        partial_pressures = {}
        for element in self.normalized_abundances.keys():  # what is the 10**-2 for?
            p_i = self.normalized_abundances[element] * (self.total_pressure / ((self.R * 10**-2) * self.temperature))
            partial_pressures.update({element: p_i})
        return partial_pressures

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
        # self.partial_pressures = self.partial_pressure()
        # return a dictionary of K values
        self.K = k.get_K(gas_molecules=self.gas_molecules_library, solid_molecules=self.solids_molecules_library,
                         temperature=self.temperature, gas_methods=self.gas_methods)
        # list of elements given in the input
        elements = self.normalized_abundances.keys()
        # a list of abundances corresponding to the elements list
        abundances = np.array([self.partial_pressures[element] for element in elements])
        args = (self.normalized_abundances, self.K, self.partial_pressures, self.gas_molecules_library,
                self.normalized_abundances.keys(), self.temperature, self.condensing_solids)
        gas_activity = root(mass_balance.mass_balance, abundances, args=args, method='lm',
                            options={'maxiter': 100000000, 'ftol': 1.e-15})


a = {'Ni': 1.66e+16, 'C': 2.692e+18, 'F': 363100000000000.0, 'H': 1e+22, 'K': 1072000000000000.0,
     'Mn': 2692000000000000.0, 'Mg': 3.981e+17, 'O': 4.898e+18, 'Ne': 8.511e+17, 'P': 2570000000000000.0,
     'S': 1.318e+17, 'Ti': 891300000000000.0, 'N': 6.761e+17, 'Co': 977200000000000.0, 'Cl': 3162000000000000.0,
     'Ca': 2.188e+16, 'Si': 3.236e+17, 'Al': 2.818e+16, 'Ar': 2.512e+16, 'Fe': 3.162e+17, 'Na': 1.738e+16,
     'Cr': 4365000000000000.0, 'He': 8.511e+20}
c = Condensation(
    start_temperature=2500,
    end_temperature=200,
    abundances=a,
    total_pressure=1 * 10**-3
)
c.sequence()
