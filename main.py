import collect_data
import k

class Condensation:

    def __init__(self, start_temperature, end_temperature, abundances, total_pressure, dt=2):
        self.R = 8.314
        self.start_temperature = start_temperature
        self.temperature = start_temperature
        self.end_temperature = end_temperature
        self.dt = dt
        self.total_pressure = total_pressure
        self.abundances = abundances
        self.normalized_abundances = self.normalize_abundances(abundances=self.abundances)
        self.gas_molecules_library = collect_data.get_atoms_from_molecule(path="data/Gasses.dat")
        self.solids_molecules_library = collect_data.get_atoms_from_molecule(path="data/Solids_Cp.dat", skiprows=1)
        self.partial_pressures = self.partial_pressure()


    def normalize_abundances(self, abundances):
        norm = {}
        total = sum(abundances.values())
        for element in abundances.keys():
            norm.update({element: self.abundances[element] / total})
        return norm

    def partial_pressure(self):
        # P_i = (n_i / n_total) * (P_total / (RT))
        partial_pressures = {}
        for element in self.normalized_abundances.keys():
            p_i = self.normalized_abundances[element] * (self.total_pressure / (self.R * self.temperature))
            partial_pressures.update({element: p_i})
        return partial_pressures

    def sequence(self):
        # ratio of element to major gasses in nebular cloud as the normalized abundances
        self.normalized_abundances = self.normalize_abundances(abundances=self.normalized_abundances)
        self.partial_pressures = self.partial_pressure()



c = Condensation(
    start_temperature=2000,
    end_temperature=200,
    abundances={"Si": 0.33, "Fe": 0.33, "Mg": 0.33},
    total_pressure=1.0
)
