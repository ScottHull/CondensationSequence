import os
import pandas as pd


class OutputHandler:

    def __init__(self, abundances, phase_fname="phases.csv",
                 condensed_element_fname="condensed_element_percentages.csv"):
        self.phase_fname = phase_fname
        self.condensed_element_fname = condensed_element_fname
        if self.phase_fname in os.listdir(os.getcwd()):
            os.remove(self.phase_fname)
        if self.condensed_element_fname in os.listdir(os.getcwd()):
            os.remove((self.condensed_element_fname))
        self.phase_file = pd.DataFrame({})
        self.condensed_element_file = pd.DataFrame({})
        self.__initialize(abundances=abundances)

    def __initialize(self, abundances):
        elements = ["temperature"] + list(abundances.keys())
        for element in elements:
            self.condensed_element_file[element] = []
        self.phase_file["temperature"] = []

    def write_condensed_element_percent(self, temperature, percentages):
        row = {'temperature': temperature}
        elements = percentages.keys()
        for element in elements:
            row.update({element: percentages[element][-1]['percent']})
        self.condensed_element_file = self.condensed_element_file.append(row, ignore_index=True)

    def write_condensed_phases(self, temperature, condensing_solids, condensing_liquids):
        condensing_phases = condensing_solids + condensing_liquids
        df_headers = self.phase_file.keys()
        number_of_rows = self.phase_file.shape[0]  # 0 gives number of rows, 1 gives number of columns
        for phase in condensing_phases:
            if phase not in df_headers:
                self.phase_file[phase] = [0 for i in range(0, number_of_rows)]
        row = {'temperature': temperature}
        for phase in list(df_headers)[1:]:
            if phase in condensing_phases:
                row.update({phase: 1})
            else:
                row.update({phase: 0})
        self.phase_file = self.phase_file.append(row, ignore_index=True)
        print(self.phase_file)

    def write_to_files(self):
        self.phase_file.to_csv(self.phase_fname, index=False)
        self.condensed_element_file.to_csv(self.condensed_element_fname, index=False)
        print("[*] Outputs written to {}".format(os.getcwd()))
