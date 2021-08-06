from src.main import Condensation
from analyze import plot

protosolar_disk = {'Ni': 1.66e+16, 'C': 2.692e+18, 'F': 363100000000000.0, 'H': 1e+22, 'K': 1072000000000000.0,
                   'Mn': 2692000000000000.0, 'Mg': 3.981e+17, 'O': 4.898e+18, 'Ne': 8.511e+17, 'P': 2570000000000000.0,
                   'S': 1.318e+17, 'Ti': 891300000000000.0, 'N': 6.761e+17, 'Co': 977200000000000.0,
                   'Cl': 3162000000000000.0,
                   'Ca': 2.188e+16, 'Si': 3.236e+17, 'Al': 2.818e+16, 'Ar': 2.512e+16, 'Fe': 3.162e+17, 'Na': 1.738e+16,
                   'Cr': 4365000000000000.0, 'He': 8.511e+20}

bulk_moon = {
    "Si": 16.02645359,
    "Mg": 18.80239392,
    "Al": 1.651637357,
    "Ti": 0.045956464,
    "Fe": 3.72636474,
    "Ca": 1.270525392,
    "Na": 0.034834503,
    "K": 0.001833649,
    # "Zn": 5.30606E-05,
    "O": 58.43994733
}  # atom%

liquid_species = [
    "SiO2", "MgO", "FeO", "Fe2O3", "Fe3O4", "CaO",
    "Al2O3", "TiO2", "Na2O", "K2O", "ZnO", "MgSiO3", "Mg2SiO4",
    "MgAl2O4", "MgTiO3", "MgTi2O5", "Mg2TiO4", "Al6Si2O13", "CaAl2O4",
    "CaAl4O7", "Ca12Al14O33", "CaSiO3", "CaAl2Si2O8", "CaMgSi2O6",
    "Ca2MgSi2O7", "Ca2Al2SiO7", "CaTiO3", "Ca2SiO4", "CaTiSiO5",
    "FeTiO3", "Fe2SiO4", "FeAl2O4", "CaAl12O19", "Mg2Al4Si5O18",
    "Na2SiO3", "Na2Si2O5", "NaAlSiO4", "NaAlSi3O8", "NaAlO2", "Na2TiO3",
    "NaAlSi2O6", "K2SiO3", "K2Si2O5", "KAlSiO4", "KAlSi3O8", "KAlO2",
    "KAlSi2O6", "K2Si4O9", "KCaAlSi2O7", "Zn2SiO4", "ZnTiO3", "Zn2TiO4",
    "ZnAl2O4"
]

gas_species = [
    "O", "O2", "Mg", "MgO", "Si", "SiO", "SiO2", "Fe", "FeO", "Al",
    "AlO", "AlO2", "Al2O", "Al2O2", "Ca", "CaO", "Na", "Na2", "NaO", "Na2O", "Na+",
    "K", "K2", "KO", "K2O", "K+", "Ti", "TiO", "TiO2", "Zn", "ZnO", "e−"
]

c = Condensation(
    start_temperature=4000,
    end_temperature=500,
    abundances=bulk_moon,
    total_pressure=1 * 10 ** -3,
    solid=False,
    liquid=True,
    gas=True
)

c.sequence()

# plot.plot_percent_condensed(path="condensed_element_percentages.csv")
# plot.plot_condensed_phases(path="phases.csv")
