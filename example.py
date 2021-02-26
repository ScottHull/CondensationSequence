from src.main import Condensation
from analyze import plot

a = {'Ni': 1.66e+16, 'C': 2.692e+18, 'F': 363100000000000.0, 'H': 1e+22, 'K': 1072000000000000.0,
     'Mn': 2692000000000000.0, 'Mg': 3.981e+17, 'O': 4.898e+18, 'Ne': 8.511e+17, 'P': 2570000000000000.0,
     'S': 1.318e+17, 'Ti': 891300000000000.0, 'N': 6.761e+17, 'Co': 977200000000000.0, 'Cl': 3162000000000000.0,
     'Ca': 2.188e+16, 'Si': 3.236e+17, 'Al': 2.818e+16, 'Ar': 2.512e+16, 'Fe': 3.162e+17, 'Na': 1.738e+16,
     'Cr': 4365000000000000.0, 'He': 8.511e+20}

c = Condensation(
    start_temperature=2500,
    end_temperature=500,
    abundances=a,
    total_pressure=1 * 10 ** -3,
    solid=True,
    liquid=False,
    gas=True
)

c.sequence()

# plot.plot_percent_condensed(path="condensed_element_percentages.csv")
# plot.plot_condensed_phases(path="phases.csv")
