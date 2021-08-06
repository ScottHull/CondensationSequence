import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fegley_1986 = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/CondensationSequence/data/liquid_pseudospecies_fegley_1986.dat"
janaf_oxides_file = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/CondensationSequence/data/Liquids.dat"
janaf_path = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/CondensationSequence/data/Liquids/"

fegley_df = pd.read_csv(fegley_1986, delimiter="\t", header=None, skiprows=1, index_col=0)
janaf_oxides_df = pd.read_csv(janaf_oxides_file, delimiter="\t", header=None)

do_not_plot = ['K2O1']
data = {}
for row in janaf_oxides_df.index:
    oxide = janaf_oxides_df[0][row]
    if oxide not in do_not_plot:
        A = fegley_df[2][oxide]
        B = fegley_df[3][oxide]
        f = janaf_path + oxide + ".dat"
        oxide_df = pd.read_csv(f, delimiter="\t", skiprows=2, header=None)
        data.update({oxide: {
            'temperature': [],
            'janaf': [],
            'fegley': []
        }})
        for t_row in oxide_df.index:
            temperature = oxide_df[0][t_row]
            logK_janaf = oxide_df[7][t_row]
            K_janaf = 10 ** logK_janaf
            logK_fegley = A + (B / temperature)
            K_fegley = 10 ** logK_fegley
            data[oxide]['temperature'].append(temperature)
            data[oxide]['janaf'].append(logK_janaf)
            data[oxide]['fegley'].append(logK_fegley)
            print(oxide, logK_janaf, logK_fegley)

fig = plt.figure()
ax = fig.add_subplot(111)
for oxide in data.keys():
    diff = [x - y for x, y in zip(data[oxide]['janaf'], data[oxide]['fegley'])]
    ax.plot(data[oxide]['temperature'], diff, linewidth=2.0, label=oxide)
ax.set_xlabel("Temperature")
ax.set_ylabel("logK_J - logK_F")
ax.set_title("logK_JANAF - logK_Fegley")
ax.grid()
ax.legend()
plt.show()
