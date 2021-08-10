import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt

test_temperature_range = np.arange(100, 4000, 10)

base_liq_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 1", index_col="Reactant")
gas_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 2", index_col="Product")
liq_df = pd.read_excel("data/MAGMA_Thermodynamic_Data.xlsx", sheet_name="Table 3", index_col="Product")

base_liq_K = {}
for species in base_liq_df.index:
    A, B = base_liq_df["A"][species], base_liq_df["B"][species]
    vals = []
    for T in test_temperature_range:
        vals.append(log10((10 ** (A + (B / T))) ** -1))  # log10(K) = A + B/T
    base_liq_K.update({species: vals})

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for species in base_liq_K.keys():
    ax.plot(
        test_temperature_range,
        base_liq_K[species],
        linewidth=2.0,
        label=species
    )
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("log10(K)")
ax.grid()
ax.legend()

plt.show()
