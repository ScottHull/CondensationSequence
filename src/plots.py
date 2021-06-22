import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'size': 14}

matplotlib.rc('font', **font)

def __read_file(path):
    data = {}
    with open(path, 'r') as infile:
        temperature = float(infile.readline().replace("\n", ""))
        lines = infile.read().splitlines()
        for i in lines:
            element, d = i.split(",")
            data.update({element: float(d)})
        infile.close()
        return temperature, data

def __get_mole_fractions(temperatures, data, solids=True):
    fractions = {}
    for t in temperatures:
        abundances_at_time = {}
        for element in data.keys():
            if solids:
                if "_s" in element:
                    if t in data[element].keys():
                        abundances_at_time.update({element: data[element][t]})
            else:
                if t in data[element].keys():
                    abundances_at_time.update({element: data[element][t]})
        for element in abundances_at_time.keys():
            if element not in fractions.keys():
                fractions.update({element: {}})
            fractions[element].update({t: abundances_at_time[element] / sum(abundances_at_time.values()) * 100.0})
    return fractions

def plot_mole_fractions(base_path):
    path = base_path + "/number_densities"
    data = {}
    temperatures = []
    for f in os.listdir(path):
        temperature, d = __read_file(path=path + "/{}".format(f))
        temperatures.append(temperature)
        for element in d.keys():
            if element not in data.keys():
                data.update({element: {}})
            data[element].update({temperature: d[element]})
    fractions = __get_mole_fractions(temperatures=temperatures, data=data)
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    for element in fractions.keys():
        temperatures = fractions[element].keys()
        d = [fractions[element][t] for t in temperatures]
        ax.plot(
            temperatures,
            d,
            linewidth=2.0,
            label=element.replace("_s", "")
        )
    ax.invert_xaxis()
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Mole Fraction")
    ax.grid()
    ax.legend()

    plt.show()


def plot_number_densities(base_path):
    path = base_path + "/number_densities"
    data = {}
    temperatures = []
    for f in os.listdir(path):
        temperature, d = __read_file(path=path + "/{}".format(f))
        temperatures.append(temperature)
        for element in d.keys():
            if element not in data.keys():
                data.update({element: {}})
            data[element].update({temperature: d[element]})
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    for element in data.keys():
        if "_s" in element:
            temperatures = data[element].keys()
            d = [data[element][t] for t in temperatures]
            ax.plot(
                temperatures,
                d,
                linewidth=2.0,
                label=element.replace("_s", "")
            )
    ax.invert_xaxis()
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Number Density")
    ax.grid()
    ax.legend()

    plt.show()
