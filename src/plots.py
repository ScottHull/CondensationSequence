import matplotlib.pyplot as plt
import numpy as np


def plot_percent_condensed(percent_condensed_dict, legend=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    temperatures = [percent_condensed_dict[i]['temperature'] for i in percent_condensed_dict.keys()]
    percentages = [percent_condensed_dict[i]['percent'] for i in percent_condensed_dict.keys()]
    for element in percent_condensed_dict.keys():
        pct = percent_condensed_dict[element]
        ax.plot(temperatures, percentages, label=element)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("% Condensed")
    ax.set_title("% Element Condensed")
    ax.grid()
    if legend:
        ax.legend()
    plt.show()
