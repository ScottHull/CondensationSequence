import matplotlib.pyplot as plt
import numpy as np


def plot_percent_condensed(start_temperature, end_temperature, dt, percent_condensed_dict):
    t_range = reversed(list(np.arange(end_temperature, start_temperature + dt, dt)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for element in percent_condensed_dict.keys():
        pct = percent_condensed_dict[element]
        ax.plot(t_range, )

