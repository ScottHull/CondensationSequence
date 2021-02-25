import pandas as pd
import matplotlib.pyplot as plt

def plot_percent_condensed(path):
    data = pd.read_csv(path)
    headers = [i for i in data.keys() if i != "temperature"]
    temperatures = data["temperature"]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    for element in headers:
        if sum(data[element]) != 0.0:
            ax.plot(temperatures, data[element], linewidth=2.0, label=element)
    ax.set_xlim(list(temperatures)[-1], 1850)
    ax.invert_xaxis()
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("% Condensed")
    ax.set_title("% Element Condensed")
    ax.grid()
    ax.legend()
    plt.savefig("percent_element_condensed.png", format='png')

def plot_condensed_phases(path):
    df = pd.read_csv(path)
    phases = list(df.keys())[1:]
    temperatures = list(df['temperature'])
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    counter = 1
    for phase in phases:
        extant = [temperatures[index] for index, i in enumerate(list(df[phase])) if i == 1]
        ax.plot(extant, [counter for i in extant], linewidth=4.0)
        ax.annotate(phase, (extant[0] - 5, counter + 0.2))
        counter += 1
    ax.set_xlabel("Temperature")
    ax.set_title("Phase Temperature Stability Range")
    ax.grid()
    plt.show()

