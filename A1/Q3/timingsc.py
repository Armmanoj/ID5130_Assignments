import matplotlib.pyplot as plt
import numpy as np
import sys

def read_file(filename):
    """Reads the file and extracts the data."""
    datas = {"delta": [], "time": []}
    datae = {"delta": [], "time": []}
    datar = {"delta": [], "time": []}
    with open(filename, 'r') as file:
        i = 0
        for line in file:
            parts = line.split()
            try:
                Delta = float(parts[0])
                time = float(parts[2])
            except:
                continue
            if i<3:
                datas["delta"].append(np.log10(Delta))
                datas["time"].append(np.log10(time))
            elif i<6:
                datae["delta"].append(np.log10(Delta))
                datae["time"].append(np.log10(time))
            else:
                datar["delta"].append(np.log10(Delta))
                datar["time"].append(np.log10(time))
            i += 1
    return datas,datae,datar

def plot_data(datas, datae,datad):
    """Plots the threads vs. time."""
    plt.plot(datas["delta"], datas["time"], marker='o', linestyle='-', label = 'serial code')
    plt.plot(datae["delta"], datae["time"], marker='o', linestyle='-', label = "red-black aprpoach")
    plt.plot(datad["delta"], datad["time"], marker='o', linestyle='-', label = "diagonal approach")
    plt.title('Delta vs Tim ')
    plt.xlabel('log(Delta)')
    plt.ylabel('log(Time)')
    plt.grid(True)
    plt.legend()
    plt.savefig(sys.argv[1])
    plt.close()

datas,datae,datar = read_file("prof.log")
plot_data(datas,datae,datar)