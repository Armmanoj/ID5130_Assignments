import matplotlib.pyplot as plt
import sys

def read_file(filename):
    """Reads the file and extracts the data."""
    datae = {"threads": [], "time": []}
    datar = {"threads": [], "time": []}
    with open(filename, 'r') as file:
        i = 0
        for line in file:
            parts = line.split()
            try:
                threads = int(parts[1])
                time = float(parts[2])
            except:
                continue
            if i<4:
                datae["threads"].append(threads)
                datae["time"].append(time)
            else:
                if threads!=16:
                    datar["threads"].append(threads)
                    datar["time"].append(time)
            i += 1
    return datae,datar

def plot_data(datae,datar):
    """Plots the threads vs. time."""
    plt.plot(datae["threads"], datae["time"], marker='o', linestyle='-', label = "red_black")
    plt.plot(datar["threads"], datar["time"], marker='o', linestyle='-', label = "diagonal")
    plt.title('Threads vs Time')
    plt.xlabel('Threads')
    plt.ylabel('Time')
    plt.grid(True)
    plt.legend()
    plt.savefig(sys.argv[1])
    plt.close()

datae,datar = read_file("prof.log")
plot_data(datae,datar)