import matplotlib.pyplot as plt
import sys

def read_file(filename):
    """Reads the file and extracts the data."""
    data = {"threads": [], "time": []}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            try:
                threads = int(parts[1])
                time = float(parts[2])
            except:
                continue
            data["threads"].append(threads)
            data["time"].append(time)
    return data

def plot_data(data):
    """Plots the threads vs. time."""
    plt.plot(data["threads"], data["time"], marker='o', linestyle='-')
    plt.title('Threads vs Time')
    plt.xlabel('Threads')
    plt.ylabel('Time')
    plt.grid(True)
    plt.savefig(sys.argv[1])
    plt.close()

data = read_file("profile.log")
plot_data(data)
