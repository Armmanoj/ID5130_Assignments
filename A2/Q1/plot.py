import matplotlib.pyplot as plt
import sys
import numpy as np

def read_data_from_file(filename):
    x_values = []
    y_values = [[] for k in range(3)]

    with open(filename, 'r') as file:
        line_number = 0
        for line in file:
            values = line.split()
            for i, value in enumerate(values):
                if line_number == 0:
                    x_values.append(i*0.002)
                y_values[line_number].append(float(value))
            line_number += 1

    return x_values, y_values

def plot_data(x_values, y_values):
    fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)  # Create a figure with 3 subplots vertically
    a,b,c = [],[],[]
    for i in range(len(x_values)):
        x = 2*i/len(x_values)
        if (x<0.5):
            a.append(np.sin(4*np.pi*(x)))
        else:
            a.append(0)
        if ((x>0.5) & (x<1)):
            b.append(np.sin(4*np.pi*(x-0.5)))
        else:
            b.append(0)
        if ((x>1) & (x<1.5)):
            c.append(np.sin(4*np.pi*(x-1)))
        else:
            c.append(0)
    axs[0].plot(x_values, y_values[0], label='numerical')
    axs[0].plot(x_values,a,label = "analytical")
    axs[0].set_title('t=0.0')

    axs[1].plot(x_values, y_values[1], label='numerical')
    axs[1].plot(x_values,b,label = "analytical")
    axs[1].set_title('t=0.5')
    axs[2].plot(x_values, y_values[2], label='numerical')
    axs[2].plot(x_values,c,label = "analytical")
    axs[2].set_title('t=1.0')

    for ax in axs:
        ax.set_ylabel('u')
        ax.legend()
        ax.grid(True)

    plt.xlabel('x')  # Set x-axis label once for all subplots
    plt.tight_layout()  # Adjust subplot parameters to give specified padding
    plt.suptitle(sys.argv[1])  # Set the title for the entire figure
    plt.title(sys.argv[1])
    plt.savefig(sys.argv[1]+".png")

filename = "out.txt"
x_values, y_values = read_data_from_file(filename)
plot_data(x_values, y_values)
