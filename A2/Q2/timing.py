import matplotlib.pyplot as plt
import numpy as np

# Read timing data from the file
with open('timing.txt', 'r') as file:
    timings = np.array([float(line.strip()) for line in file])

# Define x-axis values for the first 4 data points
x_values = [2,  4, 8,10]

# Plot the first 4 data points
plt.plot(x_values, (timings[:4]), label='red_black')

# Plot the next 4 data points on the same graph
plt.plot(x_values, (timings[4:]), label='jacobi')

# Add labels and title
plt.xlabel('Processor Count')
plt.ylabel('Run time')
plt.title('')
plt.legend()

# Display the plot
plt.savefig("timing.png")
plt.close()
plt.plot(x_values, 4*(timings[1]-timings[2])/(timings[:4]), label='red_black')

# Plot the next 4 data points on the same graph
plt.plot(x_values, 4*(timings[5]-timings[6])/(timings[4:]), label='jacobi')

# Add labels and title
plt.xlabel('Processor Count')
plt.ylabel('Speed-up')
plt.title('')
plt.legend()

# Display the plot
plt.savefig("speedup.png")
plt.close()
