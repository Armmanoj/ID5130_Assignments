import matplotlib.pyplot as plt

# Data
gangs = [10, 100, 1000]
time = [84.283646, 69.621795, 106.124363]

# Plot
plt.plot(gangs, time, marker='o', linestyle='-')

# Labels and title
plt.xlabel('Gangs')
plt.ylabel('Time')
plt.title('Gangs vs Time')

# Show plot
plt.savefig("time2.png")
