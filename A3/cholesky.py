import matplotlib.pyplot as plt

# Data
N = [10, 100, 400, 1000, 3000]
Tp = [9.91E-05, 5E-04, 2.58E-03, 0.0105, 0.28]#0.42]
Ts = [9.53E-07, 1.5E-04, 10.7E-03, 0.18,9.23]#13.87]
speedup = [9.4E-03, 0.3, 4.15, 17.14, 33.02]

# Plotting N vs serial_time and N vs parallel_time
plt.figure(figsize=(10, 5))
plt.plot(N, Ts, marker='o', label='Serial Time')
plt.plot(N, Tp, marker='o', label='Parallel Time')
plt.title('N vs Serial Time and Parallel Time')
plt.xlabel('N')
plt.ylabel('Time')
plt.legend()
plt.grid(True)
plt.savefig("Times.png")

# Plotting N vs speedup
plt.figure(figsize=(10, 5))
plt.plot(N, speedup, marker='o', color='r')
plt.title('N vs Speedup')
plt.xlabel('N')
plt.ylabel('Speedup')
plt.grid(True)
plt.savefig("speedup.png")
