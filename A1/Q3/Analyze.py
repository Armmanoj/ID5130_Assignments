import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

# needs N+1 N+1 image_filename as commnd line input
# Function to read the matrix from the binary file
def read_matrix(filename, rows, cols):
    with open(filename, "rb") as f:
        matrix = np.fromfile(f, dtype=np.float64)
        matrix = matrix.reshape((rows, cols))
    return matrix

# Function to plot matrix as colormap
def plot_colormap(matrix,err):
    plt.imshow(matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.title('Matrix Colormap error is {:.2f} %'.format(100*round(err,4)))
    plt.savefig(sys.argv[3])
    plt.close()

def plot_2d(matrix, m):
    x = np.linspace(-1, 1, 21)
    # Plot analytical solution
    plt.plot(x, matrix[:, 5].transpose(), label='Analytical')
    # Plot numerical solution
    plt.plot(x, m[:, 5].transpose(), label='Numerical')
    # Set plot title
    plt.title('Part a) Analytical vs Numerical')
    # Set axis labels
    plt.xlabel('x')
    plt.ylabel('Solution')
    # Add legend
    plt.legend()
    # Save plot to file
    plt.savefig(sys.argv[5])
    plt.close()
# Function to plot matrix as 3D graph
def plot_3d(matrix,m,err):
    fig = plt.figure(figsize=(10, 5))
    # Plotting numerical solution
    ax = fig.add_subplot(121, projection='3d')
    x, y = np.meshgrid(np.arange(matrix.shape[1]), np.arange(matrix.shape[0]))
    ax.plot_surface(x, y, matrix, cmap='viridis')
    ax.set_title('Numerical solution 3D Plot')
    # Plotting error with the analytical solution
    ax1 = fig.add_subplot(122, projection='3d')
    x1, y1 = np.meshgrid(np.arange(m.shape[1]), np.arange(m.shape[0]))
    ax1.plot_surface(x1, y1, m-matrix, cmap='viridis')
    ax1.set_title('Error in solution (analytic-numerical) 3D Plot')

    plt.tight_layout()
    plt.savefig(sys.argv[4])
    plt.close()

def main():
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    filename = "out.bin"
    # create theoretical solution
    m = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            m[i,j]=(1-(2*i/(rows-1)-1)*(2*i/(rows-1)-1))*(1-(2*j/(cols-1)-1)*(2*j/(cols-1)-1))
    # Read solutionfrom binary file
    matrix = read_matrix(filename, rows, cols)
    err = 0
    for i in range(rows):
        for j in range(cols):
            err = max(np.abs(matrix[i,j]-m[i,j]),err)
    # Plot solution as colormap
    plot_colormap(matrix, err)

    # Plot solution as 3D graph
    plot_3d(matrix,m, err)
    plot_2d(matrix,m)

if __name__ == "__main__":
    main()

    
