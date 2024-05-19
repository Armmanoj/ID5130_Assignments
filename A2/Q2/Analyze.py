import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
# needs N+1 N+1 image_filename as commnd line input
# Function to read the matrix from the binary file
def read_matrix(filename, rows, cols):
    print("His")
    with open(filename, "rb") as f:
        matrix = np.fromfile(f, dtype=np.float64)
        matrix = matrix.reshape((rows, cols))
    return matrix

# Function to plot matrix as 3D graph
def plot_3d(matrix,err, N):
    fig = plt.figure(figsize=(10, 5))
    # Plotting numerical solution
    ax = fig.add_subplot(121, projection='3d')
    x, y = np.meshgrid(np.arange(matrix.shape[1])/((N-1)/2)-1, np.arange(matrix.shape[0])/((N-1)/2)-1)
    ax.plot_surface(x, y, matrix, cmap='viridis')
    ax.set_title('Numerical solution 3D Plot')
    # Plotting error with the analytical solution
    ax1 = fig.add_subplot(122, projection='3d')
    x1, y1 = np.meshgrid(np.arange(err.shape[1])/((N-1)/2)-1, np.arange(err.shape[0])/((N-1)/2)-1)
    ax1.plot_surface(x1, y1, err, cmap='viridis')
    ax1.set_title('Error in solution 3D Plot')

    plt.tight_layout()
    plt.savefig("3d.png")
    plt.close()

def plot_mid(matrix):
    # Get the number of rows and columns in the matrix
    rows, cols = matrix.shape
    
    # Create meshgrid for plotting
    x, y = np.meshgrid(2*np.arange(cols)/(cols-1)-1, 2*np.arange(rows)/(rows-1)-1)
    
    # Extract the middle column
    middle_col = matrix[:, cols // 2]

    # Extract the middle row
    middle_row = matrix[rows // 2, :]

    # Create figure and axes
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Plot middle column
    axes[0].plot(y[:, cols // 2], middle_col)
    axes[0].set_title('x=0')

    # Plot middle row
    axes[1].plot(x[rows // 2, :], middle_row)
    axes[1].set_title('y=0')

    # Show plot
    plt.savefig("2da.png")

import numpy as np
import matplotlib.pyplot as plt

def plot_middle(matrix, matrix1):
    
    # Get the number of rows and columns in the matrix
    rows, cols = matrix.shape
    
    # Create meshgrid for plotting
    x, y = np.meshgrid(2*np.arange(cols)/(cols-1)-1, 2*np.arange(rows)/(rows-1)-1)
    
    # Extract the middle column
    middle_col = matrix[:, cols // 2]

    # Extract the middle row
    middle_row = matrix[rows // 2, :]
    
    # Create figure and axes
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Plot middle column
    axes[0].plot(y[:, cols // 2], middle_col, label='File 1')
    
    # Extract the middle column from the second file
    middle_col1 = matrix1[:, cols // 2]
    
    # Plot middle column from the second file on top of the first one
    axes[0].plot(y[:, cols // 2], middle_col1, label='File 2')
    
    axes[0].set_title('x=0')
    axes[0].legend()

    # Plot middle row
    axes[1].plot(x[rows // 2, :], middle_row, label='File 1')
    
    # Extract the middle row from the second file
    middle_row1 = matrix1[rows // 2, :]
    
    # Plot middle row from the second file on top of the first one
    axes[1].plot(x[rows // 2, :], middle_row1, label='File 2')
    
    axes[1].set_title('y=0')
    axes[1].legend()

    # Show plot
    plt.savefig(sys.argv[3])  # Save the plot as 2da.png


def main():
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    filename = sys.argv[4]
    filename1 = sys.argv[5]
    errname = "err.bin"
    # create theoretical solution
    err = read_matrix(errname, rows, cols)
    max_error = np.max(err)
    # Read solutionfrom binary file
    matrix = read_matrix(filename, rows, cols)
    matrix1 = read_matrix(filename1, rows, cols)
    plot_middle(matrix,matrix1)
    # Plot solution as 3D graph
    plot_3d(matrix, err,rows)
if __name__ == "__main__":
    main()

    
