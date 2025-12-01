import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
import pandas as pd


def draw_contour(filename):
    # Hyper parameters for the contour comparison
    threshold = 0.105


    # Load the data from the .res file assuming space-separated values
    data = pd.read_csv(filename, delim_whitespace=True, skiprows=2, header=None, names=['x', 'y', 'z', 'calculated', 'real'])


    # Extract x, y, and value columns
    x = data['x'].values
    y = data['y'].values
    values = data['calculated'].values
    reference = data['real'].values

    # Create a grid for interpolation
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]

    # Interpolate the values over the grid
    grid_values = griddata((x, y), values, (grid_x, grid_y), method='cubic')
    grid_ref = griddata((x, y), reference, (grid_x, grid_y), method='cubic')

    # Plot the contour where value equals the threshold
    
    plt.contour(grid_x, grid_y, grid_values, levels=[threshold], colors='blue', linestyles='dotted')
    plt.contour(grid_x, grid_y, grid_ref, levels=[threshold], colors='red', alpha=0.5)


if __name__ == '__main__':
    import sys
    plt.figure(figsize=(8, 8))
    filename = sys.argv[1]
    draw_contour(filename)
    red_line = Line2D([0], [0], color='red', lw=2)
    blue_line = Line2D([0], [0], color='blue', linestyle=':', lw=2)
    plt.legend([blue_line, red_line], ['Simulation', 'Reference'], loc='lower right', frameon=False)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)


    # Save the figure
    plt.savefig('contour.png', dpi=300, bbox_inches='tight')

