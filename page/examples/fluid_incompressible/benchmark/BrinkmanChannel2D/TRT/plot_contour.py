import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob

# Function to read data from a .res file
def read_res_file(file_path):
    # Load the file assuming columns are separated by spaces or tabs
    data = np.loadtxt(file_path)
    return data

def dynamic_label(x_pos, y_pos, label, offset=0.05):
    # Slightly offset y position to avoid overlap
    return plt.text(x_pos, y_pos + offset, label, fontsize=8, color="red", ha='left', va='center')

def draw2lines(file_path, offset=0):
    # Read the data from the file
    data = read_res_file(file_path)

    start_index = file_path.find('F0') + 3
    end_index = file_path.find("_p00000")

    if start_index != -1 and end_index != -1:
        F0_tag = file_path[start_index:end_index]
        print(f"Processing {file_path}...")

    # Extract X and Y values
    x_values = data[:, 1]  # X is the first column
    y_values_1 = data[:, 3]  # Y1 is the fourth column (index 3)
    y_values_2 = data[:, 4]  # Y2 is the seventh column (index 6)
    x_pos = np.max(x_values)
    y_pos = np.max(np.concatenate((y_values_1, y_values_2)))

    # Plot the data
    plt.plot(x_values, y_values_2, 'r-', label='Reference')  # Red line
    plt.plot(x_values, y_values_1, 'b--', label='Computation')  # Blue circle points
    dynamic_label(x_pos+0.005, y_pos, r"$F_0$="+str(F0_tag), offset)

if __name__ == "__main__":
    # Directory containing .res files
    directory = 'tracking/brinkman_F0_*.res'

    # Get all .res files in the directory
    res_files = glob.glob(directory)

    plt.figure(figsize=(8, 6))
    plt.xlim(left=None, right=2.1)

    offset = [0, 0, 0, 0, 0.002, -0.002, 0]
    for index, file_path in enumerate(res_files):
        draw2lines(file_path, offset[index])

        
    red_line = Line2D([0], [0], color='red', lw=2)
    blue_line = Line2D([0], [0], color='blue', linestyle=':', lw=2)
    plt.legend([blue_line, red_line], ['Computation', 'Reference'], loc='upper left', frameon=True)

    plt.grid(False)
    plt.xlabel('Y', fontsize=12)
    plt.ylabel(r'$U_x$', fontsize=12)

    # Save the figure
    plt.savefig('Lambda_0.25.png', dpi=300, bbox_inches='tight')
    plt.close()