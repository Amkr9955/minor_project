import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Function to plot the 3D chains with very small blue dots on the lattice grid and square monomers
def plot_3d_chains_with_lattice(chain1, chain2, filename):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Define the Y and X range for the wall at Z=0
    y_vals = np.linspace(0, 80, 100)  # Y-axis range
    x_vals = np.linspace(-2, 4, 100)   # X-axis range
    Y, X = np.meshgrid(y_vals, x_vals)  # Create grid for Y and X
    Z = np.zeros_like(X)  # Z=0 for the wall

    # Plot the dark grey wall at Z=0
    ax.plot_surface(X, Z, Y, color='darkgray', alpha=0.7)

    # Plot Chain 1 (X = 0) with square markers
    ax.plot(chain1['X'], chain1['Z'], chain1['Y'], color='green', marker='s', linestyle='-', 
            markerfacecolor='none', markersize=8, label='Chain 1 ZY Conformation (X=0)')

    # Plot Chain 2 (X = 2) with square markers
    ax.plot(chain2['X'], chain2['Z'], chain2['Y'], color='red', marker='s', linestyle='-', 
            markerfacecolor='none', markersize=8, label='Chain 2 ZY Conformation (X=2)')

    # Add very small blue dots at the lattice points for X=0 and X=2
    for y in range(0, 80, 4):
        for z in range(0, 80, 4):
            ax.scatter(0, z, y, color='blue', s=1)  # Very small blue dots at lattice points
            ax.scatter(2, z, y, color='blue', s=1)  # Very small blue dots at lattice points

    # Set labels and title
    ax.set_xlabel('X', fontsize=15)
    ax.set_ylabel('Z', fontsize=15)
    ax.set_zlabel('Y', fontsize=15)
    ax.set_title("3D Conformation Lattice Structure", fontsize=18)

    # Set limits to match the image provided
    ax.set_xlim([-2, 4])
    ax.set_ylim([0, 80])
    ax.set_zlim([0, 80])

    # Customize the tick marks for the lattice grid
    ax.set_xticks(range(-2, 4, 1))
    ax.set_yticks(range(0, 80, 4))
    ax.set_zticks(range(0, 80, 4))

    #  Customize axis grid lines
    ax.xaxis._axinfo["grid"].update(linewidth=0.2, linestyle=':', alpha=0.5)  
    ax.yaxis._axinfo["grid"].update(linewidth=0.2, linestyle=':', alpha=0.5) 
    ax.zaxis._axinfo["grid"].update(linewidth=0.2, linestyle=':', alpha=0.5)  

    # Highlight grid lines at x=0 and x=2 by plotting them manually with thicker, darker lines
    ax.plot([0, 0], [0, 80], [0, 0], color='black', linewidth=2, linestyle='--')  # Darker grid line at X=0
    ax.plot([2, 2], [0, 80], [0, 0], color='black', linewidth=2, linestyle='--')  # Darker grid line at X=2

    # Set box aspect ratio (ensure no indentation error here)
    ax.set_box_aspect([2, 3, 3])  # Adjust aspect ratio for proper scaling

    # Rotate the view to a desired angle while keeping the coordinates fixed
    ax.view_init(elev=30, azim=45)  # Adjust these parameters as needed

    # Show the legend
    ax.legend()

    # Save and show the plot
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.show()

# Function to plot monomers in the ZY plane for a specific X value (interchanging Y and Z axes)
def plot_zy_plane(chain, x_value, title, color, filename):
    # Filter the chain for the specified X value
    zy_chain = chain[chain['X'] == x_value]
    
    # Check if there is any data for the given x_value
    if zy_chain.empty:
        print(f"No data available for X = {x_value}")
        return

    # Create a 2D plot for the YZ plane
    plt.figure(figsize=(15, 15), dpi=100)
    
    # Plot the points (monomers) and connect them with lines
    plt.plot(zy_chain['Z'], zy_chain['Y'], color=color, marker='s', linestyle='-', 
             markerfacecolor='none', markersize=10)
    
    # Set grid and limits for better visibility
    plt.grid(True, linestyle=':', color='blue', alpha=1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xticks(range(0, 80, 1))  # Adjusting the range based on Z values
    plt.yticks(range(0, 80, 1))  # Adjusting the range based on Y values

    # Set margins to prevent the graph from looking squished
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    
    # Add labels and title (Note the swap of axes in the labels)
    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title(title)
    
    # Save and show the plot
    plt.savefig(filename, bbox_inches='tight', dpi=100)
    plt.show()

# Load the data from CSV
data_initial = pd.read_csv("/home/amitkumar/Documents/minor project/build/initial_chains.csv")
data_final_equill = pd.read_csv("/home/amitkumar/Documents/minor project/build/equill_final_chains.csv")
data_final = pd.read_csv("/home/amitkumar/Documents/minor project/build/final_chains.csv")

# Separate the chains
chain1_initial = data_initial[data_initial['Chain'] == 'Chain1']
chain2_initial = data_initial[data_initial['Chain'] == 'Chain2']

chain1_final_equill = data_final_equill[data_final_equill['Chain'] == 'Chain1']
chain2_final_equill= data_final_equill[data_final_equill['Chain'] == 'Chain2']

chain1_final = data_final[data_final['Chain'] == 'Chain1']
chain2_final= data_final[data_final['Chain'] == 'Chain2']

# Update X values to fixed points for 3D plotting
chain1_initial['X'] = 0  # Chain 1 at X=0
chain2_initial['X'] = 2  # Chain 2 at X=2

chain1_final_equill['X'] = 0  # Chain 1 at X=0
chain2_final_equill['X'] = 2  # Chain 2 at X=2

chain1_final['X'] = 0  # Chain 1 at X=0
chain2_final['X'] = 2  # Chain 2 at X=2

# Plot 3D views with very small blue dots and square monomers, with fixed coordinates
plot_3d_chains_with_lattice(chain1_initial, chain2_initial, "3d_lattice_chains_initial.png")
plot_3d_chains_with_lattice(chain1_final_equill, chain2_final_equill, "3d_lattice_chains_final_equillibrium.png")

# Plot ZY plane for initial chains
plot_zy_plane(chain1_initial, 0, "Chain 1 Initial ZY Conformation (X=0)", 'green', "chain1_initial_zy_x0.png")
plot_zy_plane(chain2_initial, 2, "Chain 2 Initial ZY Conformation (X=2)", 'red', "chain2_initial_zy_x2.png")

# Plot ZY plane for final chains
plot_zy_plane(chain1_final_equill, 0, "Chain 1 Equillibrium Final ZY Conformation (X=0)", 'green', "chain1_final_equill_zy_x0.png")
plot_zy_plane(chain2_final_equill, 2, "Chain 2 Equillibrium Final ZY Conformation (X=2)", 'red', "chain2_final_equill_zy_x2.png")

plot_3d_chains_with_lattice(chain1_final, chain2_final, "3d_lattice_chains_final.png")
plot_zy_plane(chain1_final, 0, "Chain 1 Final ZY Conformation (X=0)", 'green', "chain1_final_zy_x0.png")
plot_zy_plane(chain2_final, 2, "Chain 2 Final ZY Conformation (X=2)", 'red', "chain2_final_zy_x2.png")


