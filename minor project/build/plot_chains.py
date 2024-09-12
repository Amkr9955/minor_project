import pandas as pd
import matplotlib.pyplot as plt

# Function to plot monomers in the ZY plane for a specific X value (interchanging Y and Z axes)
def plot_zy_plane(chain, x_value, title, color, filename):
    # Filter the chain for the specified X value
    zy_chain = chain[chain['X'] == x_value]
    
    # Check if there is any data for the given x_value
    if zy_chain.empty:
        print(f"No data available for X = {x_value}")
        return
    # Create a 2D plot for the YZ plane
    plt.figure(figsize=(10, 10), dpi=80)
    
    # Plot the points (monomers) and connect them with lines
    plt.plot(zy_chain['Z'], zy_chain['Y'], color=color, marker='s', linestyle='-', 
             markerfacecolor='none', markersize=10)
    
    # Set grid and limits for better visibility
    plt.grid(True, linestyle=':', color='blue', alpha=0.5)
    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xticks(range(0, 50, 1))   # Adjusting the range based on Z values
    plt.yticks(range(0, 50, 1))  # Adjusting the range based on Y values

    # Set margins to prevent the graph from looking squished
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    
    
    # Add labels and title (Note the swap of axes in the labels)
    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title(title)
    
    # Save and show the plot
    plt.savefig(filename, bbox_inches='tight', dpi=80)
    plt.show()

# Load the data from CSV
data = pd.read_csv("/home/amitkumar/Documents/minor project/build/initial_chains.csv")
data_final = pd.read_csv("/home/amitkumar/Documents/minor project/build/final_chains.csv")


# Separate the chains
chain1 = data[data['Chain'] == 'Chain1']
chain2 = data[data['Chain'] == 'Chain2']

chain1_final = data_final[data_final['Chain'] == 'Chain1']
chain2_final = data_final[data_final['Chain'] == 'Chain2']

#plot_zy_plane(chain1, 0, "Chain 1 ZY Conformation (X=0)", 'green', "chain1_zy_x0.png")
plot_zy_plane(chain1, 0, "Chain 1 ZY Conformation (X=0)", 'green', "chain1_initial_x0.png")
plot_zy_plane(chain2, 2, "Chain 2 ZY Conformation (X=2)", 'red', "chain2_intial_x2.png")


# Plot YZ plane for X=0 and X=2 (Final)
plot_zy_plane(chain1_final, 0, "Chain 1 ZY Conformation (X=0)", 'green', "chain1_final_x0.png")
plot_zy_plane(chain2_final, 2, "Chain 2 ZY Conformation (X=2)", 'red', "chain2_final_x2.png")