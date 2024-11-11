import pandas as pd
import matplotlib.pyplot as plt

# Load data from CSV file
data = pd.read_csv('/home/amitkumar/Documents/minor project/build/results.csv')  # Replace with your actual CSV filename

# Extract columns
temperature = data['Temperature']
energy = data['E']
energy_squared = data['Esquared']
distance_d1 = data['d1']
distance_d2 = data['d2']
heat_capacity = data['Cv']

# Plot T vs E (in fourth quadrant by setting y-axis to negative values only)
plt.figure()
plt.plot(temperature, energy, label='E', color='blue')
plt.xlabel('Temperature (T)')
plt.ylabel('Average Energy (E)')
plt.title('T vs E ')
plt.ylim(bottom=min(energy) , top=1)  # Set y-axis to show only negative values, extending slightly beyond the minimum
plt.grid(True)
plt.legend()
plt.savefig('T_vs_E.png')

# Plot T vs Esquared
plt.figure()
plt.plot(temperature, energy_squared, label='Esquared', color='orange')
plt.xlabel('Temperature (T)')
plt.ylabel('Average Energy Squared (Esquared)')
plt.title('T vs Esquared')
plt.grid(True)
plt.legend()
plt.savefig('T_vs_Esquared.png')

# Plot T vs d1
plt.figure()
plt.plot(temperature, distance_d1, label='d1', color='green')
plt.xlabel('Temperature (T)')
plt.ylabel('Average Distance d1')
plt.title('T vs d1')
plt.grid(True)
plt.legend()
plt.savefig('T_vs_d1.png')

# Plot T vs d2
plt.figure()
plt.plot(temperature, distance_d2, label='d2', color='red')
plt.xlabel('Temperature (T)')
plt.ylabel('Average Distance d2')
plt.title('T vs d2')
plt.grid(True)
plt.legend()
plt.savefig('T_vs_d2.png')

# Plot T vs Cv
plt.figure()
plt.plot(temperature, heat_capacity, label='Cv', color='purple')
plt.xlabel('Temperature (T)')
plt.ylabel('Heat Capacity (Cv)')
plt.title('T vs Cv')
plt.grid(True)
plt.legend()
plt.savefig('T_vs_Cv.png')

print("Plots saved as PNG files.")
