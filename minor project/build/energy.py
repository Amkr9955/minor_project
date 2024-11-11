import matplotlib.pyplot as plt

# File path 
file_path = '/home/amitkumar/Documents/minor project/build/energy_history.txt'

# Initialize lists for steps and energy values
steps = []
energy = []

# Reading the file and storing values
with open(file_path, 'r') as file:
    for line in file:
        step, e = map(int, line.strip().split(','))
        steps.append(step)
        energy.append(e)

# Filter steps and energy for every 100th step
filtered_steps = steps[::50]  # Take every 100th step
filtered_energy = energy[::50]  # Corresponding energy values

# Plotting
plt.figure(figsize=(12, 8))
plt.plot(filtered_steps, filtered_energy, marker='.', linestyle='-', color='b')
plt.xlabel('Steps', fontsize=12)
plt.ylabel('Energy', fontsize=12)
plt.title('Energy vs Steps (Plotted every 100 steps)', fontsize=16)
plt.grid(True)

# Show plot
plt.show()
