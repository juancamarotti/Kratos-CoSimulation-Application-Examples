import numpy as np
import matplotlib.pyplot as plt

# Read the file and parse the lines
with open('output_reaction.txt', 'r') as f:
    lines = f.readlines()

# Initialize lists
reac = []
total_reac = []

# Process each line in the file
for line in lines:
    split = line.split()
    if split[0] in {'1', '2','3'}:
        reac.append(float(split[1]))
    if split[0] == '3':  # Assuming '4' signifies the end of a block
        total_reac.append(sum(reac))
        reac.clear()  # Clear the list for the next block

# Define time parameters
end_time = 20
time_step = 0.01

# Create time array
time = np.arange(0.0, end_time, time_step).tolist()

# Ensure time and total_reac have the same length
if len(time) != len(total_reac):
    print("Warning: The length of time and total_reac do not match.")

# Plotting
print(total_reac)
print(len(total_reac))
plt.plot(time[:len(total_reac)], total_reac)  # Use only the matching length
plt.title('Time-Reaction Curve')
plt.xlabel('Time (s)')
plt.ylabel('Reaction')
plt.show()
