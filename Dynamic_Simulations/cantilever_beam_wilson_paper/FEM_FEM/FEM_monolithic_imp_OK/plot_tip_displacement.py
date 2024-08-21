import numpy as np
import matplotlib.pyplot as plt 

with open ('tip_y_displacement.json','r') as f:
    lines = f.readlines()
    list = [entry.strip() for entry in lines]

y_disp = []

for line in lines:
    split = line.split()
    print
    if '"DISPLACEMENT_Y":' in split:
        continue
    if "{" in split:
        continue
    if "[" in split:
        continue
    if "]" in split:
        continue
    if "}" in split:
        continue
    y_disp.append(float(split[0].replace(',', '')))


time_step = 0.035
end_time = 2.765

time = np.arange(0, end_time, time_step).tolist()

plt.plot(time, y_disp)  
plt.title('Tip y-displacement curve')
plt.xlabel('time [s]')
plt.ylabel('displacement [m]')

# Add a grid
plt.grid(True)

# Optional: Customize grid appearance (e.g., color, line style)
plt.grid(color='gray', linestyle='--', linewidth=0.5)

plt.show()