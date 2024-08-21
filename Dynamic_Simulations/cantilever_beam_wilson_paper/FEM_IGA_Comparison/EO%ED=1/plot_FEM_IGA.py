import numpy as np
import matplotlib.pyplot as plt 

# Process the FEM Monolithic data
with open ('tip_y_displacement_monolithic.json','r') as f_fem_mono:
    lines_fem_mono = f_fem_mono.readlines()

y_disp_fem_mono = []

for line in lines_fem_mono:
    split_fem_mono = line.split()
    if '"DISPLACEMENT_Y":' in split_fem_mono:
        continue
    if "{" in split_fem_mono:
        continue
    if "[" in split_fem_mono:
        continue
    if "]" in split_fem_mono:
        continue
    if "}" in split_fem_mono:
        continue
    y_disp_fem_mono.append(float(split_fem_mono[0].replace(',', '')))

# Process the IGA monolithic data
with open('output_monolithic.txt', 'r') as f_iga_mono:
    lines_iga_mono = f_iga_mono.readlines()
    list_iga_mono = [entry.strip() for entry in lines_iga_mono]
    
y_disp_iga_mono = []
    
for line in lines_iga_mono:
    split_iga_mono = line.split()
    if(split_iga_mono[0] == '2'):
    	y_disp_iga_mono.append(float(split_iga_mono[2]))

# Process the IGA partitioned data
with open('output.txt', 'r') as f_iga_part:
    lines_iga_part = f_iga_part.readlines()
    list_iga_part = [entry.strip() for entry in lines_iga_part]
    
y_disp_iga_part = []
    
for line in lines_iga_part:
    split_iga_part = line.split()
    if(split_iga_part[0] == '2'):
    	y_disp_iga_part.append(float(split_iga_part[2]))


# Process the FEM Partitioned data
with open ('tip_y_displacement_partitioned.json','r') as f_fem_part:
    lines_fem_part = f_fem_part.readlines()

y_disp_fem_part = []

for line in lines_fem_part:
    split_fem_part = line.split()
    if '"DISPLACEMENT_Y":' in split_fem_part:
        continue
    if "{" in split_fem_part:
        continue
    if "[" in split_fem_part:
        continue
    if "]" in split_fem_part:
        continue
    if "}" in split_fem_part:
        continue
    y_disp_fem_part.append(float(split_fem_part[0].replace(',', '')))


time_step = 0.035
end_time = 2.765

time = np.arange(0, end_time, time_step).tolist()

plt.plot(time, y_disp_fem_mono, label='FEM Monolithic')  
plt.plot(time, y_disp_iga_part, label='IGA Partitioned')  
plt.plot(time, y_disp_fem_part, label='FEM Partitioned') 
plt.plot(time, y_disp_iga_mono, label='IGA Monolithic')   
plt.title('y-displacement curve of the Cantilever Beam Tip')
plt.xlabel('time [s]')
plt.ylabel('v [m]')
plt.legend()

plt.ylim(bottom=-0.30, top=0.02)  # y-axis limits

# Add a grid
plt.grid(True)

# Optional: Customize grid appearance (e.g., color, line style)
plt.grid(color='gray', linestyle='--', linewidth=0.5)

plt.show()
