import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16}) 


with open('coupling_iterations.txt', 'r') as f:
#with open('output_reaction.txt', 'r') as f:
    lines = f.readlines()
    list = [entry.strip() for entry in lines]
    
coupling_iterations = []
    
for line in lines:
    split = line.split()
    coupling_iterations.append(float(split[0]))
    
average = sum(coupling_iterations)/len(coupling_iterations)

    	
#print(z_disp)
#print(len(z_disp))

time_step = 0.035
end_time = 2.765


time = np.arange(0.0, end_time, time_step).tolist()
#print(time)
#print(len(time))


plt.plot(time, coupling_iterations,label=r'IGA-IGA Partitioned')  
plt.axhline(y=average, color='r', linestyle='--', label=r'Average value')
plt.title('Number of coupling iterations for Partitioned IGA/IGA Simulation ($E_O/E_D=1$)')
plt.xlabel('time (s)')
plt.ylabel('Number of coupling iterations')


plt.legend()


plt.ylim(bottom=0, top=11)  # y-axis limits

# Add a grid
plt.grid(True)

# Optional: Customize grid appearance (e.g., color, line style)
plt.grid(color='gray', linestyle='--', linewidth=0.5)

#plt.title('time-displacement curve')
#plt.xlabel('time (s)')
#plt.ylabel('displacement z')

plt.show()  	    
