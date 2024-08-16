import numpy as np
import matplotlib.pyplot as plt


with open('output.txt', 'r') as f:
#with open('output_reaction.txt', 'r') as f:
    lines = f.readlines()
    list = [entry.strip() for entry in lines]
    
z_disp = []
    
for line in lines:
    split = line.split()
    if(split[0] == '1'):
    	z_disp.append(float(split[1]))
    	
#print(z_disp)
#print(len(z_disp))

time_step = 0.0001
end_time = 61636*time_step


time = np.arange(0, end_time, time_step).tolist()
#print(time)
#print(len(time))


plt.plot(time, z_disp)  
plt.title('time-displacement curve')
plt.xlabel('time (s)')
plt.ylabel('reaction')

#plt.title('time-displacement curve')
#plt.xlabel('time (s)')
#plt.ylabel('displacement z')

plt.show()  	    
