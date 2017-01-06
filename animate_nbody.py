import sys,os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

if len(sys.argv) != 3:
    print 'incorrect usage'
    sys.exit(1)

inputf = sys.argv[1]
frameskip = np.int(sys.argv[2])

time_list = []
pos_list = []

with open(inputf, 'r') as inpf:
    lines = inpf.readlines()
    
    lidx = 0
    
    while lidx < len(lines):

        if lines[lidx].startswith('t ='):
            time = np.float(lines[lidx].split()[2])
            time_list.append(time)
            pos_list.append([])

        else:
            tok = lines[lidx].split()
            mass, x, y = np.float(tok[1]), np.float(tok[2]), np.float(tok[3])
            pos_list[-1].append([mass, x, y])

        lidx = lidx + 1
    
times = np.asarray(time_list)
xypos = np.asarray(pos_list)

fig, ax = plt.subplots()

masses = xypos[0][:,0]
#scaled_masses = ( (masses - np.amin(masses)) / (np.amax(masses) - np.amin(masses)) + 1.0) * 50.0

scat = ax.scatter(xypos[0][:,1], xypos[0][:,2])

ax.set_title('{0}'.format(times[0]))

xlims = [-10.0, 10.0]
ylims = [-10.0, 10.0]

ax.set_xlim(xlims)
ax.set_ylim(ylims)

def update(i):
    N = xypos[i][:,0].size
    xy = np.zeros((N,2))
    for j in range(N):
        xy[j,0] = xypos[i][:,1][j]
        xy[j,1] = xypos[i][:,2][j]
       
    scat.set_offsets(xy)
    ax.set_title('{0}'.format(times[i]))

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)


    return scat

anim = animation.FuncAnimation(fig, update, frames=np.arange(0,times.size,frameskip), interval=1)
plt.show()
    


