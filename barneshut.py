##################
## barneshut.py 
## ------------
## implements the Barnes-Hut quadtree n-body gravity solver.
## also implements direct pair-wise force calculation for 
## comparision.
## created: 09-31-2016
## modified:
## Christopher Mauney
##################


import sys, os, time
import numpy as np
from quadtree import cell
from quadtree import particle

G = 1.0             # Newton's constant. set to 1.0 for convenience
EPS  = 1.0E-1       # softening factor to prevent singularities
EPS2 = EPS * EPS    # softening term used in force calculations

def accel_xy(xi, yi, xj, yj, mj):
    dx = xj - xi
    dy = yj - yi
    deno = np.power(dx * dx + dy * dy + EPS2, 1.5)
    fact = G * mj / deno
    ax = fact * dx
    ay = fact * dy

    return ax, ay

def tf(cell, particle, threshold):
    ax, ay = 0.0, 0.0
    if cell.check_threshold(particle, threshold):
        ax, ay = accel_xy(particle.x, particle.y, cell.xcom, cell.ycom, cell.mcom)
        return ax, ay
    else:
        for child in cell.children:
            cax, cay = tf(child, particle, threshold)
            ax = ax + cax
            ay = ay + cay

    return ax, ay

def tree_force_recursive(particles, xmin, xmax, ymin, ymax, threshold):
    N = len(particles)
    a = np.zeros((N,2))

    root = cell(xmin, xmax, ymin, ymax, 0)
    for particle in particles:
        root.add(particle)

    for particle in root.particles():
        ax, ay = tf(root, particle, threshold)
        a[particle.pid,0] = ax
        a[particle.pid,1] = ay
        
    return a

def tree_force_iterative(particles, xmin, xmax, ymin, ymax, threshold):
    N = len(particles)
    a = np.zeros((N,2))

    root = cell(xmin, xmax, ymin, ymax, 0)

    for particle in particles:
        root.add(particle)

    for particle in particles:
        cells = [root]
        while cells:
            c = cells.pop()
            if c.check_threshold(particle, threshold):
                if c.np > 0:
                    ax, ay = accel_xy(particle.x, particle.y, c.xcom, c.ycom, c.mcom)
                    a[particle.pid,0] += ax
                    a[particle.pid,1] += ay
            else:
                cells.extend(c.children)

    return a

def pair_force(particles):
    N = len(particles)
    a = np.zeros((N,2))
    for i in range(N):
        for j in range(i+1,N):
            ax, ay = accel_xy(particles[i].x, particles[i].y, particles[j].x, particles[j].y, particles[j].mass)
            a[i,0] += ax
            a[i,1] += ay
            a[j,:] -= (particles[i].mass / particles[j].mass * a[i,:])
    return a

# get input
if len(sys.argv) != 5:
    print 'incorrect usage'
    print '{0} nparticle nstep timestep method[p/ti/tr]'.format(sys.argv[0])
    sys.exit(1)

input_np = np.int(sys.argv[1])
input_ns = np.int(sys.argv[2])
input_ts = np.float(sys.argv[3])
input_mc = sys.argv[4]

if input_mc != 'p' and input_mc != 'ti' and input_mc != 'tr':
    print 'method {0} unknown/not implemented'.format(input_mc)
    sys.exit(1)

# nbody paratmeters
xmin = -5.0
xmax = 5.0               
ymin = -5.0
ymax = 5.0

nsteps = input_ns
timestep = input_ts
threshold = 1.0
method = input_mc

# list of particles
lp = []

# setup initial conditions
sim_time = 0.0
nparticles = input_np
for i in range(nparticles):
    mass = 1.0
    x = xmin + np.random.random() * (xmax - xmin)
    y = ymin + np.random.random() * (ymax - ymin)
    lp.append(particle(x, y, 0.0, 0.0, mass, i))

force_time_total = 0.0
evolve_time_total = 0.0

print 'beginning nbody simulation'
print '# of particles = {0}'.format(nparticles)
print 'nsteps = {0}, timestep = {1}'.format(nsteps, timestep)

with open('output.dat','w') as outf:

    nbody_time_start = time.time()
    for i in xrange(nsteps):

        force_time_start = time.time()

        if method == 'tr':
            a = tree_force_recursive(lp,xmin,xmax,ymin,ymax,threshold)
        elif method == 'ti':
            a = tree_force_iterative(lp,xmin,xmax,ymin,ymax,threshold)
        elif method == 'p':
            a = pair_force(lp)

        force_time_end = time.time()

        evolve_time_start = time.time()

        for i, particle in enumerate(lp):
            particle.vx += a[i,0] * timestep
            particle.vy += a[i,1] * timestep 

        for particle in lp:
            particle.x += particle.vx * timestep
            particle.y += particle.vy * timestep

        evolve_time_end = time.time()

        outf.write('t = {0}\n'.format(sim_time))
        for particle in lp:
            outf.write('{0:<5} {1:10.4f} {2:10.4f} {3:10.4f}\n'.format(particle.pid, particle.mass, particle.x, particle.y))

        # update simulation time
        sim_time = sim_time + timestep

        # update timing totals
        force_time_total += (force_time_end - force_time_start)
        evolve_time_total += (evolve_time_end - evolve_time_start)

    nbody_time_end = time.time()

print 'nbody simulation ended.'
print 'total time = {0} s'.format(nbody_time_end - nbody_time_start)
print 'total time in force = {0}, average time in force = {1} s'.format(force_time_total, force_time_total / np.float(nsteps))
print 'total time in evolve = {0}, average time in evolve = {1} s'.format(evolve_time_total, evolve_time_total / np.float(nsteps))
