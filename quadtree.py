##################
## quadtree.py
## ----------- 
## implements a custom quadtree data structure 
## to use in the Barnes-Hut n-body algorithm
## created: 09-31-2016
## modified: 
## Christopher Mauney
##################

import numpy as np

class particle(object):
    '''
        particle represents the dynamical objects with 
        gravitational mass.
    '''
    def __init__(self, x, y, vx, vy, mass, pid):
        self.x, self.y = x, y       # positions
        self.vx, self.vy = vx, vy   # velocities
        self.mass = mass            # mass
        self.pid = pid              # unique identifier

class cell(object):
    '''
        cell is a tree node representing the quadrature.
        it will hold the center of mass (com) information
        as well as being a quadtree.
    '''
    def __init__(self, xmin, xmax, ymin, ymax, cid):
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.cid = cid

        # cells are created with no particles
        # com data is updated when particles get added.
        self.np = 0
        self.mcom = 0.0E0
        self.xcom, self.ycom = 0.0E0, 0.0E0

        self.children = []
        self.particle = None

    def point_in_cell(self, x, y):
        '''
            checks if the point (x,y) is in this cell
            using the range open [xmin, xmax)
        '''
        if x >= self.xmin and x < self.xmax and \
            y >= self.ymin and y < self.ymax:
            return True

        return False

    def add(self, particle):
        '''
            recursively add this particle to the tree
        '''
        # ignore points outside of cell
        if not self.point_in_cell(particle.x, particle.y): return

        # if this subtree is non-empty, recurse down into an emtpy leaf
        if self.np > 0:
            # if this subtree is a 'leaf' (no branches) with data, 
            # create child nodes and recurse the data
            # into the appropriate child node
            if self.np == 1:
                self.make_children()
                for child in self.children:
                    child.add(self.particle)
                # this subtree is no longer a leaf, so clear the particle
                self.particle = None
            # add new particle to subtree
            # note that, at this point, we always have children
            for child in self.children:
                child.add(particle)
        # otherwise, this subtree is an empty leaf, so add it
        else:
            self.particle = particle
        # update center of mass of the cell
        # this is setup as a running average
        oldm = self.mcom
        self.mcom += particle.mass
        self.xcom = (oldm * self.xcom + particle.x * particle.mass) / self.mcom
        self.ycom = (oldm * self.ycom + particle.y * particle.mass) / self.mcom
        # increment the number of particles in subtree
        self.np = self.np + 1
    
    def make_children(self):
        '''
            insert 4 children into this tree
        '''
        # find the midpoints of this cell
        xmid = self.xmin + 0.5 * (self.xmax - self.xmin)
        ymid = self.ymin + 0.5 * (self.ymax - self.ymin)

        # attach new cells
        self.children.append(cell(self.xmin, xmid, ymid, self.ymax, 0))
        self.children.append(cell(xmid, self.xmax, ymid, self.ymax, 0))
        self.children.append(cell(xmid, self.xmax, self.ymin, self.ymax, 0))
        self.children.append(cell(self.xmin, xmid, self.ymin, ymid, 0))

    def particles(self):
        '''
            recursively returns a list of particles in a depth-first
            approach.
        '''
        # if there is particle data, this is a leaf
        if self.particle: return [self.particle]
        # add child nodes to the list
        if self.children:
            plist = []
            for child in self.children:
                plist.extend(child.particles())
            return plist

        return []

    def check_threshold(self, particle, threshold):
        '''
            the criteria (D / r < theta) is used to check if we take the com
            as the source, or if we need to go deeper into the tree.
            note: if theta -> 0.0, then we are simply doing pair-wise n-body
        '''
        # if this is a node with children (and thus has com), determine
        # wether the algorithm should use the node as a gravitational source
        if self.children:
            d = self.xmax - self.xmin
            dx = particle.x - self.xcom
            dy = particle.y - self.ycom
            r  = np.sqrt(dx * dx + dy * dy)

            # return truth value of the criteria 
            return (d / r) < threshold
        # otherwise, the node is just a particle, and we will perform pairwise
        # gravity calculation provided that the two particles are seperate
        else:
            return self.particle != particle
