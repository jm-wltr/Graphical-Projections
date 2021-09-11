#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 14:22:08 2021

@author: jaimewalter
"""

## Visualize a cube that rotates in 3D, and its 2D projection. Made with linear transformations.
## This file only generates the coordinates of the vertices of the rotated cube and plots them in 3D.

import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.pyplot as plt

# Coordinates of vertices
points = np.array([[-1, -1, -1],
                  [1, -1, -1 ],
                  [1, 1, -1],
                  [-1, 1, -1],
                  [-1, -1, 1],
                  [1, -1, 1 ],
                  [1, 1, 1],
                  [-1, 1, 1]])
                       
# Define rotation matrices
phi = math.pi / 9
theta = 19 * math.pi / 90
psi = 2 * math.pi /45

XrollMatrix = np.mat([[1, 0, 0], \
                      [0, math.cos(phi), math.sin(phi)], \
                      [0, -math.sin(phi), math.cos(phi)]])
    
YpitchMatrix = np.mat([[math.cos(theta), 0, -math.sin(theta)], \
                      [0, 1, 0], \
                      [math.sin(theta), 0, math.cos(theta)]])
                       
ZyawMatrix = np.mat([[math.cos(psi), math.sin(psi), 0], \
                      [-math.sin(psi), math.cos(psi), 0], \
                      [0, 0, 1]])
                       
# Rotate the cube
Z = np.zeros((8,3))
for i in range(8): Z[i,:] = np.dot(points[i,:],XrollMatrix)
for i in range(8): Z[i,:] = np.dot(Z[i,:],YpitchMatrix)
for i in range(8): Z[i,:] = np.dot(Z[i,:],ZyawMatrix)


## Make the plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1,1,1])

#Plot the points
ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2])

# Plot the faces
verts = [[Z[0],Z[1],Z[2],Z[3]],
  [Z[4],Z[5],Z[6],Z[7]], 
  [Z[0],Z[1],Z[5],Z[4]], 
  [Z[2],Z[3],Z[7],Z[6]], 
  [Z[1],Z[2],Z[6],Z[5]],
  [Z[4],Z[7],Z[3],Z[0]]]

ax.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='blue', alpha=.25))

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])

ax.locator_params(axis='x', nbins=3)
ax.locator_params(axis='y', nbins=3)
ax.locator_params(axis='z', nbins=3)

plt.show()
