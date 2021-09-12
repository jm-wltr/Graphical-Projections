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
from matplotlib.patches import Polygon
from matplotlib.widgets import Slider, Button

# Define rotation matrices
phi = 0
theta = 0
psi = 0

def rotate(phi, theta, psi):
    # Coordinates of vertices
    vertices = np.array([[-1, -1, -1],
                      [1, -1, -1 ],
                      [1, 1, -1],
                      [-1, 1, -1],
                      [-1, -1, 1],
                      [1, -1, 1 ],
                      [1, 1, 1],
                      [-1, 1, 1]])
    
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
    for i in range(8): Z[i,:] = np.dot(vertices[i,:],XrollMatrix)
    for i in range(8): Z[i,:] = np.dot(Z[i,:],YpitchMatrix)
    for i in range(8): Z[i,:] = np.dot(Z[i,:],ZyawMatrix)
    
    # Plot the edges and faces
    faces = [[Z[0],Z[1],Z[2],Z[3]],
      [Z[4],Z[5],Z[6],Z[7]], 
      [Z[0],Z[1],Z[5],Z[4]], 
      [Z[2],Z[3],Z[7],Z[6]], 
      [Z[1],Z[2],Z[6],Z[5]],
      [Z[4],Z[7],Z[3],Z[0]]]
    
    # Orthographic views
    ortXY = np.mat([[1, 0, 0], \
                   [0, 1, 0]])
    ortXZ = np.mat([[1, 0, 0], \
                    [0, 0, 1]])
    ortYZ = np.mat([[0, 1, 0], \
                    [0, 0, 1]])
    
    XY = np.zeros((8,2))
    for i in range(8): XY[i,:] = np.dot(ortXY, Z[i,:])
    #XY = np.delete(XY, 2, 1)

    XZ = np.zeros((8,2))
    for i in range(8): XZ[i,:] = np.dot(ortXZ, Z[i,:])
    #XZ = np.delete(XZ, 1, 1)

    YZ = np.zeros((8,2))
    for i in range(8): YZ[i,:] = np.dot(ortYZ, Z[i,:])
    #YZ = np.delete(YZ, 0, 1)
    
    facesXY = [[XY[7],XY[4],XY[5],XY[6]],
      [XY[1],XY[2],XY[3],XY[0]], 
      [XY[4],XY[0],XY[3],XY[7]], 
      [XY[5],XY[1],XY[2],XY[6]], 
      [XY[7],XY[3],XY[2],XY[6]],
      [XY[0],XY[4],XY[5],XY[1]]]
    facesXZ = [[XZ[7],XZ[4],XZ[5],XZ[6]],
      [XZ[1],XZ[2],XZ[3],XZ[0]], 
      [XZ[4],XZ[0],XZ[3],XZ[7]], 
      [XZ[5],XZ[1],XZ[2],XZ[6]], 
      [XZ[7],XZ[3],XZ[2],XZ[6]],
      [XZ[0],XZ[4],XZ[5],XZ[1]]]
    facesYZ = [[YZ[7],YZ[4],YZ[5],YZ[6]],
      [YZ[1],YZ[2],YZ[3],YZ[0]], 
      [YZ[4],YZ[0],YZ[3],YZ[7]], 
      [YZ[5],YZ[1],YZ[2],YZ[6]], 
      [YZ[7],YZ[3],YZ[2],YZ[6]],
      [YZ[0],YZ[4],YZ[5],YZ[1]]]
    
    return Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces
    

results = rotate(phi, theta, psi)
Z = results[0]
XY = results[1]
XZ = results[2]
YZ = results[3]
facesXY = results[4]
facesXZ = results[5]
facesYZ = results[6]
faces = results[7]

 ## Make the plot
fig = plt.figure()
ax = fig.add_subplot(221, projection='3d', proj_type = "ortho", azim = -45, elev = 35, title="Isometric")
ax.set_box_aspect([1,1,1])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])

ax.locator_params(axis='x', nbins=3)
ax.locator_params(axis='y', nbins=3)
ax.locator_params(axis='z', nbins=3)


ax2 = fig.add_subplot(222, title="Top")
ax2.set_aspect('equal', adjustable='box')
ax3 = fig.add_subplot(224, title="Front")
ax3.set_aspect('equal', adjustable='box')
ax4 = fig.add_subplot(223, title="Side")
ax4.set_aspect('equal', adjustable='box')

ax2.set_xlim([-2, 2])
ax2.set_ylim([-2, 2])
ax3.set_xlim([-2, 2])
ax3.set_ylim([-2, 2])
ax4.set_xlim([-2, 2])
ax4.set_ylim([-2, 2])

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax3.set_xlabel('X')
ax3.set_ylabel('Z')
ax4.set_xlabel('Y')
ax4.set_ylabel('Z')

plt.subplots_adjust(bottom=0.3)

# Show plot
plt.subplots_adjust(wspace = 0.5, hspace = 0.5)
plt.show()

def plot(Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces):
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    
    
    ax.set_box_aspect([1,1,1])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])
    
    ax.locator_params(axis='x', nbins=3)
    ax.locator_params(axis='y', nbins=3)
    ax.locator_params(axis='z', nbins=3)
    
    ax2.set_xlim([-2, 2])
    ax2.set_ylim([-2, 2])
    ax3.set_xlim([-2, 2])
    ax3.set_ylim([-2, 2])
    ax4.set_xlim([-2, 2])
    ax4.set_ylim([-2, 2])
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')
    ax4.set_xlabel('Y')
    ax4.set_ylabel('Z')
        
    
    
    ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2])

    ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=1, edgecolors='blue', alpha=.25))

    ax2.scatter(XY[:, 0], XY[:, 1])
    ax3.scatter(XZ[:, 0], XZ[:, 1])
    ax4.scatter(YZ[:, 0], YZ[:, 1])
    
    for i in range(6):
        poly = Polygon(facesXY[i], facecolor='cyan', edgecolor='blue', alpha = 0.25)
        ax2.add_patch(poly)
    for i in range(6):
        poly = Polygon(facesXZ[i], facecolor='cyan', edgecolor='blue', alpha = 0.25)
        ax3.add_patch(poly)
    for i in range(6):
        poly = Polygon(facesYZ[i], facecolor='cyan', edgecolor='blue', alpha = 0.25)
        ax4.add_patch(poly)
    
plot(Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces)
    
axphi = plt.axes([0.19, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')
phi_slider = Slider(
    ax=axphi,
    label='Phi [rad]',
    valmin=0,
    valmax=1.57,
    valinit=phi,
)

axtheta = plt.axes([0.19, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
theta_slider = Slider(
    ax=axtheta,
    label='Theta [rad]',
    valmin=0,
    valmax=1.57,
    valinit=theta,
)

axpsi = plt.axes([0.19, 0.05, 0.65, 0.03], facecolor='lightgoldenrodyellow')
psi_slider = Slider(
    ax=axpsi,
    label='Psi [rad]',
    valmin=0,
    valmax=1.57,
    valinit=psi,
)

# The function to be called anytime a slider's value changes
def update(val):
    phi = phi_slider.val
    theta = theta_slider.val
    psi = psi_slider.val
    results = rotate(phi, theta, psi)
    Z = results[0]
    XY = results[1]
    XZ = results[2]
    YZ = results[3]
    facesXY = results[4]
    facesXZ = results[5]
    facesYZ = results[6]
    faces = results[7]
    plot(Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces)


# register the update function with each slider
phi_slider.on_changed(update)
theta_slider.on_changed(update)
psi_slider.on_changed(update)
