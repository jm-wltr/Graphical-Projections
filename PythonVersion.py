#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:55:44 2021
@author: Jaimew
"""

## Visualize a cube that rotates in 3D, and its 2D projection. Made with linear transformations. 
## This file only generates the coordinates of the vertices of the rotated cube.

import numpy as np
import math

# Coordinates of vertices
points = np.mat([[0.5, 0.5, -0.5], \
                 [0.5, 0.5, 0.5], \
                 [0.5, -0.5, -0.5], \
                 [0.5, -0.5,  0.5], \
                 [-0.5, 0.5, -0.5], \
                 [-0.5, 0.5, 0.5], \
                 [-0.5, -0.5, -0.5], \
                 [-0.5, -0.5, 0.5]])
                       
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

for n in range(0,8):
    points[n] = np.transpose(XrollMatrix * np.transpose(points[n]))
    points[n] = np.transpose(YpitchMatrix * np.transpose(points[n]))
    points[n] = np.transpose(ZyawMatrix * np.transpose(points[n]))
