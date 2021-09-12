# Graphical-Projections
Orthographic and isometric views of a 3D cube that rotates. Made using Matplotlib and linear algebra.

The whole Python code is found in 'PythonVersion.py'. When you run it, you can use three sliders to adjust the values of:
 - phi (angle of rotation around the X axis in degrees)
 - theta (angle of rotation around the Y axis in degrees)
 - psi (angle of rotation around the Z axis in degrees)
Then you will see the result in 3D isometric view and its orthographic projection onto the XY, XZ and YZ planes. 

There is also a Matlab version of the same program ('MatlabVersion.mlx'), although it is more outdated. It only generates the 3D view and the 2D view from the top.

How it works:
The code takes the coordinates of all the vertices of the cube and applies a linear transformation for each rotation (roll, pitch and yaw) to get a new set of coordinates which it plots in 3D. Then, these are multiplied again separately by three matrices to get the three sets of coordinates in 2D, which correspond to the top, front and side orthographic views. Then, each one of these is plotted in a different graph, all within the same figure. You can find more details on the steps followed for this in the comments of the Python code.
