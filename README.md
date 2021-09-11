# Graphical-Projections
Project a 3D cube to a 2D plane and render an image using Matplotlib graphs.

The Matlab version of this project (MatlabVersion.mlx) is working perfectly although I'm not improving it anymore. When you run it, you can use the sliders to change the value of
 - phi (angle of rotation around the X axis in degrees)
 - theta (angle of rotation around the Y axis in degrees)
 - psi (angle of rotation around the Z axis in degrees)
Then you will see the result in 3D and its orthographic projection onto the XY plane. 

The Python version (PythonVersion.py) is supposed to be more or less the same as MatlabVersion.mlx, but in Python. For now, it has all the same features except for the sliders to change the value of phi, theta and psi (in radians). Therefore, you have to change these variables by manually editing the code. Then, you will get the same graphs as in Matlab, but with a slightly different format.

How it works:
The code takes the coordinates of all the vertices of the cube and applies a linear transformation for each rotation (roll, pitch and yaw) to get a new set of coordinates which it plots in 3D. Then, these are multiplied again by another matrix (ortXY = np.mat([[1, 0, 0], [0, 1, 0]]) to get the coordinates in 2D, which are plotted in a new graph.
