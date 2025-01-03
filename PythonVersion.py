## Visualize a cube that rotates in 3D, and its 2D orthographic projections. Made with linear transformations.

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.widgets import Slider

# Assign angles of rotation
phi = 0
theta = 0
psi = 0

# Define rotation matrices
def rotation_matrix(phi, theta, psi):
    x_roll_matrix = np.array([[1, 0, 0],
                   [0, np.cos(phi), np.sin(phi)],
                   [0, -np.sin(phi), np.cos(phi)]])

    y_pitch_matrix = np.array([[np.cos(theta), 0, -np.sin(theta)],
                   [0, 1, 0],
                   [np.sin(theta), 0, np.cos(theta)]])

    z_yaw_matrix = np.array([[np.cos(psi), np.sin(psi), 0],
                   [-np.sin(psi), np.cos(psi), 0],
                   [0, 0, 1]])

    return x_roll_matrix @ y_pitch_matrix @ z_yaw_matrix  # Apply rotations in order


# Create function that generates coordinates of the vertices of the cube after rotating.
def rotate(phi, theta, psi):
    # Vertex coordinates of original cube
    vertices = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                         [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])

    # Generate vertex coords of rotated cube
    Z = vertices @ rotation_matrix(phi, theta, psi)


    # Vertices in each face of the cube
    faces = [[Z[0], Z[1], Z[2], Z[3]],
             [Z[4], Z[5], Z[6], Z[7]],
             [Z[0], Z[1], Z[5], Z[4]],
             [Z[2], Z[3], Z[7], Z[6]],
             [Z[1], Z[2], Z[6], Z[5]],
             [Z[4], Z[7], Z[3], Z[0]]]

    # Projection matrices
    ort_xy = np.array([[1, 0, 0],
                    [0, 1, 0]])
    ort_xz = np.array([[1, 0, 0],
                    [0, 0, 1]])
    ort_yz = np.array([[0, 1, 0],
                    [0, 0, 1]])

    # Generate new coordinates and faces
    XY = np.zeros((8, 2))
    for i in range(8): XY[i, :] = np.dot(ort_xy, Z[i, :])
    # XY = np.delete(XY, 2, 1)

    XZ = np.zeros((8, 2))
    for i in range(8): XZ[i, :] = np.dot(ort_xz, Z[i, :])
    # XZ = np.delete(XZ, 1, 1)

    YZ = np.zeros((8, 2))
    for i in range(8): YZ[i, :] = np.dot(ort_yz, Z[i, :])
    # YZ = np.delete(YZ, 0, 1)

    facesXY = [[XY[7], XY[4], XY[5], XY[6]],
               [XY[1], XY[2], XY[3], XY[0]],
               [XY[4], XY[0], XY[3], XY[7]],
               [XY[5], XY[1], XY[2], XY[6]],
               [XY[7], XY[3], XY[2], XY[6]],
               [XY[0], XY[4], XY[5], XY[1]]]
    facesXZ = [[XZ[7], XZ[4], XZ[5], XZ[6]],
               [XZ[1], XZ[2], XZ[3], XZ[0]],
               [XZ[4], XZ[0], XZ[3], XZ[7]],
               [XZ[5], XZ[1], XZ[2], XZ[6]],
               [XZ[7], XZ[3], XZ[2], XZ[6]],
               [XZ[0], XZ[4], XZ[5], XZ[1]]]
    facesYZ = [[YZ[7], YZ[4], YZ[5], YZ[6]],
               [YZ[1], YZ[2], YZ[3], YZ[0]],
               [YZ[4], YZ[0], YZ[3], YZ[7]],
               [YZ[5], YZ[1], YZ[2], YZ[6]],
               [YZ[7], YZ[3], YZ[2], YZ[6]],
               [YZ[0], YZ[4], YZ[5], YZ[1]]]

    # Return all important values
    return Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces


# Get values from rotate()
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

ax = fig.add_subplot(221, projection='3d', proj_type="ortho", azim=-45, elev=35, title="Isometric")
ax.set_box_aspect([1, 1, 1])

ax2 = fig.add_subplot(222, title="Top")
ax2.set_aspect('equal', adjustable='box')
ax3 = fig.add_subplot(224, title="Front")
ax3.set_aspect('equal', adjustable='box')
ax4 = fig.add_subplot(223, title="Side")
ax4.set_aspect('equal', adjustable='box')

plt.subplots_adjust(bottom=0.3)
plt.subplots_adjust(wspace=0.5, hspace=0.5)


# Define the function for updating the plot
def plot(Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces):
    # Clear previous points
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()

    ax.set_box_aspect([1, 1, 1])

    ax.set_title("Isometric")

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

    ax2.set_title("Top")
    ax3.set_title("Front")
    ax4.set_title("Side")

    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')
    ax4.set_xlabel('Y')
    ax4.set_ylabel('Z')

    # Plot new graphs
    ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2])

    ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=1, edgecolors='blue', alpha=.25))

    ax2.scatter(XY[:, 0], XY[:, 1])
    ax3.scatter(XZ[:, 0], XZ[:, 1])
    ax4.scatter(YZ[:, 0], YZ[:, 1])

    for i in range(6):
        poly = Polygon(facesXY[i], facecolor='cyan', edgecolor='blue', alpha=0.25)
        ax2.add_patch(poly)
    for i in range(6):
        poly = Polygon(facesXZ[i], facecolor='cyan', edgecolor='blue', alpha=0.25)
        ax3.add_patch(poly)
    for i in range(6):
        poly = Polygon(facesYZ[i], facecolor='cyan', edgecolor='blue', alpha=0.25)
        ax4.add_patch(poly)


# Make the first plot
plot(Z, XY, XZ, YZ, facesXY, facesXZ, facesYZ, faces)

# Create sliders
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
    plt.draw()
    fig.canvas.draw_idle()  #Ensures the figure updates properly


# Register the update function with each slider
phi_slider.on_changed(update)
theta_slider.on_changed(update)
psi_slider.on_changed(update)

plt.show()
