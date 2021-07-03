'''plot the configuration of a dislocation as constructed
fig must be defined in advance. plt.show() is needed after calling these functions'''

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def disl_config_plot(atom_coor, atom_type, box_boundary):
    '''
    plot atomic structure of a dislocation (before relaxation)
    '''

    plt.scatter(atom_coor[:,0], atom_coor[:,1], s=0.5, c=atom_type)
    plt.xlabel(r'$x$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$y$ ($\mathrm{\AA}$)')
    plt.xlim(box_boundary[0])
    plt.ylim(box_boundary[1])
    plt.gca().set_aspect('equal', adjustable='box')


def region_boundary_plot(fig, sample_center, sample_radius, shell_radius):
    '''
    add the boundaries of the region where atoms will be relaxed and the frozen region
    fig must be defined in advance
    '''

    axes = fig.gca()
    sample_circle = plt.Circle((sample_center[0], sample_center[1]), sample_radius,
                               linestyle='--', color='k', fill=False)
    axes.add_artist(sample_circle)
    shell_circle = plt.Circle((sample_center[0], sample_center[1]), shell_radius,
                              linestyle='--', color='k', fill=False)
    axes.add_artist(shell_circle)


def disl_position_plot(disl_center):
    '''
    label the position of dislocation in the x-y plane
    '''

    plt.plot(*disl_center, marker="x")


def disl_config3d_plot(fig, atom_coor, atom_type, box_boundary):
    '''
    plot atomic structure of a dislocation in 3d space
    fig must be defined in advance
    '''

    axes = fig.gca(projection='3d')
    axes.scatter(atom_coor[:,0], atom_coor[:,1], atom_coor[:,2], s=8, c=atom_type)
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.set_zlabel('z')
    # matplotlib does not support 'equal' aspect in 3d plot yet
    # Create cubic bounding box to simulate equal aspect ratio
    xlo, xhi = box_boundary[0]
    ylo, yhi = box_boundary[1]
    zlo, zhi = box_boundary[2]
    max_range = np.array([xhi - xlo, yhi - ylo,  zhi - zlo]).max()
    corner_x_arr = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xhi + xlo)
    corner_y_arr = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(yhi + ylo)
    corner_z_arr = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zhi + zlo)
    for corner_x, corner_y, corner_z in zip(corner_x_arr, corner_y_arr, corner_z_arr):
        axes.plot([corner_x], [corner_y], [corner_z], 'w')
