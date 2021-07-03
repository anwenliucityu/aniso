'''wrap all atoms into one period by pbc'''

import numpy as np

def pbc_wrap_orthogonal(atom_coor, box_boundary):
    '''
    the box is orthogonal
    '''

    # size of box (Angstrom)
    box_size = [box_boundary[0][1] - box_boundary[0][0],
                box_boundary[1][1] - box_boundary[1][0],
                box_boundary[2][1] - box_boundary[2][0]]

    for axis in range(3):
        tempt = atom_coor[:, axis]
        tempt[np.where(tempt > box_boundary[axis][1])] -= box_size[axis]
        tempt[np.where(tempt < box_boundary[axis][0])] += box_size[axis]
        atom_coor[:, axis] = tempt

    return atom_coor
