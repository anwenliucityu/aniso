'''select atoms in a region of sample (return indices of the selected atoms)'''

import numpy as np

def cylinder_z(atom_coor, center, radius, side):
    '''
    cylindrical region with its axis in z idrection
    atom_coor: an array of size num_atom x 3
    side: 'in' select atoms inside region; 'out' select atoms outside region
    '''

    max_radius_allowed = min(np.abs(center[0] - np.min(atom_coor[:,0])),
                             np.abs(center[0] - np.max(atom_coor[:,0])),
                             np.abs(center[1] - np.min(atom_coor[:,1])),
                             np.abs(center[1] - np.max(atom_coor[:,1])))

    radius_sq = radius**2
    num_atom = atom_coor.shape[0]
    atom_index_selected = []
    if side in ('out', 'o'):
        for i in range(num_atom):
            try_x, try_y = atom_coor[i, 0:2] - center
            if try_x**2 + try_y**2 > radius_sq:
                atom_index_selected.append(i)
    elif side in ('in', 'i'):
        for i in range(num_atom):
            try_x, try_y = atom_coor[i, 0:2] - center
            if try_x**2 + try_y**2 < radius_sq:
                atom_index_selected.append(i)

    return atom_index_selected
