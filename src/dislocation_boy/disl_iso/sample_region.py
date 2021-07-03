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
    print('The maximum cylindrical radius allowed = ', max_radius_allowed)

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
    
def upper_lower_z(atom_coor, box_boundary, distance, side):
    '''
    select region that higher than Zhi and lower than Zlo 
    distance : unit: Angstrom
    '''
    num_atom = atom_coor.shape[0]
    atom_index_selected = []
    zlo = box_boundary[2][0]
    zhi = box_boundary[2][1]
    z_select_lo = zlo + distance
    z_select_hi = zhi - distance
    if side in ('up_down','ud'):
        for i in range(num_atom):
            try_z = atom_coor[i,2]
            if try_z > z_select_hi or try_z < z_select_lo :
                atom_index_selected.append(i)
    if side in ('in','i'):
        for i in range(num_atom):
            try_z = atom_coor[i,2]
            if try_z < z_select_hi and try_z > z_select_lo :
                atom_index_selected.append(i)
    
    return atom_index_selected

def crack_tip_region(atom_coor, x_range, y_range, z_range):
    '''
    add a rectangle crack tip
    '''
    num_atom = atom_coor.shape[0]
    atom_index_selected = []
    for i in range(num_atom):
        try_x, try_y, try_z = atom_coor[i, 0:3]
        if try_x > x_range[0] and try_x < x_range[1] and try_y > y_range[0] and try_y < y_range[1] and try_z > z_range[0] and try_z < z_range[1]:
            atom_index_selected.append(i)
            
    return atom_index_selected

def crack_tip_angle_region(atom_coor, degree, start_point, catalog):
    '''
    crack size looks in x-z plane and periodic in y dierction, SHAPE is triangle
    '''
    num_atom = atom_coor.shape[0]
    atom_index_selected = []
    for i in range(num_atom):
        try_y = atom_coor[i, 1] - start_point[0]
        try_z = atom_coor[i, 2] - start_point[1]
        if try_z/try_y < np.tan(degree[0]) and try_z/try_y > np.tan(degree[1]) and try_y >0 and catalog == 'up':
            atom_index_selected.append(i)
        if try_z/try_y > np.tan(-degree[0]) and try_z/try_y < np.tan(-degree[1]) and try_y >0 and catalog == 'down':
            atom_index_selected.append(i)

    return atom_index_selected   
    
def spherical_crack(atom_coor, center, length):
    num_atom = atom_coor.shape[0]
    atom_index_selected = []
    for i in range(num_atom):
        try_x = atom_coor[i, 0] - center[0]
        try_z = atom_coor[i, 2] - center[1]
        if try_x**2 + try_z**2 < length**2:
            atom_index_selected.append(i)
    return atom_index_selected   
    '''  [example]
    center = [box_boundary[0][1]-10, 0]
    length = 8
    atom_index_selected = sample_region.spherical_crack(atom_coor, center, length)
    atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    atom_type = np.delete(atom_type, atom_index_selected)
    '''
                
                
def compress(atom_coor, compress_factor):
   num_atom = atom_coor.shape[0]
   for i in range(num_atom):
       atom_coor[i,0] = atom_coor[i, 0] * compress_factor[0]
       atom_coor[i,1] = atom_coor[i, 1] * compress_factor[1]
       atom_coor[i,2] = atom_coor[i, 2] * compress_factor[2]
   return atom_coor

def cone_crack(atom_coor, center, direction):
    num_atom = atom_coor.shape[0]
    x0 = center[0]
    y0 = center[1]
    z0 = center[2]
    atom_index_selected = []
    if direction == "left": #the crack tip direction is pointed to left
        for i in range(num_atom):
            try_x, try_y, try_z = atom_coor[i, 0:3]
            if try_y > y0 and (try_x-x0)**2+(try_z-z0)**2 < 0.4*(try_y-y0)**2:
                atom_index_selected.append(i)
    return atom_index_selected
            
