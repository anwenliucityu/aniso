'''wrap all atoms into one period by pbc'''

import numpy as np

def pbc_wrap_orthogonal(atom_coor, box_boundary,):
    '''
    the box is orthogonal
    '''

    # size of box (Angstrom)
    box_size = [box_boundary[0][1] - box_boundary[0][0],
                box_boundary[1][1] - box_boundary[1][0],
                box_boundary[2][1] - box_boundary[2][0]]
    
    for axis in range(3):
    	for i in range(1,50):
            tempt = atom_coor[:, axis]
            tempt[np.where(tempt > box_boundary[axis][1])] -= box_size[axis]
            tempt[np.where(tempt < box_boundary[axis][0])] += box_size[axis]
            atom_coor[:, axis] = tempt

    return atom_coor

def pbc_wrap_tilt(atom_coor, box_boundary,tilt):
    # size of box (Angstrom)
    xlen = box_boundary[0][1] - box_boundary[0][0]
    xlo = box_boundary[0][0]
    xhi = box_boundary[0][1]
    ylen = box_boundary[1][1] - box_boundary[1][0]
    ylo = box_boundary[1][0]
    yhi = box_boundary[1][1]
    zlen = box_boundary[2][1] - box_boundary[2][0]
    zlo = box_boundary[2][0]
    zhi = box_boundary[2][1]
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2] 
    xaxis = atom_coor[:,0]
    yaxis = atom_coor[:,1]
    zaxis = atom_coor[:,2]
    p1 = [xlo, ylo, zlo]
    p2 = [xlo+xy, yhi, zlo]
    p3 = [xlo+xz, ylo, zhi]
    p4 = [xhi, ylo, zlo]
    p5 = [xhi+xy, yhi, zlo]
    p6 = [xhi+xz, ylo, zhi]
    p7 = [xlo+xy+xz, yhi, zhi]
    p8 = [xhi+xy+xz, yhi, zhi]

    # x axis
    # right bound: point(xhi, zlo,) k=zlen/xz
    # left bound : point(xlo, zlo), k = zlen/xz
    for i in range(1,5):
        # x axis
        a1,b1,c1,d1 = get_panel(p1,p2,p3)
        a2,b2,c2,d2 = get_panel(p4,p5,p6)
        x_left = np.where(a1*atom_coor[:,0]+b1*atom_coor[:,1]+c1*atom_coor[:,2]+d1<0)
        x_right = np.where(a2*atom_coor[:,0]+b2*atom_coor[:,1]+c2*atom_coor[:,2]+d2>0)

        atom_coor[:,0][x_left] =  atom_coor[:,0][x_left] + xlen
        atom_coor[:,0][x_right] = atom_coor[:,0][x_right] - xlen


    return atom_coor


def get_panel(p1,p2,p3):
    a = (p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1])
    b = (p2[2]-p1[2])*(p3[0]-p1[0])-(p2[0]-p1[0])*(p3[2]-p1[2])
    c = (p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0])
    d = -(a*p1[0]+b*p1[1]+c*p1[2])
    return a,b,c,d

def pbc_wrap_tilt_box(box_boundary, tilt):
    wrap = 1/2
    xy = tilt[0]*wrap
    xz = tilt[1]*wrap
    yz = tilt[2]*wrap
    box_boundary = [[box_boundary[0][0] -xy-xz, box_boundary[0][1]-xy-xz],
                    [box_boundary[1][0] -yz, box_boundary[1][1]-yz],
                    [box_boundary[2][0] , box_boundary[2][1]]]
    return box_boundary
    
def duplicate_z(atom_coor, atom_type, box_boundary, duplicate_time):
    # duplicate box boundary in z direction
    num_atom = atom_coor.shape[0]
    zlo = box_boundary[2][0]
    zhi = box_boundary[2][1]
    repeat_unit = zhi - zlo
    new_zhi = zlo + repeat_unit * duplicate_time
    new_box_boundary = [[box_boundary[0][0], box_boundary[0][1]],
                        [box_boundary[1][0], box_boundary[1][1]],
                        [zlo, new_zhi]]
    
    new_coor = []
    new_type = []
    shift = np.array([0,0,repeat_unit])
    for i in range(duplicate_time):
        new_coor = np.append(new_coor, np.array(atom_coor) + shift*i)
        new_type = np.append(new_type, atom_type)
    return new_coor.reshape((num_atom*duplicate_time,3)), new_type.reshape((num_atom*duplicate_time,1)), new_box_boundary
