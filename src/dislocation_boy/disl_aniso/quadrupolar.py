import numpy as np

def quadrupolar_params(disl_center, box_boundary, repeat_para, unit_cell_size):

    qua_rot = np.array([[0,0,1.],[1.,0,0],[0,1.,0]])
    box_boundary_lo = np.array([box_boundary[0][0],box_boundary[1][0],box_boundary[2][0]])
    box_boundary_hi = np.array([box_boundary[0][1],box_boundary[1][1],box_boundary[2][1]])
    box_boundary_lo = np.inner(box_boundary_lo, np.linalg.inv(qua_rot))
    box_boundary_hi = np.inner(box_boundary_hi, np.linalg.inv(qua_rot))
    Lx = box_boundary_hi[0] - box_boundary_lo[0]
    Ly = box_boundary_hi[1] - box_boundary_lo[1]
    Lz = box_boundary_hi[2] - box_boundary_lo[2]
    start_x = box_boundary_lo[0]
    start_y = box_boundary_lo[1]
    start_z = box_boundary_lo[2]
    disl_center_1 = [disl_center[0]*Lx+start_x-0.00001,
                      disl_center[1]*Ly+start_y-0.00001]
    if disl_center[0] < 0.5:
        new_center_x = disl_center[0]+0.5
    elif disl_center[0] < 1.0:
        new_center_x = disl_center[0]-0.5
    if disl_center[1] < 0.5:
        new_center_y = disl_center[1]+0.5
    elif disl_center[1] < 1.0:
        new_center_y = disl_center[1]-0.5
    S_area = repeat_para[1,1] *repeat_para[2,2]* unit_cell_size[1] *unit_cell_size[2]
    disl_center_2 = [new_center_x*Lx+start_x-0.00001, 
                     new_center_y*Ly+start_y-0.00001]
    disl_centers = []
    image = [-1,0,1]
    for i in image:
        for j in image:
            disl_center = [disl_center_1[0] + i*Lx, disl_center_1[1] + j*Ly]
            disl_centers.append(disl_center)
            disl_center = [disl_center_2[0] + i*Lx, disl_center_2[1] + j*Ly]
            disl_centers.append(disl_center)

    d = [0, Lx/2, Ly/2]
    l = [1,0,0]
    A = np.cross(l,d)

    # find four points for calculating err plane
    plane_point = np.array([[start_x, start_y, 0.],
                            [start_x+Lx, start_y, 0.],
                            [start_x, start_y+Ly, 0.],
                            [start_x+Lx,start_y+Ly,0.]])
    return disl_centers, d, l, A, S_area, plane_point
