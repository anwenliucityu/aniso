'''initialize a perfect crystal prepared for introduction of dislocation'''

import numpy as np
from euler_rotation import Rx, Ry, Rz

def perfect_cryst_constructor(config, config_style, element_struct, dislocation_type, latt_const,
                              unit_cell_size, quadru=False, atom_shift=1./4):
    '''
    construct a unit cell, then replicate it in 3 directions
    available: fcc, hcp_basal_a and hcp_prism1_a
    '''
    # fcc
    # edge length of a unit cell (Angstrom)
    cell = [config["cell_x"] * latt_const,
            config["cell_y"] * latt_const,
            config["cell_z"] * latt_const]
    # unit cell basis vectors (row)
    box = [[cell[0], 0,       0      ],
           [0,       cell[1], 0      ],
           [0,       0,       cell[2]]]
            
    # shift atoms to box center
    basis_atom = config["basis_atoms"]
    basis_atom = np.array(basis_atom) + atom_shift + 0.000001
    
    # coordinates of atoms in one unit cell
    atom_coor_cell = np.matmul(basis_atom, box)

    if element_struct == 'bcc' and dislocation_type == 'screw':
        repeat_para = np.array([[cell[0],0,0],
                                [0,cell[1],-cell[2]/3],
                                [0, 0, cell[2]]])
    else:
        repeat_para = None
    
    # rotate for quadru
    if quadru == True:
        qua_rot = np.array([[0,0,1.],[1.,0,0],[0,1.,0]])
        atom_coor_cell = np.inner(atom_coor_cell, qua_rot)
        repeat_para    = np.array([repeat_para[2], repeat_para[0], repeat_para[1]])
        repeat_para   = np.inner(repeat_para, qua_rot) 

    # number of atoms in one unit cell
    num_atom_cell = atom_coor_cell.shape[0]
    # box boundaries (Angstrom)
    if config_style == 'cylinder':
        box_num_cell = [unit_cell_size[0],
                        unit_cell_size[1],
                        1]
    else:
        box_num_cell = [unit_cell_size[0],
                        unit_cell_size[1],
                        unit_cell_size[2]]
    box_num_cell_x_lo = -(box_num_cell[0] // 2)
    box_num_cell_x_hi = box_num_cell[0] + box_num_cell_x_lo
    box_num_cell_y_lo = -(box_num_cell[1] // 2)
    box_num_cell_y_hi = box_num_cell[1] + box_num_cell_y_lo
    box_num_cell_z_lo = -(box_num_cell[2] // 2)
    box_num_cell_z_hi = box_num_cell[2] + box_num_cell_z_lo
    # xlo, xhi, ylo, yhi, zlo, zhi (Angstrom)
    box_boundary = [[box_num_cell_x_lo * cell[0], box_num_cell_x_hi * cell[0] ],
                    [box_num_cell_y_lo * cell[1], box_num_cell_y_hi * cell[1] ],
                    [box_num_cell_z_lo * cell[2], box_num_cell_z_hi * cell[2] ]]
    if element_struct == 'bcc' and dislocation_type == 'screw':
        box_boundary = [[box_num_cell_x_lo * repeat_para[0][0], box_num_cell_x_hi * repeat_para[0][0] ],
                    [box_num_cell_y_lo * repeat_para[1][1], box_num_cell_y_hi * repeat_para[1][1] ],
                    [box_num_cell_z_lo * repeat_para[2][2], box_num_cell_z_hi * repeat_para[2][2] ]]
    # total number of atoms
    num_atom = box_num_cell[0] * box_num_cell[1] * box_num_cell[2] * num_atom_cell
    # replicate unit cell in three directions
    atom_coor = np.empty(shape=(num_atom, 3))
    atom_type = np.empty(num_atom)
    num_replicate = 0
    for k in range(box_num_cell_z_lo, box_num_cell_z_hi):
        for j in range(box_num_cell_y_lo, box_num_cell_y_hi):
            for i in range(box_num_cell_x_lo, box_num_cell_x_hi):
            # shift all atoms in unit cell by the same displacement
                shift = np.array([[i * cell[0], j * cell[1],  k * cell[2]]
                                 * num_atom_cell]).reshape(num_atom_cell,3)
                if element_struct == 'bcc' and dislocation_type == 'screw':
                    repeat_num = np.array([i, j, k])
            # shift all atoms in unit cell by the same displacement
                    shift = np.array([[np.dot(repeat_para[:,0], repeat_num), 
                                       np.dot(repeat_para[:,1], repeat_num), 
                                       np.dot(repeat_para[:,2], repeat_num)]
                                  * num_atom_cell]).reshape(num_atom_cell,3)
                atom_coor[num_replicate : num_replicate + num_atom_cell, :] = atom_coor_cell + shift
                atom_type[num_replicate : num_replicate + num_atom_cell] = config["type_basis_atoms"]
                num_replicate = num_replicate + num_atom_cell

    return (atom_coor, atom_type, box_boundary, repeat_para)

def perfect_hcp_constructor(config, config_style, latt_const, pot_element, dislocation_type, unit_cell_size, 
                            quadru=False, atom_shift=3./8):
    a = latt_const["a"]
    c = latt_const["c"]
    case_name = pot_element + '_' + dislocation_type
    basis_atom = config["basis_atoms"]
    
    # shift atoms to box center
    basis_atom = np.array(basis_atom) + atom_shift + 0.000001
    
    cell = [config["cell_x"] * latt_const[config["cell_x_latt_const"]],
            config["cell_y"] * latt_const[config["cell_y_latt_const"]],
            config["cell_z"] * latt_const[config["cell_z_latt_const"]]]
    box = [[cell[0], 0,       0      ],
           [0,       cell[1], 0      ],
           [0,       0,       cell[2]]]

    if "screw_a_basal" in case_name:  #CORRECT
        cell_rotation_matrix = np.eye(3)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[a*np.sqrt(3), 0, 0],
                                [0, c, 0],
                                [0, 0, a]])

    if "edge_a_basal" in case_name:  #CORRECT
        cell_rotation_matrix = np.eye(3)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[a, 0, 0],
                                [0, c, 0],
                                [0, 0, a*np.sqrt(3)]])

    if "screw_a_prismI" in case_name: #Strange
        cell_rotation_matrix = np.eye(3)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[c, 0, 0],
                                [0, a*np.sqrt(3), 0],
                                [0, 0, a]])

    if "edge_a_prismI" in case_name: #CORRECT
        cell_rotation_matrix = np.eye(3)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[a, 0, 0],
                                [0, c, 0],
                                [0, 0, a*np.sqrt(3)]])
    
    if "edge_a_pyrI" in case_name: #CORR
        cell_rotation_angle  = np.arctan2(2*c,a*np.sqrt(3))
        cell_rotation_matrix = Rx(-cell_rotation_angle)
        ela_const_rot_mat = Rx(cell_rotation_angle)
        repeat_para = np.array([[a, 0, 0],
        	                [0, np.sqrt(3)*a*np.sin(cell_rotation_angle), np.sqrt(3)*a*np.cos(cell_rotation_angle)],
        	                [0, 0, np.sqrt(3)*a/np.cos(cell_rotation_angle)]])
        	                
    if "screw_ca_pyrI" in case_name:
        cell_rotation_angle = np.arctan2(c,a)
        cell_rotation_angle_z = np.arctan2(1/2*a*np.sin(cell_rotation_angle), np.sqrt(3)/2*a)
        cell_rotation_matrix = np.dot(Rz(cell_rotation_angle_z),Rx(-cell_rotation_angle)) 
        ela_const_rot_mat = np.dot(Rz(-cell_rotation_angle_z),Rx(cell_rotation_angle))
        angle = np.arctan2(1/2*a, np.sqrt((np.sqrt(3)/2*a)**2 + c**2))
        ang = np.arctan2(c, np.sqrt(3)/2*a)
        repeat_para = np.array([[np.sqrt((1/2*a*np.sin(cell_rotation_angle))**2 + (np.sqrt(3)/2*a)**2), 0, 1/2*a*np.cos(cell_rotation_angle)],
        	                [c*np.sin(ang)*np.sin(angle), c*np.cos(ang), a/np.cos(cell_rotation_angle)-c*np.sin(ang)*np.cos(angle)],
        	                [0, 0, a/np.cos(cell_rotation_angle)]])
        	               
    if "screw_ca_pyrII" in case_name: # CORRECT
        cell_rotation_angle = np.arctan2(c,a)
        cell_rotation_matrix = Rx(-cell_rotation_angle) 
        ela_const_rot_mat = Rx((cell_rotation_angle))
        repeat_para = np.array([[a*np.sqrt(3), 0, 0],
        	                [0, a*np.sin(cell_rotation_angle), a*np.cos(cell_rotation_angle)],
        	                [0, 0, a/np.cos(cell_rotation_angle)]])
        
        	                
    if "edge_ca_pyrII" in case_name: #CORRECT
        cell_rotation_angle = np.arctan2(c,a)
        cell_rotation_matrix = Rz(-cell_rotation_angle)
        ela_const_rot_mat = Ry(-cell_rotation_angle)
        repeat_para = np.array([[a/np.cos(cell_rotation_angle), 0, 0],
        	                [-a*np.cos(cell_rotation_angle), a*np.sin(cell_rotation_angle), 0],
        	                [0, 0, a*np.sqrt(3)]])
    
    if "mixed_ca_prismI" in case_name: #CORRECT
        cell_rotation_matrix = np.eye(3)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[a, 0, 0],
        	                [0, a*np.sqrt(3), 0],
        	                [0, 0, c]])  
    
    if "mixed_ca_pyrI" in case_name:
        cell_rotation_angle = np.arctan2(2*c, np.sqrt(3)*a)
        cell_rotation_matrix = Rz(cell_rotation_angle)
        ela_const_rot_mat = Ry(cell_rotation_angle)
        repeat_para = np.array([[np.sqrt(3)*a/np.cos(cell_rotation_angle), 0, 0],
        	                [np.sqrt(3)*a*np.cos(cell_rotation_angle), np.sqrt(3)*a*np.sin(cell_rotation_angle), 0],
        	                [0, 0, a]]) 
        
    if "screw_c_prismI" in case_name:  #CORRECT  
        cell_rotation_angle = 0
        cell_rotation_matrix = Ry(cell_rotation_angle)   
        ela_const_rot_mat = np.eye(3)     
        repeat_para = np.array([[a, 0, 0],
        	                [0, a*np.sqrt(3), 0],
        	                [0, 0, c]])  
    
    if "edge_c_prismI" in case_name:   #CORRECT 
        cell_rotation_angle = 0
        cell_rotation_matrix = Ry(cell_rotation_angle)
        ela_const_rot_mat = np.eye(3)
        repeat_para = np.array([[c, 0, 0],
        	                [0, a*np.sqrt(3), 0],
        	                [0, 0, a]]) 
        	                
    if "edge_c_prismII" in case_name:    # CORRECT
        cell_rotation_angle = 0
        cell_rotation_matrix = Ry(cell_rotation_angle) 
        ela_const_rot_mat = np.eye(3)       
        repeat_para = np.array([[c, 0, 0],
        	                [0, a, 0],
        	                [0, 0, a*np.sqrt(3)]]) 
    
        	                        
    # lattice_const
    atom_coor_cell = np.matmul(basis_atom, box)
    num_atom_cell = atom_coor_cell.shape[0]
    frame_new = np.dot(ela_const_rot_mat, config["frame_new"])

    # rotation unit cell
    atom_coor_cell = np.inner(atom_coor_cell, cell_rotation_matrix)

    lattice_const = [repeat_para[0][0], repeat_para[1][1], repeat_para[2][2]]
    if quadru == True:
        qua_rot = np.array([[0,0,1.],[1.,0,0],[0,1.,0]])
        atom_coor_cell = np.inner(atom_coor_cell, qua_rot)
        repeat_para    = np.array([repeat_para[2], repeat_para[0], repeat_para[1]])
        repeat_para   = np.inner(repeat_para, qua_rot)
        latt_c_rot    = np.array([[0,1.,0], [0,0,1.],[1.,0,0]])
        frame_new     = np.dot(frame_new, latt_c_rot)	                

    # box boundaries (Angstrom)
    if config_style == 'cylinder':
        box_num_cell = [unit_cell_size[0],
                        unit_cell_size[1],
                        1]
    else:
        box_num_cell = [unit_cell_size[0],
                        unit_cell_size[1],
                        unit_cell_size[2]]
    box_num_cell_x_lo = -(box_num_cell[0] // 2)
    box_num_cell_x_hi = box_num_cell[0] + box_num_cell_x_lo
    box_num_cell_y_lo = -(box_num_cell[1] // 2)
    box_num_cell_y_hi = box_num_cell[1] + box_num_cell_y_lo
    box_num_cell_z_lo = -(box_num_cell[2] // 2)
    box_num_cell_z_hi = box_num_cell[2] + box_num_cell_z_lo
    # xlo, xhi, ylo, yhi, zlo, zhi (Angstrom)
    box_boundary = [[box_num_cell_x_lo * repeat_para[0][0], box_num_cell_x_hi * repeat_para[0][0] ],
                    [box_num_cell_y_lo * repeat_para[1][1], box_num_cell_y_hi * repeat_para[1][1] ],
                    [box_num_cell_z_lo * repeat_para[2][2], box_num_cell_z_hi * repeat_para[2][2] ]]                     
    # total number of atoms
    num_atom = box_num_cell[0] * box_num_cell[1] * box_num_cell[2] * num_atom_cell        
    # replicate unit cell in three directions
    atom_coor = np.empty(shape=(num_atom, 3))
    atom_type = np.empty(num_atom)
    num_replicate = 0
    for k in range(box_num_cell_z_lo, box_num_cell_z_hi):
        for j in range(box_num_cell_y_lo, box_num_cell_y_hi):
            for i in range(box_num_cell_x_lo, box_num_cell_x_hi):
                repeat_num = np.array([i, j, k])
            # shift all atoms in unit cell by the same displacement
                shift = np.array([[np.dot(repeat_para[:,0], repeat_num), 
                                   np.dot(repeat_para[:,1], repeat_num), 
                                   np.dot(repeat_para[:,2], repeat_num)]
                                  * num_atom_cell]).reshape(num_atom_cell,3)
                atom_coor[num_replicate : num_replicate + num_atom_cell, :] = atom_coor_cell + shift
                atom_type[num_replicate : num_replicate + num_atom_cell] = config["type_basis_atoms"]
                num_replicate = num_replicate + num_atom_cell
    return (atom_coor, atom_type, box_boundary, lattice_const, frame_new, repeat_para)
     
    
