'''input file'''

import math

Al_edge_config = {
    "case_name"         : "Al_edge",
    "elem"              : "EAM_Mishin_Al",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 70,#44,
    "box_num_cell_y"    : 40,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(6.), -1./math.sqrt(6.), 2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Al_screw_config = {
    "case_name"         : "Al_screw",
    "elem"              : "EAM_Mishin_Al",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 40,#52,
    "box_num_cell_y"    : 27,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.5,
    "disl_center_y"     : 0.5,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
}

Al_mix_config = {
    "case_name"         : "Al_mix",
    "elem"              : "EAM_Mishin_Al",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 88,#44,
    "box_num_cell_y"    : 36,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Cu_edge_config = {
    "case_name"         : "Cu_edge",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 140,#44,
    "box_num_cell_y"    : 58,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : -1./4,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Cu_screw_config = {
    "case_name"         : "Cu_screw",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM1_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 103,#52,
    "box_num_cell_y"    : 72,#36,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(2.), -1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
}

Cu_mix_config = {
    "case_name"         : "Cu_mix",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM1_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 176,#44,
    "box_num_cell_y"    : 80,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Ni_mix_config = {
    "case_name"         : "Ni_mix",
    "elem"              : "EAM_Mishin_Ni",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Ni", 2 : "EAM_Mishin_Ni"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 132,#44,
    "box_num_cell_y"    : 54,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}
Al_screw_300K_config = {
    "case_name"         : "Al_screw",
    "elem"              : "EAM_Mishin_Al_300K",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 40,#52,
    "box_num_cell_y"    : 30,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.5,
    "disl_center_y"     : 0.5,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
                           }
W_screw_300K_config = {
    "case_name"         : "W_screw",
    "elem"              : "EAM_Ma_W_300K",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : 1,
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [0.5  ,  0.5,    1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Ma_W", 2 : "EAM_Ma_W"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 55,#52,
    "box_num_cell_y"    : 55,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0,0,1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.25,
    "disl_center_y"     : 0.25,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],

                           }
Fe_screw_300K_config = {
    "case_name"         : "Fe_screw",
    "elem"              : "EAM_Mendelev_Fe_300K",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : 1,
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2  ,  1./2,    1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mendelev_Fe", 2 : "EAM_Mendelev_Fe"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 70,#52,
    "box_num_cell_y"    : 70,#36,
    "box_num_cell_z"    : 200,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0,0,1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.25,
    "disl_center_y"     : 0.25,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],

                           }

Nb_screw_config = {
    "case_name"         : "Nb_screw",
    "elem"              : "Zhang_Nb_EAM",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : 1,
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [0.5  ,  0.5,    1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "Zhang_Nb_EAM", 2 : "Zhang_Nb_EAM"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 55,#52,
    "box_num_cell_y"    : 55,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0,0,1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.25,
    "disl_center_y"     : 0.25,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],}

Ta_screw_config = {
    "case_name"         : "Ta_screw",
    "elem"              : "Ravelo_Ta_EAM",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : 1,
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [0.5  ,  0.5,    1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "Ta_EAM", 2 : "Ta_EAM"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 55,#52,
    "box_num_cell_y"    : 55,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0,0,1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.25,
    "disl_center_y"     : 0.25,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],}
Pt_screw_config = {
    "case_name"         : "Pt_screw",
    "elem"              : "Zhou_Pt_EAM",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 40,#52,
    "box_num_cell_y"    : 30,#36,
    "box_num_cell_z"    : 150,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.5,
    "disl_center_y"     : 0.5,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
                           }

