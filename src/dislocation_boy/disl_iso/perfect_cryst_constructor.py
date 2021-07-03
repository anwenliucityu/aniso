'''initialize a perfect crystal prepared for introduction of dislocation'''

import numpy as np

def perfect_cryst_constructor(config, elem_property):
    '''
    construct a unit cell, then replicate it in 3 directions
    '''

    # lattice constant (Angstrom)
    latt_const = elem_property["latt_const"]
    # edge length of a unit cell (Angstrom)
    cell = [config["cell_x"] * latt_const,
            config["cell_y"] * latt_const,
            config["cell_z"] * latt_const]
    # unit cell basis vectors (row)
    box = [[cell[0], 0,       0      ],
           [0,       cell[1], 0      ],
           [0,       0,       cell[2]]]
    # coordinates of atoms in one unit cell
    atom_coor_cell = np.matmul(config["basis_atoms"], box)
    # number of atoms in one unit cell
    num_atom_cell = atom_coor_cell.shape[0]
    # box boundaries (Angstrom)
    box_num_cell = [config["box_num_cell_x"],
                    config["box_num_cell_y"],
                    config["box_num_cell_z"]]
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
                shift = np.array([[i * cell[0], j * cell[1], k * cell[2]]
                                  * num_atom_cell]).reshape(num_atom_cell,3)
                atom_coor[num_replicate : num_replicate + num_atom_cell, :] \
                        = atom_coor_cell + shift
                atom_type[num_replicate : num_replicate + num_atom_cell] \
                        = config["type_basis_atoms"]
                num_replicate = num_replicate + num_atom_cell

    return (atom_coor, atom_type, box_boundary)
