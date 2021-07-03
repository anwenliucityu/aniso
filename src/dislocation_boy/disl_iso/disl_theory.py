'''displace atoms by displacement field of a dislocation line'''

import numpy as np

def initialize_disl_config(config, elem_property):
    '''
    pre-processing according to settings of dislocaiton configuration
    output: position of disloction (Angstrom),
    screw and edge components of Burgers vector (Angstrom)
    '''

    latt_const = elem_property["latt_const"]
    burgers = [i * latt_const for i in config["burgers"]]
    disl_line_direction = config["disl_line_direction"]
    disl_center = [config["disl_center_x"] * config["cell_x"] * latt_const,
                   config["disl_center_y"] * config["cell_y"] * latt_const]
    b_screw = np.dot(burgers, disl_line_direction)
    b_edge = -np.linalg.norm(burgers - np.dot(b_screw, disl_line_direction))
    print(b_screw,b_edge)
    return (disl_center, b_screw, b_edge)


def iso_disl_displ(atom_coor, elem_property,
                  disl_center, b_screw, b_edge):
    '''
    calculate displacement on each atom due to introduction of a dislocation in z axis
    based on isotropic elasticity (see geometry in Hirth & Lothe)
    dislocation line direction is a unit vector along z axis
    if edge, extra half plane is at x = 0, y > 0
    slip plane is x-z plane
    atom_coor: initial coordinates of atoms, size is num_atom x 3
    '''

    c12 = elem_property["C12"]
    c44 = elem_property["C44"]
    poisson_ratio = c12/(2 * (c44 + c12))
    x_coor = atom_coor[:, 0] - disl_center[0]
    y_coor = atom_coor[:, 1] - disl_center[1]
    xsq = x_coor**2
    ysq = y_coor**2
    xsq_ysq_sum = xsq + ysq
    arctan_y_x = np.arctan2(y_coor, x_coor)
    for i in range(arctan_y_x.shape[0]):
        if arctan_y_x[i] < 0:
            arctan_y_x[i] += np.pi*2
    atom_displ = np.zeros(shape=(atom_coor.shape[0], 3))

    if b_edge != 0.:
        # dislplacement due to edge component
        atom_displ[:,0] += b_edge/(2*np.pi) \
                * (arctan_y_x
                   + x_coor * y_coor / (2*(1-poisson_ratio) * xsq_ysq_sum))
        atom_displ[:,1] += -b_edge/(2*np.pi) \
                * ((1-2*poisson_ratio)/(4*(1-poisson_ratio)) * np.log(xsq_ysq_sum)
                   + (xsq - ysq) / (4*(1-poisson_ratio) * xsq_ysq_sum))

    if b_screw != 0.:
        # displacement due to screw component
        atom_displ[:,2] += b_screw/(2*np.pi) * arctan_y_x

    return atom_displ

def iso_disl_displ_x(atom_coor, elem_property,
                  disl_center, b_screw, b_edge):
    '''
    calculate displacement on each atom due to introduction of a dislocation in z axis
    based on isotropic elasticity (see geometry in Hirth & Lothe)
    dislocation line direction is a unit vector along z axis
    if edge, extra half plane is at x = 0, y > 0
    slip plane is x-z plane
    atom_coor: initial coordinates of atoms, size is num_atom x 3
    '''

    c12 = elem_property["C12"]
    c44 = elem_property["C44"]
    poisson_ratio = c12/(2 * (c44 + c12))
    y_coor = atom_coor[:, 1] - disl_center[1]
    z_coor = atom_coor[:, 2] - disl_center[0]
    ysq = y_coor**2
    zsq = z_coor**2
    ysq_zsq_sum = ysq + zsq
    arctan_z_y = np.arctan2(z_coor, y_coor)
    for i in range(arctan_z_y.shape[0]):
        if arctan_z_y[i] < 0:
            arctan_z_y[i] += np.pi*2
    atom_displ = np.zeros(shape=(atom_coor.shape[0], 3))

    if b_edge != 0.:
        # dislplacement due to edge component
        atom_displ[:,0] += b_edge/(2*np.pi) \
                * (arctan_z_y
                   + y_coor * z_coor / (2*(1-poisson_ratio) * ysq_zsq_sum))
        atom_displ[:,1] += -b_edge/(2*np.pi) \
                * ((1-2*poisson_ratio)/(4*(1-poisson_ratio)) * np.log(ysq_zsq_sum)
                   + (ysq - zsq) / (4*(1-poisson_ratio) * ysq_zsq_sum))

    if b_screw != 0.:
        # displacement due to screw component
        atom_displ[:,2] += b_screw/(2*np.pi) * arctan_z_y

    return atom_displ


def iso_disl_displ_oppo_x(atom_coor, elem_property,
                  disl_center, b_screw, b_edge):
    '''
    calculate displacement on each atom due to introduction of a dislocation in z axis
    based on isotropic elasticity (see geometry in Hirth & Lothe)
    dislocation line direction is a unit vector along z axis
    if edge, extra half plane is at x = 0, y > 0
    slip plane is x-z plane
    atom_coor: initial coordinates of atoms, size is num_atom x 3
    '''

    c12 = elem_property["C12"]
    c44 = elem_property["C44"]
    poisson_ratio = c12/(2 * (c44 + c12))
    y_coor = - atom_coor[:, 1] - disl_center[1]
    z_coor =  atom_coor[:, 2] - disl_center[0]
    ysq = y_coor**2
    zsq = z_coor**2
    ysq_zsq_sum = ysq + zsq
    arctan_z_y = np.arctan2(z_coor, y_coor)
    for i in range(arctan_z_y.shape[0]):
        if arctan_z_y[i] < 0:
            arctan_z_y[i] += np.pi*2
    atom_displ = np.zeros(shape=(atom_coor.shape[0], 3))

    if b_edge != 0.:
        # dislplacement due to edge component
        atom_displ[:,0] += b_edge/(2*np.pi) \
                * (arctan_z_y
                   + y_coor * z_coor / (2*(1-poisson_ratio) * ysq_zsq_sum))
        atom_displ[:,1] += -b_edge/(2*np.pi) \
                * ((1-2*poisson_ratio)/(4*(1-poisson_ratio)) * np.log(ysq_zsq_sum)
                   + (ysq - zsq) / (4*(1-poisson_ratio) * ysq_zsq_sum))

    if b_screw != 0.:
        # displacement due to screw component
        atom_displ[:,2] += b_screw/(2*np.pi) * arctan_z_y

    return atom_displ

def iso_disl_stress(atom_coor, elem_property,
                  disl_center, b_screw, b_edge):
    '''
    calculate stress on each atom due to introduction of a dislocation in z axis
    based on isotropic elasticity (see geometry in Hirth & Lothe)
    output 6 components of stress of each atom; the order is:
    0: xx, 1: yy, 2: zz, 3: yz, 4: zx, 5: xy
    '''

    c12 = elem_property["C12"]
    c44 = elem_property["C44"]
    poisson_ratio = c12/(2 * (c44 + c12))
    shear_modulus = elem_property["shear_modulus"]   # (Pa)
    x_coor = atom_coor[:, 0] - disl_center[0]
    y_coor = atom_coor[:, 1] - disl_center[1]
    xsq = x_coor**2
    ysq = y_coor**2
    xsq_ysq_sum = xsq + ysq
    xsq_ysq_sum_sq = xsq_ysq_sum**2
    mu_b_edge_pi_nu = shear_modulus * b_edge / (2. * np.pi * (1. - poisson_ratio))
    mu_b_screw_pi = shear_modulus * b_screw / (2. * np.pi)
    atom_sigma = np.zeros(shape=(atom_coor.shape[0], 6))

    if b_edge != 0.:
        atom_sigma[:,0] += - mu_b_edge_pi_nu \
                * y_coor * (3.*xsq + ysq) / xsq_ysq_sum_sq
        atom_sigma[:,1] += mu_b_edge_pi_nu \
                * y_coor * (xsq - ysq) / xsq_ysq_sum_sq
        atom_sigma[:,2] += - 2. * poisson_ratio * mu_b_edge_pi_nu \
                * y_coor / xsq_ysq_sum
        atom_sigma[:,5] += mu_b_edge_pi_nu \
                * x_coor * (xsq - ysq) / xsq_ysq_sum_sq

    if b_screw != 0.:
        #atom_sigma[:,3] += mu_b_screw_pi * x_coor / xsq_ysq_sum
        # anisotropic stress field sigma 23
        atom_sigma[:,3] += 2.46666667e9*b_screw/(2*np.pi)*y_coor/xsq_ysq_sum + 2.91333333e10*b_screw/(2*np.pi)*x_coor/xsq_ysq_sum
        atom_sigma[:,4] += - mu_b_screw_pi * y_coor / xsq_ysq_sum

    return atom_sigma
