from voigt_notation import voigt_one_to_two, voigt_two_to_one, \
                           four_tensor_one_to_two, four_tensor_two_to_one,\
                           four_tensor_matrix_to_list,\
                           four_tensor_list_to_matrix
import numpy as np
from sympy import Eijk
import itertools
from integrate import zero_theta_integrate

def initialize_disl_config(config, config_style, latt_const, element_struct, frame_new, lattice_const, disl_center, 
                           repeat_para=None, ):
    '''
    pre-processing according to settings of dislocaiton configuration
    output: position of disloction (Angstrom),
    e'1, e'2, e'3 components of Burgers vector (Angstrom)
    '''
    disl_line_direction = config["disl_line_direction"]
    if element_struct == 'fcc' or element_struct == 'bcc':
        burgers = [i * latt_const for i in config["burgers"]]
        disl_center = [disl_center[0] * config["cell_x"] * latt_const,
                   disl_center[1] * config["cell_y"] * latt_const]
    if element_struct == 'hcp': 
        burgers = np.multiply(config["burgers"], lattice_const)
        disl_center = [repeat_para.sum(axis=0)[0]*disl_center[0], repeat_para.sum(axis=0)[1]*disl_center[1]]
    b_vector = []
    for i in range(3):
        b = np.dot(burgers , np.transpose(frame_new)[i])
        b_vector.append(b)
    return disl_center, b_vector


def stress_field (new_elastic_constant, atom_coor, disl_center, S, B, s_theta, q_theta, b_vector):
    num_atom = atom_coor.shape[0]
    x_coor = atom_coor[:, 0] - disl_center[0]
    y_coor = atom_coor[:, 1] - disl_center[1]
    theta  = np.arctan2(y_coor, x_coor).reshape((num_atom,1))
    R  = np.sqrt(x_coor**2+y_coor**2).reshape((num_atom,1))
    iec      = new_elastic_constant
    sigma    = np.zeros([num_atom,9])
    for i,j,k,l,s in itertools.product(range(3),repeat=5):
        ij      = four_tensor_matrix_to_list[(i, j)]
        ks      = four_tensor_matrix_to_list[(k, s)]
        ij__    = voigt_two_to_one[(i, j)]
        __kl    = voigt_two_to_one[(k, l)]
        C_ijkl  = four_tensor_two_to_one[(ij__, __kl)]
    if l == 0:
        sigma[:,ij] = sigma[:,ij] + 1/(2*np.pi*R[:,0])*iec[C_ijkl]*b_vector[s]*(-np.cos(theta[:,0])*S[:,ks])
        for r in range(3):
            kr = four_tensor_matrix_to_list[(k, r)]
            rs = four_tensor_matrix_to_list[(r, s)]
            sigma[:,ij] = sigma[:,ij] + 1/(2*np.pi*R[:,0])*iec[C_ijkl]*b_vector[s]*\
                     (-np.sin(theta[:,0])*(s_theta[:,kr]*S[:,rs]+q_theta[:,kr]*B[:,rs]))
          
    if l == 1:
        sigma[:,ij] = sigma[:,ij] + 1/(2*np.pi*R[:,0])*iec[C_ijkl]*b_vector[s]*(-np.sin(theta[:,0])*S[:,ks])
        for r in range(3):
            kr = four_tensor_matrix_to_list[(k, r)]
            rs = four_tensor_matrix_to_list[(r, s)]
            sigma[:,ij] = sigma[:,ij] + 1/(2*np.pi*R[:,0])*iec[C_ijkl]*b_vector[s]*\
                     (np.cos(theta[:,0])*(s_theta[:,kr]*S[:,rs]+q_theta[:,kr]*B[:,rs]))
    return sigma


def displacement_field(atom_coor, disl_center, S, B, S_theta, Q_theta, b_vector):
    atom_num = atom_coor.shape[0]
    x_coor = atom_coor[:, 0] - disl_center[0]
    y_coor = atom_coor[:, 1] - disl_center[1]
    r = np.sqrt(x_coor**2+y_coor**2).reshape((atom_num,1))
    u        = np.zeros([atom_num,3])
    for k, j in itertools.product(range(3), repeat = 2):
        kj = four_tensor_matrix_to_list[(k, j)]
        u[:,k] = u[:,k] + 1/(2*np.pi) * (-S[0,kj]*np.log(r[:,0]))*b_vector[j]
        for i in range(3):
            ki = four_tensor_matrix_to_list[(k, i)]
            ij = four_tensor_matrix_to_list[(i, j)]
            u[:,k] = u[:,k] + 1/(2*np.pi) * (S_theta[:,ki] * S[0,ij] + Q_theta[:,ki] * B[0,ij])*b_vector[j]

    return u

def partial_init(disl_center, partial_split_dis, b_vector, dislocation_type):
    if dislocation_type == 'screw_ca_pyrII':
        b_vectors = [np.array(b_vector)/2, np.array(b_vector)/2]
        disl_centers = [[disl_center[0]+partial_split_dis/2, disl_center[1]], 
                        [disl_center[0]-partial_split_dis/2, disl_center[1]]]
    return b_vectors, disl_centers
        


# define a function for multiprocessing to calculate displacement for different dislocation line postions
def displacement_fields(i, s, q, S, B, atom_coor_ref, disl_center, S_list, Q_list, b_vector):
    # initialize b for pair dislocations
    if i%2 == 1:
        b_vector = -np.array(b_vector)

    # main func
    S_theta, s_theta = zero_theta_integrate(s, atom_coor_ref, disl_center, S_list)
    Q_theta, q_theta = zero_theta_integrate(q, atom_coor_ref, disl_center, Q_list)
    u_displacement_field = displacement_field(atom_coor_ref, disl_center, S, B, S_theta, Q_theta, b_vector)
    return u_displacement_field

