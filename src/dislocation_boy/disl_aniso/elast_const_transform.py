'''
transform of elastic constant C_ijkl
due to transform of coordinate system
'''

import itertools
import numpy as np
from voigt_notation import voigt_one_to_two, voigt_two_to_one, \
                           four_tensor_one_to_two, four_tensor_two_to_one

def elast_const_transform(elast_const_initial, frame_initial, frame_new):
    '''
    transform of elastic constant C_ijkl
    representation of 4th-rank tensor: C_ijkl -> C_IJ -> C_A (21)
    '''
    # rotation matrix Q_ij
    rotation_mat = np.dot(np.transpose(frame_new), frame_initial)
    elast_const_new = np.zeros(21)
    for c_new_idx in four_tensor_one_to_two:
        idx_i, idx_j = voigt_one_to_two[four_tensor_one_to_two[c_new_idx][0]]
        idx_k, idx_l = voigt_one_to_two[four_tensor_one_to_two[c_new_idx][1]]
        for idx_g, idx_h, idx_m, idx_n in itertools.product(range(3), repeat=4):
            c_initial_idx = four_tensor_two_to_one[(voigt_two_to_one[(idx_g, idx_h)],
                                                    voigt_two_to_one[(idx_m, idx_n)])]
            elast_const_new[c_new_idx] += rotation_mat[idx_i, idx_g] *\
                                          rotation_mat[idx_j, idx_h] *\
                                          rotation_mat[idx_k, idx_m] *\
                                          rotation_mat[idx_l, idx_n] *\
                                          elast_const_initial[c_initial_idx]
    return elast_const_new

