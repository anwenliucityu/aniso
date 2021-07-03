import numpy as np
from sympy import Eijk
import itertools
import datetime

def zero_pi_integrate (integrand_matrix):
  '''
  --------------------------------
  zero_pi_integrate function returns the S_ij, Q_ij and B_ij which are integrated in the range [0, pi),
  and the integral result for each theta.
  -----------------------
  varible definition:
  integrand_matrix: should be the integrand function calculated s_mat, q_mat and b_mat
  '''
  timestep = integrand_matrix.shape[0]
  d_theta = np.pi/(timestep-1)
  integration_matrix = np.zeros([1,9])
  integration_theta_list = np.zeros([1,9])
  for i in range(timestep-1):
    # [0, pi) integral result, which in the form of 1*9 list
    integration_matrix[0,:] = integration_matrix[0,:] -1/np.pi * (integrand_matrix[i,:] + integrand_matrix[i+1,:]) * d_theta/2 
    # [0, theta) integral result, which in the form of of N*9 list, where N reps the number of theta in the theta_list
    integration_theta_list = np.r_[integration_theta_list, integration_matrix] 
  return integration_matrix, integration_theta_list



def zero_theta_integrate (integrand_matrix, atom_coor, disl_center, integration_theta_list):
  '''
  ---------------
  function aim:
  return the S_hat (or Q_hat) list for each theta, and s_theta (or q_theta) list for each theta
  -----------------
  varible definition:
  integrand_matrix: s_mat, q_mat, b_mat
  R_coor:                 r for each atom position with regards to dislocation line position
  theta:                  theta list
  integration_theta_list: the list calculated from zero_pi_fun (the second return value)
  '''
  atom_num = atom_coor.shape[0]
  x_coor = atom_coor[:, 0] - disl_center[0]
  y_coor = atom_coor[:, 1] - disl_center[1]
  r = np.sqrt(x_coor**2+y_coor**2).reshape((atom_num,1))
  theta = np.arctan2(y_coor, x_coor)
  theta_list = theta.reshape((atom_num,1))*1
  integrate_stepnumber = integrand_matrix.shape[0]
  d_theta = np.pi/(integrate_stepnumber-1)
  y = np.linspace(0, np.pi, integrate_stepnumber)
  captial_matrix_theta = np.zeros([atom_num,9])
  lowercase_tensor_theta = np.zeros([atom_num,9])
  
  for i in range(atom_num):
    if theta_list[i] < 0:
      theta_list[i] = theta_list[i] + np.pi
      captial_matrix_theta[i,:] = -1*np.pi*integration_theta_list[integrate_stepnumber-1,:]
      lowercase_tensor_theta[i,:] = integrand_matrix[-1,:]
      
    for j in range(integrate_stepnumber):
      if theta_list[i] < y[j]:
        remainder = theta_list[i]-y[j-1]
        lowercase_tensor_theta[i,:] = integrand_matrix[j-1,:] + ((integrand_matrix[j,:] - integrand_matrix[j-1,:])/d_theta*remainder)
        captial_matrix_theta[i,:] = captial_matrix_theta[i,:] -1*np.pi*integration_theta_list[j-1,:] + (integrand_matrix[j-1,:]+ lowercase_tensor_theta[i,:])*remainder/2
        break
        
  return captial_matrix_theta, lowercase_tensor_theta

  
