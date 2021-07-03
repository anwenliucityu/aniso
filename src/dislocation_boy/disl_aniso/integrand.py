from voigt_notation import voigt_one_to_two, voigt_two_to_one, \
                           four_tensor_one_to_two, four_tensor_two_to_one,\
                           four_tensor_matrix_to_list,\
                           four_tensor_list_to_matrix
import numpy as np
from sympy import Eijk
import itertools

def integrand (elastic_const, stepnumber):
  '''
  integrand function calculates the mm, nn, nm , q, s, b matrix with regard to
  a list of theta in the range [0, pi)
  ----------------------
  varible definition
  elastic_const: the calculated elastic constant after coordinate system rotation
                 the shape is (1,21)
  stepnumber: the stepnumber decides how many intervals in the range [0, pi)
              with bigger stepnumber, we get more accurate integral results and
              longer time for calculation
  '''
  idx11 = []
  idx12 = []
  idx21 = []
  idx22 = []
  for j,k in itertools.product(range(3), repeat=2):
    idx11.append(four_tensor_two_to_one[(voigt_two_to_one[(0, j)], voigt_two_to_one[(k, 0)])])
    idx12.append(four_tensor_two_to_one[(voigt_two_to_one[(0, j)], voigt_two_to_one[(k, 1)])])
    idx21.append(four_tensor_two_to_one[(voigt_two_to_one[(1, j)], voigt_two_to_one[(k, 0)])])
    idx22.append(four_tensor_two_to_one[(voigt_two_to_one[(1, j)], voigt_two_to_one[(k, 1)])])

  # list all c_1jk1, c_1jk2, c_2jk1, c_2jk2
  c_1jk1 = []
  c_1jk2 = []
  c_2jk1 = []
  c_2jk2 = []
  for i in range(9):
    c_1jk1.append(elastic_const[idx11[i]])
    c_1jk2.append(elastic_const[idx12[i]])
    c_2jk1.append(elastic_const[idx21[i]])
    c_2jk2.append(elastic_const[idx22[i]])
  c_1jk1 = np.array(c_1jk1).reshape(1,9)
  c_1jk2 = np.array(c_1jk2).reshape(1,9)
  c_2jk1 = np.array(c_2jk1).reshape(1,9)
  c_2jk2 = np.array(c_2jk2).reshape(1,9)

  # list cos^2(theta), sin^(theta), sin(theta)*cos(theta)
  theta_list = np.linspace(0, np.pi, stepnumber)
  cos2 = (np.cos(theta_list)*np.cos(theta_list)).reshape([stepnumber,1])
  sin2 = (np.sin(theta_list)*np.sin(theta_list)).reshape([stepnumber,1])
  sincos = (np.sin(theta_list)*np.cos(theta_list)).reshape([stepnumber,1])

  # matrix in form M * 9, every M represent the components value for theta, 9 reps the 9 components (i=1,2,3, j=1,2,3)
  mm_mat = c_1jk1 * cos2 + (c_1jk2 + c_2jk1) * sincos + c_2jk2 * sin2
  nn_mat = c_1jk1 * sin2 - (c_1jk2 + c_2jk1) * sincos + c_2jk2 * cos2
  nm_mat =-c_1jk2 * sin2 + (c_2jk2 - c_1jk1) * sincos + c_2jk1 * cos2

  q_mat = np.zeros(shape=(stepnumber,9))
  denominator = np.zeros(shape=(stepnumber,1))
  for p,g,n in itertools.product(range(3), repeat = 3):
    p1 = four_tensor_matrix_to_list[(0, p)]
    g2 = four_tensor_matrix_to_list[(1, g)]
    n3 = four_tensor_matrix_to_list[(2, n)]
    denominator[:,0] = denominator[:,0] + 2 * Eijk(p,g,n) * nn_mat[:,p1] * nn_mat[:,g2] * nn_mat[:,n3]

  for i,s,m,j,r,w in itertools.product(range(3), repeat = 6):
    ij = four_tensor_matrix_to_list[(i, j)]
    sr = four_tensor_matrix_to_list[(s, r)]
    mw = four_tensor_matrix_to_list[(m, w)]
    q_mat[:,ij] = q_mat[:,ij] + (Eijk(i,s,m) * Eijk(j,r,w) * nn_mat[:,sr] * nn_mat[:,mw]) / denominator[:,0]
  
  s_mat = np.zeros(shape = (stepnumber,9))
  for i,j,k in itertools.product(range(3), repeat = 3):
    ij = four_tensor_matrix_to_list[(i, j)]
    ik = four_tensor_matrix_to_list[(i, k)]
    kj = four_tensor_matrix_to_list[(k, j)]
    jk = four_tensor_matrix_to_list[(j, k)]
    s_mat[:, ij] = s_mat[:, ij] + q_mat[:, ik] *nm_mat[:, kj]

  b_mat = np.zeros(shape = (stepnumber, 9))
  for i,j in itertools.product(range(3), repeat = 2):
    ij = four_tensor_matrix_to_list[(i, j)]
    b_mat[:,ij] = b_mat[:,ij] - mm_mat[:, ij]
    for k in range(3):
      kj = four_tensor_matrix_to_list[(k, j)]
      ki = four_tensor_matrix_to_list[(k, i)]
      b_mat[:,ij] = b_mat[:,ij] + nm_mat[:,ki] * s_mat[:, kj]

  return theta_list, q_mat, s_mat, b_mat
