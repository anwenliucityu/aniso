import numpy as np
import matplotlib.pyplot as plt
from input_dict import Al_screw_config as config
from elem_dict import elem_dict
import disl_theory


def contourf_field():
  timestep = 100
  timestep_sq = timestep * timestep
  x = np.linspace(-63*0.5 * np.sqrt(6),63*0.5 * np.sqrt(6),timestep)
  y = np.linspace(-36*np.sqrt(3),36*np.sqrt(3),timestep)
  X, Y = np.meshgrid(x, y)
  x_coor = X.reshape([timestep_sq,1])
  y_coor = Y.reshape([timestep_sq,1])
  atom_coor = np.zeros([timestep_sq,2])
  atom_coor[:,0] += x_coor[:,0]
  atom_coor[:,1] += y_coor[:,0]
  r = np.sqrt(atom_coor[:,1]**2+atom_coor[:,0]**2).reshape([timestep_sq,1])
  theta = np.arctan2(y_coor, x_coor).reshape([timestep_sq,1])
  '''
  elast_const_initial = np.array([1.001,  1,  1,   0,   0,   0,
                                      1.001,  1,   0,   0,   0,
                                          1.001,   0,   0,   0,
                                                  20,   0,   0,
                                                       20,   0,
                                                            20])
  '''
  elem_property = elem_dict[config["elem"]]
  c11 = elem_property["C11"]
  c12 = elem_property["C12"]
  c44 = elem_property["C44"]
  elast_const_initial = np.array([2*c12+c44,  c12,  c12,   0,   0,   0,
                                        2*c12+c44,  c12,   0,   0,   0,
                                              2*c12+c44,   0,   0,   0,
                                                   c44,   0,   0,
                                                        c44,   0,
                                                             c44])  
  
 

  elem_property = elem_dict[config["elem"]]
  disl_center, b_screw, b_edge = \
            disl_theory.initialize_disl_config(config, elem_property)
  
  sigma = disl_theory.iso_disl_stress(atom_coor, elem_property,
                                       disl_center, b_screw, b_edge)
                                       
  sigma11 = sigma[:,3].reshape(timestep,timestep)
  u_displacement_field = disl_theory.iso_disl_displ(atom_coor, elem_property,
                                            disl_center, b_screw, b_edge)
  u1 = u_displacement_field[:,2].reshape(timestep,timestep)
  
  print(np.min(np.abs(u1)))
  plt.contourf(X,Y, u1, cmap='rainbow')#,levels = np.linspace(-3e8,3e8,200), extend = 'both')
  plt.colorbar()

  atom_sigma11  =  np.sqrt(2)/2*4.05/2/np.pi/0.75 \
                * y_coor * (3.*x_coor**2 + y_coor**2) / r**4
  atom_sigma111 = atom_sigma11.reshape(100,100)
  #print(np.max(atom_sigma111))

  aaaa = np.max(atom_sigma111 - sigma11)
  #print(aaaa)
  #plt.contourf(X, Y, atom_sigma111, cmap = 'rainbow',levels = np.linspace(-0.6,0.6,20), extend = 'both')
  #plt.colorbar()
  plt.show()


def difference():  
  timestep = 100
  timestep_sq = timestep * timestep
  x = np.linspace(-5,5.0001,timestep)
  y = np.linspace(-5,5.0001,timestep)
  X, Y = np.meshgrid(x, y)
  x_coor = X.reshape([timestep_sq,1])
  y_coor = Y.reshape([timestep_sq,1])
  atom_coor = np.zeros([timestep_sq,2])
  atom_coor[:,0] += x_coor[:,0]
  atom_coor[:,1] += y_coor[:,0]
  r = np.sqrt(atom_coor[:,1]**2+atom_coor[:,0]**2).reshape([timestep_sq,1])
  theta = np.arctan2(y_coor, x_coor).reshape([timestep_sq,1])
  frame_initial = config["frame_initial"]
  frame_new = config["frame_new"]
  frame_initial /= np.linalg.norm(frame_initial, axis=0)
  frame_new /= np.linalg.norm(frame_new, axis=0)

  elast_const_initial = np.array([3,  1,  1,   0,   0,   0,
                                      3,  1,   0,   0,   0,
                                          3,   0,   0,   0,
                                               1,   0,   0,
                                                    1,   0,
                                                         1])
  elastic_const_new = elast_const_initial
  elem_property = elem_dict[config["elem"]]
  disl_center, b_vector = initialize_disl_config(config, elem_property, frame_new)
  integrate_stepnumber = 100
  theta_list, q, s, b = integrand(elastic_const_new, integrate_stepnumber)
  S, S_list = zero_pi_integrate(s)
  Q, Q_list = zero_pi_integrate(q)
  B, B_list = zero_pi_integrate(b)
  S_theta, s_theta = zero_theta_integrate(s, r, theta, S_list)
  Q_theta, q_theta = zero_theta_integrate(q, r, theta, Q_list)

  sigma = stress_field(elast_const_initial, atom_coor, r, theta, S, B, s_theta, q_theta, b_vector) 
  atom_sigma11  =  np.sqrt(2)/2*4.05/2/np.pi/0.75 \
                * y_coor * (3.*x_coor**2 + y_coor**2) / r**4  
  x = np.linspace(-1,1,timestep_sq).reshape(10000,1)
  plt.plot(x, sigma[:,0]-atom_sigma11)
  plt.show()

if __name__ == '__main__':
  contourf_field()
  #difference()
