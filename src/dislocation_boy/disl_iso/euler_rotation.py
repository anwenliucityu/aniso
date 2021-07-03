import numpy as np
import math as m
  
def Rx(theta):
  return np.array([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])
  
def Ry(theta):
  return np.array([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]])
  
def Rz(theta):
  return np.array([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

def atom_torsion(box_boundary_z_len, rotation_angle, sample_center, atom_coor):
  z_coor = atom_coor[2] + box_boundary_z_len/2
  atom_coor = atom_coor - [sample_center[0], sample_center[1], 0]
  atom_rotation_angle = z_coor/box_boundary_z_len * rotation_angle
  atom_coor = np.dot(Rz(atom_rotation_angle), atom_coor) + [sample_center[0], sample_center[1], 0]
  return atom_coor
  
             
if __name__ == "__main__":
  a = np.array([1,0,0])
  c = np.dot(Rz(np.pi/2), a)
  print(Ry(np.pi))
