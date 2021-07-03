import numpy as np
import sympy as sym

def homo_strain(S, b_vector, A):
    strain = np.empty([3,3])
    for i in range(3):
        for j in range(3):
            strain[i,j] = -(b_vector[i]*A[j]+b_vector[j]*A[i])/(2*S)
    return strain

def deformation_gradient(
        strain,
        F11 = sym.symbols('F11'), F12 = sym.symbols('F12'), F13 = sym.symbols('F13'),
        F21 = sym.symbols('F21'), F22 = sym.symbols('F22'), F23 = sym.symbols('F23'),
        F31 = sym.symbols('F31'), F32 = sym.symbols('F32'), F33 = sym.symbols('F33'),
        ):
    F = sym.Matrix(
            [[F11, F12, F13],
             [F21, F22, F23],
             [F31, F32, F33]])
    equation =1/2 *(F.T + F) - sym.eye(3) - np.array(strain)
    para = [F11, F12, F13, F21, F22, F23, F31, F32, F33]
    var = []
    for i in para:
        if isinstance(i, int) == False:
            if isinstance(i,float) == False:
                var.append(i)
    solution = sym.solve(equation,var)
    
    deform_grad = []
    for i in range(len(para)):
        if isinstance(para[i],int) == True or isinstance(para[i], float) == True:
            deform_grad.append(para[i])
        else:
            if i == 0:
                deform_grad.append(float(solution[F11]))
            if i == 1:
                deform_grad.append(float(solution[F12])) 
            if i == 2:
                deform_grad.append(float(solution[F13]))
            if i == 3:
                deform_grad.append(float(solution[F21]))
            if i == 4:
                deform_grad.append(float(solution[F22]))
            if i == 5:
                deform_grad.append(float(solution[F23]))
            if i == 6:
                deform_grad.append(float(solution[F31]))
            if i == 7:
                deform_grad.append(float(solution[F32]))
            if i == 8:
                deform_grad.append(float(solution[F33]))
        
    deform_grad = np.array(deform_grad).reshape((3,3))
    return deform_grad

def deform_atom_coor(deform_grad, atom_coor):
    atom_coor = np.inner(atom_coor, deform_grad)
    return atom_coor

def deform_box_boundary(deform_grad, box_boundary, tilt):
    xlo, xhi = box_boundary[0]
    ylo, yhi = box_boundary[1]
    zlo, zhi = box_boundary[2]
    x = [xhi-xlo, 0, 0]
    y = [0, yhi-ylo, 0]
    z = [0, 0, zhi-zlo]
    new_x = np.dot(deform_grad, x)
    new_y = np.dot(deform_grad, y)
    new_z = np.dot(deform_grad, z)
    new_boundary = [[xlo, xlo+new_x[0]],
                    [ylo, ylo+new_y[1]],
                    [zlo, zlo+new_z[2]]]
    if new_y[0] == 0:
        new_y[0] = 0.
    if new_z[0] == 0:
        new_z[0] = 0.
    if new_z[1] == 0:
        new_z[1] = 0.
    tilt_para = [new_y[0]+tilt[0], new_z[0]+tilt[1], new_z[1]+tilt[2]]
    return new_boundary, tilt_para

def rotate_atom_coor(atom_coor, rotation_matrix):
    atom_coor = np.inner(atom_coor, rotation_matrix)
    return atom_coor

def rotate_box_boundary(box_boundary, rotation_matrix):
    xlo, xhi = box_boundary[0]
    ylo, yhi = box_boundary[1]
    zlo, zhi = box_boundary[2]
    box_size = [[xhi-xlo, 0, 0],
   	        [0, yhi-ylo, 0],
                [0, 0, zhi-zlo]]
    new_box = np.dot(rotation_matrix, box_size)
    new_boundary = [[0, new_box.sum(axis=0)[2]],
                    [-new_box.sum(axis=0)[0]/2, new_box.sum(axis=0)[0]/2],
                    [-new_box.sum(axis=0)[1]/2, new_box.sum(axis=0)[1]/2]]
    return new_boundary
    
if __name__ =='__main__':
    import numpy as np
    strain = [[2,0,0],
              [0,1,0],
              [0,0,1]]
    deform_grad = deformation_gradient(strain,F21=0,F31=0,F32=0)
    print(deform_grad)
