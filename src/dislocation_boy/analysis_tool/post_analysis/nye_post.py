import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import multiprocessing as mp
from scipy.optimize import fminbound

sym = {
       'screw_a_basal'  :   4,
        }

def nye_standard_deviation(start_temperature, path, dislocation_type, latt_const, shift=5):
    info = pd.read_csv(path, delim_whitespace = True).sort_values(by='id')
    symmetry = sym[dislocation_type]
    a = latt_const['a']
    c = latt_const['c']
    xx = np.array(info['xx'])
    yy = np.array(info['yy'])
    xy = np.array(info['xy'])
    rotate_theta = []
    for i in range(xx.shape[0]):
        matrix = np.array([[xx[i],xy[i]],[xy[i],yy[i]]])
        values, vectors = np.linalg.eig(matrix)
        maxvalue_index = np.argmax(values)
        vector = abs(vectors[maxvalue_index])
        angle = np.arctan2(vector[1],vector[0])
        rotate_theta.append(np.degrees(angle))
    if dislocation_type == 'screw_a_basal':
        planes = [0, np.arctan(c/a)*180/np.pi, 90]
        plane_name = ['basal', 'pyrI', 'prism']
    write_path = os.path.abspath(os.path.join(path,'../analysis_sd.dat'))
    file = open(write_path,'w')
    file.write(f'temperarure  {start_temperature}\n')
    file.write('spread_orientation    number     possibility\n')
    c = 0
    for i in range(len(planes)):
        counts = 0
        for j in range(len(rotate_theta)):
            if rotate_theta[j]<planes[i]+shift and rotate_theta[j]>planes[i]-shift:
                counts +=1
        file.write(f'{plane_name[i]:<10} {planes[i]:>3.3f}     {counts:>5}       {counts/len(rotate_theta)*100:>3.3f}%\n')
        c+=counts
    file.write(f'total           {c}        {c/len(rotate_theta)*100:>3.3f}%\n')
    file.close()

    


