import atomman as am
import multiprocessing as mp
from functools import partial
import re
import os
import sys
import numpy as np

def output_dir(main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
               simulation_type, unit_cell_size, partial_dislocation=False, global_emin=True, T=0, temp=300):
    temp = str(int(temp)) + 'K'
    T    = str(int(T)) + 'K'
    box  = str(unit_cell_size[0])+'_'+str(unit_cell_size[1])+'_'+ str(unit_cell_size[2])
    dir = os.path.join(main_path, pot_element, pot_name, simulation_type, element_struct, config_style, dislocation_type)
                 # e.g. main/Ti/Ti_kevin_2020/energy_minimization/hcp/cylinder/screw_a_basal/
    #if config_style == 'quadrupolar':
     #   global_emin = False
    if simulation_type == 'energy_minimization':
        if global_emin == True:
            dir = os.path.join(dir, T+'_'+temp+'_'+T,)
            # e.g. energy_minimization/0_300_0
        if global_emin == False:
            dir = os.path.join(dir, T)
    elif simulation_type == 'metastable':
        dir = os.path.join(dir, T+'_'+temp)
    elif simulation_type == 'finite_T':
        dir = os.path.join(dir, T)
    if partial_dislocation == False:
        par = 'full'
    elif partial_dislocation == True:
        par = 'partial'
    dir = os.path.join(dir, box, par)
    print("DDplot and Nye is parepared ......")
    print("work_path = " + dir)
    return dir

# extract the number of files and name of files
def fileprocess(path):
    filenames = []
    filenum = 0
    for lists in os.listdir(path):
        sub_path = os.path.join(path, lists)
        if os.path.isfile(sub_path):
            filenum = filenum+1
            filenames.append(lists)
    # sort the list according to its dump sequence
    filenames.sort(key=lambda x: int(x.split('_')[-1][:-5]))
    return filenum, filenames

def ave_atom_pos_write(file_paths, path, alat, cut_off, plot_name):
    file_num = len(file_paths)
    atompos = 0
    for i in range(file_num):
        aim_path = os.path.join(path, file_paths[i])
        atompos += am.load('atom_dump', aim_path).atoms.pos
    atompos /= file_num

    # read box
    config = am.load('atom_dump', aim_path)
    box = [[config.box.origin[0], config.box.origin[0]+config.box.avect[0]],
           [config.box.origin[1], config.box.origin[1]+config.box.bvect[1]],
           [config.box.origin[2], config.box.origin[2]+config.box.cvect[2]]]
    natoms = atompos.shape[0]
    atom_type = config.atoms.atype
    # pbc
    
    repeat_time = 1
    # repeat z 3 times for nye tensor plot
    if config.box.cvect[2] < cut_off*alat:
        repeat_time = 3
    repeat_unit = box[2][1] - box[2][0]
    new_zhi = box[2][0] + repeat_unit * repeat_time
    box = [[box[0][0], box[0][1]],
               [box[1][0], box[1][1]],
               [box[2][0], new_zhi]]
    new_coor = []
    new_type = []
    shift = np.array([0,0,repeat_unit])
    for j in range(repeat_time):
        new_coor = np.append(new_coor, np.array(atompos) + shift*j)
        new_type = np.append(new_type, atom_type)
    atompos = new_coor.reshape((natoms*repeat_time,3))
    atom_type = new_type.reshape((natoms*repeat_time,1))
    natoms *= repeat_time

    # write ave_data_file
    ave_file_path = os.path.abspath(os.path.join(path, f'../../ave_file_{plot_name}'))
    if not os.path.exists(ave_file_path):
        os.makedirs(ave_file_path)
    file_name = file_paths[i]+ '_ave.dat'
    file_path_name = os.path.join(ave_file_path, file_name)
    file = open(file_path_name,'w')

    file.write(f"\n\n"
               f"{natoms} atoms\n"
               f"1 atom types\n"
               f"{box[0][0]:>14.12f} {box[0][1]:>14.12f} xlo xhi\n"
               f"{box[1][0]:>14.12f} {box[1][1]:>14.12f} ylo yhi\n"
               f"{box[2][0]:>14.12f} {box[2][1]:>14.12f} zlo zhi\n\n")
    file.write("Atoms\n\n")
    for i in range(natoms):
        file.write(f"{i+1} {int(atom_type[i])} {atompos[i,0]} {atompos[i,1]} {atompos[i,2]}\n")
    file.close()
    
def average(latt_const, cutoff_factor, partial_dislocation, plot_name,
            main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
            simulation_type, unit_cell_size ,global_emin=True, T=0, temp=300,
            ave_file_num=1, ave_interval=1, plotting_region=[80,80], ):
    if element_struct =='hcp':
        alat = latt_const['a']
    else:
        alat = latt_const
    path = output_dir(main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
               simulation_type, unit_cell_size ,
               partial_dislocation=partial_dislocation, global_emin=global_emin, T=T, temp=temp)
    if simulation_type == 'finite_T' or simulation_type == 'metastable':
        path = os.path.join(path, f'{plot_name}/emin')
        per_file = os.path.abspath(os.path.join(path, '../../'+pot_element+'_'+dislocation_type+'_perfect_ref.dat'))
    filenum, filename = fileprocess(path)
    atompos = am.load('atom_data', per_file).atoms.pos
    natoms = atompos.shape[0]
  
    pool = mp.Pool()
    for i in range((filenum-ave_file_num)//ave_interval+1):
        rep = ave_interval*i
        file_paths = filename[0+rep:ave_file_num+rep]
        #ave_atom_pos_write(file_paths, path, alat, cutoff_factor, atom_index_selected)
        pool.apply_async(ave_atom_pos_write, (file_paths, path, alat, cutoff_factor, plot_name),)
    pool.close()
    pool.join()
    ave_file_path = os.path.abspath(os.path.join(path, f'../../ave_file_{plot_name}'))
    return ave_file_path
    
def ave_file_path(partial_dislocation, plot_name,
            main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
            simulation_type, unit_cell_size ,global_emin=True, T=0, temp=300,):
    path = output_dir(main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
               simulation_type, unit_cell_size ,
               partial_dislocation=partial_dislocation, global_emin=global_emin, T=T, temp=temp)
    if simulation_type == 'finite_T' or simulation_type == 'metastable':
        path = os.path.join(path, f'{plot_name}/emin')
    ave_file_path = os.path.abspath(os.path.join(path, f'../../ave_file_{plot_name}'))
    return ave_file_path
