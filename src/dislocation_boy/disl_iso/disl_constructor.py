'''Construct a single straight dislocation'''

import numpy as np
import matplotlib.pyplot as plt
import os
from elem_dict import elem_dict
from perfect_cryst_constructor import perfect_cryst_constructor
import sample_region
import disl_theory
import pbc_wrap
import disl_plot
from write_output import write_datafile, write_cfg_auxiliary, write_in_file, write_job, write_info, output_directory
from input_dict import Pt_screw_config as config
from check_lmp import check_lmp
import time
from euler_rotation import atom_torsion
import multiprocessing

def disl_constructor():
    '''
    main function for constructing a dislocation
    '''

    # material parameters
    elem_property = elem_dict[config["elem"]]

    # construct a perfect crystal
    atom_coor, atom_type, box_boundary = \
            perfect_cryst_constructor(config, elem_property)
    
    # select a region of sample
    sample_center = [0,0]
    sample_radius = np.min(np.abs([box_boundary[0], box_boundary[1]])) - 10.
    atom_index_selected = sample_region.cylinder_z(atom_coor, sample_center,
                                                   sample_radius, 'out')

      
    # delete the selected atoms
    atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    atom_type = np.delete(atom_type, atom_index_selected)
    
    sample_radius_new = 12
    atom_index_selected = sample_region.cylinder_z(atom_coor, sample_center, sample_radius_new, 'in')

    atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    atom_type = np.delete(atom_type, atom_index_selected)
    # select crack region
    '''
    x_range = [83, 100]
    y_range = [float("-inf"), float("inf")]
    z_range = [-3.5, 3.5]
    atom_index_selected = sample_region.crack_tip_region(atom_coor, x_range, y_range, z_range)
    '''
    '''
    # triangle region
    degree = (np.pi/6) # pi/6 [63,0.01]  ## pi/9 [70, 0.01]
    start_point = [63,0.01]
    atom_index_selected = sample_region.crack_tip_angle_region(atom_coor, degree, start_point)
    atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    atom_type = np.delete(atom_type, atom_index_selected)
    '''
    # introduce a dislocation
    # set dislocation configuration
    disl_center, b_screw, b_edge = \
            disl_theory.initialize_disl_config(config, elem_property)
    print(disl_center, b_screw, b_edge)
    # calculate displacement by isotropic elasticity
    # then, displace atoms by the displacement predicted by theory
    atom_coor += disl_theory.iso_disl_displ(atom_coor, elem_property,
                                            disl_center, b_screw, b_edge)
    
    '''
    # create one edge dislocation
    disl_center = [0.01,60.0]

    atom_displ = disl_theory.iso_disl_displ_x(atom_coor, elem_property,
                                            disl_center, b_edge, b_screw)
       
    # PBC
    
    for i in range(atom_displ.shape[0]):
    	for j in range(0,3):
    		atom_displ[i, j] = atom_displ[i,j] * np.exp(-np.abs(atom_coor[i,2])/25)
    
    	
    
    atom_coor += atom_displ
    '''
    # wrap all atoms into one period by PBC
    atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)

    
    '''
    # calculate and plot stress field based on isotropic elasticity theory
    atom_sigma = disl_theory.iso_disl_stress(atom_coor, elem_property,
                                       disl_center, b_screw, b_edge)
    '''
    
    '''
    # change the atom type in a selected region
    shell_radius = sample_radius - 4* elem_property["pot_cutoff"]
    #atom_index_selected = sample_region.cylinder_z(atom_coor, sample_center,
     #                                              shell_radius, 'out')
    atom_index_selected = sample_region.upper_lower_z(atom_coor, box_boundary, 20, 'ud')
    new_type = 2
    atom_type[atom_index_selected] = new_type
    
    #add crack
    degree = [np.pi/6, 0] 
    start_point = [60,0.1]
    # change the crack direction with "down or up"
    atom_index_selected = sample_region.crack_tip_angle_region(atom_coor, degree, start_point, 'down')
    '''
    '''
    center = [0.,58.0, 0.1]
    direction = "left"
    atom_index_selected = sample_region.cone_crack(atom_coor, center, direction)
    atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    atom_type = np.delete(atom_type, atom_index_selected)
    '''
    # Creat a new directory path for installing output file
    directory = output_directory(config, elem_property, '/gauss12/home/cityu/anwenliu/dislocation/non_torsion/')
    
    # write datafile
    box_boundary[0][0] =-200
    box_boundary[0][1] = 200
    box_boundary[1][0] =-200
    box_boundary[1][1] = 200
    write_datafile(config, directory, atom_coor, atom_type, box_boundary)
    
    '''
    # write cfg file with atomic stress as the auxiliary property
    write_cfg_auxiliary(atom_coor, config, elem_property, box_boundary, directory,
                          s_11=atom_sigma[:,0], s_12=atom_sigma[:,1], s_13=atom_sigma[:,2],
                          s_22=atom_sigma[:,3], s_23=atom_sigma[:,4], s_33=atom_sigma[:,5])
                       
    # write in file for lammps and return the directory path
    write_in_file(config, elem_property, directory)

    # write job.sh file for submitting job
    write_job(config, directory, 12) # core number
  
    # write basic information about the system
    write_info(config, elem_property, atom_coor, atom_index_selected, sample_radius, shell_radius, box_boundary, directory)
    
    # cd job.sh directory to submit job, and calculate
    os.chdir(directory)
    os.system('sbatch job.sh')
    
    # check if job is finished, once finished, open with ovito
    start = time.time()
    check_lmp(config, 5) # check if finished every 5 secs.
    end = time.time()
    timecost = str(end-start)
    print("lammps finished! The processing time is " + timecost + ' secs')
    '''
    '''
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.scatter3D(atom_coor[:,0],atom_coor[:,1],atom_coor_check[:,1],'gray') 
    #ax1.plot3D(atom_coor[0,:],atom_coor[1,:],atom_coor[2,:],'gray')    #绘制空间曲线
    plt.show()
    '''
if __name__ == '__main__':
    disl_constructor()

