import input_dict
from perfect_cryst_constructor import perfect_cryst_constructor, perfect_hcp_constructor
import sample_region
from elast_const_transform import elast_const_transform
from integrand import integrand
from integrate import zero_pi_integrate, zero_theta_integrate
from aniso_disl_theory import initialize_disl_config, stress_field, displacement_field, partial_init
import pbc_wrap
import numpy as np
import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__,'../../write_output')))
from write_output import write_datafile, write_cfg_auxiliary, write_in_file, write_job, output_dir, write_in_finite_T_meta_file, write_emin_py
import matplotlib.pyplot as plt
from euler_rotation import atom_torsion
import multiprocessing

def aniso_disl_constructor(main_path, config_style, dislocation_type, 
                           pot_element, latt_const, pot_name, element_struct, unit_cell_size,
                           simulation_type, start_temp,
                           elastic_const, mass, pot_cutoff, in_pot, pot_type, pot_path,
                           disl_center,
                           partition, mem, module_load, appexe, ncpu,
                           running_steps=None, dump_interval=None,
                           boundary_freeze_width = 3.0, calc_atomic_stress=True, calc_aniso_stress=False,
                           global_emin=False, temp = 300, msd=None, output_perfect=True,
                           cooling_rate=500, sbatch_job=False, ovito=True, min_step=2000,
                           partial_dislocation=False, partial_split_dis=6, zeroK_latt_const=None):

    '''
    main function for constructing an anisotropic dislocation configuration
    '''
    # import info module
    local = locals()
    exec(f'config = input_dict.{element_struct}_{dislocation_type}')
    config = local['config']

    # Creat a new directory path for installing output file
    directory = output_dir(main_path, pot_element, pot_name, element_struct, config_style,
            dislocation_type, simulation_type, unit_cell_size, 
            partial_dislocation=partial_dislocation, global_emin=global_emin, T=start_temp, temp=temp)
    
    # construct a perfect crystal
    if element_struct == 'hcp':
        atom_coor, atom_type, box_boundary, new_box_lattice_const, frame_new, repeat_para = \
                perfect_hcp_constructor(config, config_style, latt_const, pot_element, 
                        dislocation_type, unit_cell_size)
    if element_struct == 'fcc' or element_struct == 'bcc':
        atom_coor, atom_type, box_boundary, repeat_para = perfect_cryst_constructor(config, config_style,
                element_struct, dislocation_type, latt_const, unit_cell_size,)
        frame_new = config["frame_new"]
  
    # wrap all atoms into one period by PBC
    atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
  
    # select a region of sample
    sample_center = [0, 0]
    sample_radius = np.min(np.abs([box_boundary[0], box_boundary[1]])) - 10
    atom_index_selected = sample_region.cylinder_z(atom_coor, sample_center, sample_radius, 'out')
    
    # delete the selected atoms
    if True:
    	atom_coor = np.delete(atom_coor, atom_index_selected, 0)
    	atom_type = np.delete(atom_type, atom_index_selected)
  
    # write perfect configuration for dislocation analysis 
    if output_perfect == True:
        if unit_cell_size[2] == 1:
            write_datafile(pot_element, dislocation_type, mass, directory, 
        	atom_coor, atom_type, box_boundary, suffix='_perfect_ref.dat')
        else:
            dup_atom_coor, dup_atom_type, dup_box_boundary = pbc_wrap.duplicate_z(atom_coor, 
                    atom_type, box_boundary, unit_cell_size[2])
            write_datafile(pot_element, dislocation_type, mass, directory,
            dup_atom_coor, dup_atom_type, dup_box_boundary, suffix='_perfect_ref.dat')
       
    # normalize the initial and new coordinate
    frame_initial = config["frame_initial"]
    frame_initial /= np.linalg.norm(frame_initial, axis=0)
    frame_new /= np.linalg.norm(frame_new, axis=0)

    # new elastic constant under new rotated coordinate system
    elastic_const_new = elast_const_transform(elastic_const, frame_initial, frame_new)
  
    # dislocation parameter
    frame_crystal = config["frame_initial"]
    if element_struct == 'hcp':
        disl_center, b_vector = initialize_disl_config(config, config_style, latt_const, element_struct, 
                frame_crystal, new_box_lattice_const, disl_center, repeat_para=repeat_para) 
    else:
        disl_center, b_vector = initialize_disl_config(config, config_style, latt_const, element_struct,
                frame_new, latt_const, disl_center) 
  
    # calculate s_theta, q_theta, S, Q, B, S_theta, Q_theta
    integrate_stepnumber = 100
    theta_list, q, s, b = integrand(elastic_const_new, integrate_stepnumber)
    S, S_list = zero_pi_integrate(s)
    Q, Q_list = zero_pi_integrate(q)
    B, B_list = zero_pi_integrate(b)
    
    # stress field
    # write cfg file with atomic stress as the auxiliary property
    if calc_aniso_stress == True:
        atom_sigma = stress_field(elastic_const_new, atom_coor, disl_center, S, B, s_theta, q_theta, b_vector)
        write_cfg_auxiliary(atom_coor, pot_element, dislocation_type, mass, box_boundary, directory,
                            s_11=atom_sigma[:,0], s_12=atom_sigma[:,1], s_13=atom_sigma[:,2],
                            s_22=atom_sigma[:,4], s_23=atom_sigma[:,5], s_33=atom_sigma[:,8])
    
    # displacement
    if partial_dislocation == True:
        # initialize
        b_vectors, disl_centers = partial_init(disl_center, partial_split_dis, b_vector, dislocation_type)
        S_theta, s_theta = zero_theta_integrate(s, atom_coor, disl_centers[0], S_list)
        Q_theta, q_theta = zero_theta_integrate(q, atom_coor, disl_centers[0], Q_list)
        u_displacement_field_1 = displacement_field(atom_coor, disl_centers[0], S, B, S_theta, Q_theta, b_vectors[0])
        # add one partial
        S_theta, s_theta = zero_theta_integrate(s, atom_coor, disl_centers[1], S_list)
        Q_theta, q_theta = zero_theta_integrate(q, atom_coor, disl_centers[1], Q_list)
        u_displacement_field_2 = displacement_field(atom_coor, disl_centers[1], S, B, S_theta, Q_theta, b_vectors[1])
        # add two partials
        u_displacement_field = u_displacement_field_1 + u_displacement_field_2
    if partial_dislocation == False:
        S_theta, s_theta = zero_theta_integrate(s, atom_coor, disl_center, S_list)
        Q_theta, q_theta = zero_theta_integrate(q, atom_coor, disl_center, Q_list)
        u_displacement_field = displacement_field(atom_coor, disl_center, S, B, S_theta, Q_theta, b_vector)
    atom_coor += u_displacement_field
  
    # wrap all atoms into one period by PBC
    atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
  
    # write datafile
    if True:
        if unit_cell_size[2] == 1:
            write_datafile(pot_element, dislocation_type, mass, directory, 
                atom_coor, atom_type, box_boundary)
        else:
            dup_atom_coor, dup_atom_type, dup_box_boundary = pbc_wrap.duplicate_z(atom_coor,
                    atom_type, box_boundary, unit_cell_size[2])
            write_datafile(pot_element, dislocation_type, mass, directory,
            dup_atom_coor, dup_atom_type, dup_box_boundary)
  
    # write in file for lammps and return the directory path
    shell_radius = sample_radius - boundary_freeze_width * pot_cutoff # the radius of region of free-moved atoms
  
    if simulation_type == 'energy_minimization':
        write_in_file(pot_path, config_style, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, sample_center, disl_center, calc_atomic_stress=calc_atomic_stress, global_emin=global_emin,temp=temp, shell_radius = shell_radius,
                  cooling_rate=cooling_rate)
    if simulation_type == 'metastable' or simulation_type == 'finite_T':
        write_in_finite_T_meta_file(pot_path, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, sample_center, disl_center, simulation_type, start_temp, config_style,
                  msd = msd, running_steps=running_steps, 
                  dump_interval=dump_interval, calc_atomic_stress=False, temp = temp, shell_radius = shell_radius,)
        atom_num = atom_coor.shape[0]
        if simulation_type == 'metastable':
            write_emin_py(directory, config_style, mass, in_pot, mem, 
                          partition, module_load, appexe, pot_path, pot_name, atom_num,
                          pot_element, dislocation_type, start_temp, simulation_type,
                          shell_radius=shell_radius, sample_center = sample_center, suffix='quench')
        if simulation_type == 'finite_T':
            if element_struct == 'hcp':
                _, _, zeroK_box, _, _, _ = \
                             perfect_hcp_constructor(config, config_style, zeroK_latt_const, pot_element, 
                                  dislocation_type, unit_cell_size)
            if element_struct == 'fcc' or element_struct == 'bcc':
                _, _, zeroK_box, _ = perfect_cryst_constructor(config, config_style,
                                        element_struct, dislocation_type, zeroK_latt_const, unit_cell_size,)
            if unit_cell_size[2] != 1:
                _, _, zeroK_box = pbc_wrap.duplicate_z(atom_coor, atom_type, zeroK_box, unit_cell_size[2])
            write_emin_py(directory, config_style, mass, in_pot, mem, 
                          partition, module_load, appexe, pot_path, pot_name, atom_num, 
                          pot_element, dislocation_type, start_temp, simulation_type, 
                          min_step=min_step, shell_radius=shell_radius, sample_center = sample_center, suffix='emin')
            write_emin_py(directory, config_style, mass, in_pot, mem, 
                          partition, module_load, appexe, pot_path, pot_name, atom_num, 
                          pot_element, dislocation_type, start_temp, simulation_type, 
                          shell_radius=shell_radius, sample_center = sample_center, suffix = 'quench', zeroK_box = zeroK_box)
  
    # write job.sh file for submitting job
    write_job(pot_element, dislocation_type, partition, mem, module_load, appexe, directory, ncpu=ncpu)

    # print burgers length
    init_path = os.getcwd()
    os.chdir(directory)
    os.system(f'printf "{np.linalg.norm(b_vector)}" > b_length.dat')
    os.chdir(init_path)
    
    # ovito
    if ovito == True:
        os.chdir(directory)
        os.system(f"ovito {pot_element}_{dislocation_type}_disl.dat")  
        os.chdir(init_path)
    
    # cd job.sh directory to submit job
    if sbatch_job == True:
        os.chdir(directory)
        os.system('sbatch job.sh')
        os.chdir(init_path)
    
    return directory

    
if __name__ == '__main__':
    aniso_disl_constructor([0.25,0.25])
