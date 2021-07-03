import input_dict
from perfect_cryst_constructor import perfect_cryst_constructor, perfect_hcp_constructor
import sample_region
from elast_const_transform import elast_const_transform
from integrand import integrand
from integrate import zero_pi_integrate, zero_theta_integrate
from aniso_disl_theory import initialize_disl_config, stress_field, displacement_field, displacement_fields, partial_init
import pbc_wrap
import numpy as np
import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__,'../../write_output')))
from write_output import write_datafile, write_cfg_auxiliary, write_in_file, write_job, output_dir, write_in_finite_T_meta_file, write_emin_py
import matplotlib.pyplot as plt
from euler_rotation import atom_torsion
from quadrupolar import quadrupolar_params
import deform_gradient
import multiprocessing

def aniso_disl_constructor(main_path, config_style, dislocation_type,                  
                           pot_element, latt_const, pot_name, element_struct, unit_cell_size,
                           simulation_type, start_temp,
                           elastic_const, mass, pot_cutoff, in_pot, pot_type, pot_path,
                           disl_center,
                           partition, mem, module_load, appexe, ncpu,
                           running_steps=None, dump_interval=None,
                           boundary_freeze_width = 3.0, calc_atomic_stress=True, calc_aniso_stress=False,
                           global_emin=False, temp = 300, output_perfect=True,
                           cooling_rate=500, sbatch_job=False, ovito=True, msd=None, min_step=2000,
                           partial_dislocation=False, partial_split_dis=6):

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
                perfect_hcp_constructor(config, config_style, latt_const, pot_element, dislocation_type, unit_cell_size, quadru=True)
    if element_struct == 'fcc' or element_struct == 'bcc':
        atom_coor, atom_type, box_boundary, repeat_para = perfect_cryst_constructor(config, config_style,
                element_struct, dislocation_type, latt_const, unit_cell_size, quadru=True)
        frame_new = config["frame_new"]

    xy = unit_cell_size[1]*repeat_para[1,0] + unit_cell_size[0]*repeat_para[0,1]
    yz = unit_cell_size[2]*repeat_para[2,1] + unit_cell_size[1]*repeat_para[1,2]
    xz = unit_cell_size[2]*repeat_para[2,0] + unit_cell_size[0]*repeat_para[0,2]
    tilt = [xy,xz,yz]  
    
    # arrange box to satisfy pbc
    box_boundary = pbc_wrap.pbc_wrap_tilt_box(box_boundary, tilt)
    
    # write perfect configuration for dislocation analysis 
    if output_perfect == True:
        write_datafile(pot_element, dislocation_type, mass, directory, 
            atom_coor, atom_type, box_boundary, suffix='_perfect_ref.dat', tilt_para=tilt)
    
    # normalize the initial and new coordinate
    frame_initial = config["frame_initial"]
    frame_initial /= np.linalg.norm(frame_initial, axis=0)
    frame_new /= np.linalg.norm(frame_new, axis=0)
    
    # new elastic constant under new rotated coordinate system
    elastic_const_new = elast_const_transform(elastic_const, frame_initial, frame_new)

    # dislocation parameter
    frame_crystal = config["frame_initial"]
    if element_struct == 'hcp':
        __, b_vector = initialize_disl_config(config, config_style, latt_const, element_struct, 
                frame_crystal, new_box_lattice_const, disl_center, repeat_para=repeat_para) 
    else:
        __, b_vector = initialize_disl_config(config, config_style, latt_const, element_struct,
                frame_new, latt_const, disl_center, repeat_para=repeat_para) 
    
    # calculate s_theta, q_theta, S, Q, B, S_theta, Q_theta
    integrate_stepnumber = 100
    theta_list, q, s, b = integrand(elastic_const_new, integrate_stepnumber)
    S, S_list = zero_pi_integrate(s)
    Q, Q_list = zero_pi_integrate(q)
    B, B_list = zero_pi_integrate(b)
    
    # rotate coor with z // dislocation line
    qua_rot = np.array([[0,0,1.],[1.,0,0],[0,1.,0]])
    atom_coor_ref = np.inner(atom_coor, np.linalg.inv(qua_rot))

    # initialize quadrupolar params (add neighbour image dislocation line to satisfy pbc)
    disl_centers, d, l, A, S_area, plane_point = quadrupolar_params(disl_center, box_boundary, repeat_para, unit_cell_size)
    if partial_dislocation == True:
        disl_list = []
        for i in range(len(disl_centers)):
            b_vec, disl_cens = partial_init(disl_centers[i], partial_split_dis, b_vector, dislocation_type)
            disl_list.append(disl_cens)
    # calculate displacements according to PBC dislocation centers
    pool = multiprocessing.Pool()
    u_displacements = []
    u_errors = []
    for i in range(len(disl_centers)): 
        if partial_dislocation == False:
            u_displacements.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, atom_coor_ref, disl_centers[i], S_list, Q_list, b_vector,)))
            u_errors.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, plane_point, disl_centers[i], S_list, Q_list, b_vector,)))
        if partial_dislocation == True:
            u_displacements.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, atom_coor_ref, disl_list[i][0], S_list, Q_list, b_vec[0],)))
            u_displacements.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, atom_coor_ref, disl_list[i][1], S_list, Q_list, b_vec[1],)))
            u_errors.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, plane_point, disl_list[i][0], S_list, Q_list, b_vec[0],)))
            u_errors.append(pool.apply_async(displacement_fields, 
                    (i, s, q, S, B, plane_point, disl_list[i][1], S_list, Q_list, b_vec[1],)))
    pool.close()
    pool.join()
    
    # calculate u_error and fit an error plane, Ax = b, solute X
    u_err = 0
    for err in u_errors:
        u_err += np.array(err.get())
    d = np.linalg.lstsq(plane_point, u_err, rcond=None)[0]

    # add displacement and rotate back with x // dislocation line
    atom_coor = atom_coor_ref
    for u_displacement_field in u_displacements:
        atom_coor += np.array(u_displacement_field.get())
    atom_coor -= np.inner(atom_coor, d.T)
    atom_coor = np.inner(atom_coor,qua_rot)
    
    # add homogeneous strain to stablized dislocations
    b_vector = np.inner(b_vector, qua_rot)
    
    strain = deform_gradient.homo_strain(S_area,b_vector,A)
    deform_grad = deform_gradient.deformation_gradient(strain, F21=0, F31=0, F32=0)
    atom_coor = deform_gradient.deform_atom_coor(deform_grad, atom_coor)
    box_boundary, tilt_para = deform_gradient.deform_box_boundary(deform_grad, box_boundary, tilt)
    tilt_para += tilt

    # wrap atom coordinate to pbc
    atom_coor = pbc_wrap.pbc_wrap_tilt(atom_coor, box_boundary, tilt_para)

    # write datafile
    write_datafile(pot_element, dislocation_type, mass, directory, 
            atom_coor, atom_type, box_boundary, tilt_para=tilt_para)
    
    # write in file for lammps and return the directory path
    sample_center = None
    if simulation_type == 'energy_minimization':
        write_in_file(pot_path, config_style, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, sample_center, disl_center, calc_atomic_stress=calc_atomic_stress, global_emin=global_emin,temp=temp,
                  cooling_rate=cooling_rate, box_boundary = box_boundary, tilt = tilt_para)
    if simulation_type == 'metastable' or simulation_type == 'finite_T':
        write_in_finite_T_meta_file(pot_path, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, sample_center, disl_center, simulation_type, start_temp, config_style,
                  running_steps=running_steps,
                  dump_interval=dump_interval, calc_atomic_stress=False, temp = temp, 
                  box_boundary = box_boundary, tilt = tilt_para)
        atom_num = atom_coor.shape[0]
        if simulation_type == 'metastable':
            write_emin_py(directory, config_style, mass, in_pot, mem, 
                          partition, module_load, appexe, pot_path, pot_name, atom_num,
                          pot_element, dislocation_type, start_temp, simulation_type,
                          box_boundary = box_boundary, tilt = tilt_para)
        if simulation_type == 'finite_T':
            write_emin_py(directory, config_style, mass, in_pot, mem, 
                          partition, module_load, appexe, pot_path, pot_name, atom_num, 
                          pot_element, dislocation_type, start_temp, simulation_type,
                          min_step=min_step, box_boundary = box_boundary, tilt = tilt_para)
   
    # write job.sh file for submitting job
    write_job(pot_element, dislocation_type, partition, mem, module_load, appexe, directory, ncpu=ncpu)
    
    # print burgers length
    init_path = os.getcwd()
    os.chdir(directory)
    os.system(f'printf "{np.linalg.norm(b_vector)}" > b_lenght.dat')
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
    aniso_disl_constructor()
