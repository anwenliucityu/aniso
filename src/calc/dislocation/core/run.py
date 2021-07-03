#!/usr/bin/python3.8
import os
import sys
import importlib
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../.."))) #yes
import default_path as path
# Read information in spec file and path
main_path = path.output_path
import multiprocessing
spec      = __import__(sys.argv[2].replace('.py',''))
pot_path  = path.pot_path
sys.path.append(os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id)))
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/disl_aniso")))
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../potential",spec.pot_element,spec.pot_id)))
import pot_mod
# find potential file, if file is in package then use it, otherwise find outer file
for i in range(len(pot_mod.pot_path)):
    if os.path.exists(pot_mod.pot_path[i]) == False:
        __, pot_name = os.path.split(pot_mod.pot_path[i])
        pot_mod.pot_path[i] = os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id, pot_name))


#####################
#   Main Function   #
#####################

# Customize
params = {}
if 'boundary_freeze_width' in dir(spec):
    params['boundary_freeze_width'] = spec.boundary_freeze_width
if 'calc_atomic_stress' in dir(spec):
    params['calc_atomic_stress']    = spec.calc_atomic_stress
if 'global_emin' in dir(spec):
    params['global_emin']           = spec.global_emin
if 'cooling_rate' in dir(spec):
    params['cooling_rate']          = spec.cooling_rate
if 'calc_aniso_stress' in dir(spec):
    params['calc_aniso_stress']     = spec.calc_aniso_stress
if 'temp' in dir(spec):
    params['temp']                  = spec.temp
if 'running_steps' in dir(spec):
    params['running_steps']         = spec.running_steps
if 'dump_interval' in dir(spec):
    params['dump_interval']         = spec.dump_interval
if 'sbatch_job' in dir(spec):
    params['sbatch_job']            = spec.sbatch_job
if 'output_perfect' in dir(spec):
    params['output_perfect']        = spec.output_perfect
if 'ovito' in dir(spec):
    params['ovito']                 = spec.ovito
if 'min_step' in dir(spec):
    params['min_step']              = spec.min_step
params['partial_dislocation']                   = spec.partial_dislocation
if 'partial_split_dis' in dir(spec):
    params['partial_split_dis']         = spec.partial_split_dis

if isinstance(spec.start_temperature, int) == True:
            spec.start_temperature = [spec.start_temperature]  
unit_cell_size = [spec.num_unit_cell_x, spec.num_unit_cell_y, spec.num_unit_cell_z]
disl_center    = [spec.disl_center_x, spec.disl_center_y]

# import the module to constructure the final atomic structure
if spec.config_style == 'quadrupolar':
    import quadru_aniso_disl_config_constructor as aniso_disl
else:
    import aniso_disl_config_constructor as aniso_disl
module_path = os.path.abspath(aniso_disl.__file__)

if __name__ == "__main__":
    if spec.simulation_type == 'metastable':
        spec.start_temperature = [0]
    if sys.argv[1] == 'calc' and spec.config_style == 'cylinder':
        exec(f'import info_{spec.element_struct}_0K as mi_zero')
        params['zeroK_latt_const']             = mi_zero.latt_const
        pool = multiprocessing.Pool()
        for i in range(len(spec.start_temperature)):
            exec(f'import info_{spec.element_struct}_{spec.start_temperature[i]}K as mi')
            if 'msd' in dir(mi):
                params['msd']                   = mi.msd
            # convert elastic constant unit from GPa to Pa
            elastic_const = np.array(mi.elastic_const)*1e9  
            pool.apply_async(aniso_disl.aniso_disl_constructor, (main_path, spec.config_style, spec.dislocation_type, 
                           spec.pot_element, mi.latt_const, spec.pot_id, spec.element_struct, unit_cell_size,
                           spec.simulation_type, spec.start_temperature[i],
                           elastic_const, mi.mass, pot_mod.pot_cutoff, pot_mod.in_pot, 
                           pot_mod.pot_type, pot_mod.pot_path,
                           disl_center,
                           spec.partition, spec.mem, pot_mod.module_load, pot_mod.appexe, spec.ncpu,), 
                           params,)

        pool.close()
        pool.join()
	
    if sys.argv[1] == 'calc' and spec.config_style == 'quadrupolar':
        for i in range(len(spec.start_temperature)):
            exec(f'import info_{spec.element_struct}_{spec.start_temperature[i]}K as mi')
            if 'msd' in dir(mi):
                params['msd']                   = mi.msd
            # convert elastic constant unit from GPa to Pa
            elastic_const = np.array(mi.elastic_const)*1e9
            aniso_disl.aniso_disl_constructor(main_path, spec.config_style, spec.dislocation_type, 
                           spec.pot_element, mi.latt_const, spec.pot_id, spec.element_struct, unit_cell_size,
                           spec.simulation_type, spec.start_temperature[i],
                           elastic_const, mi.mass, pot_mod.pot_cutoff, pot_mod.in_pot, 
                           pot_mod.pot_type, pot_mod.pot_path,
                           disl_center,
                           spec.partition, spec.mem, pot_mod.module_load, pot_mod.appexe, spec.ncpu,
                           **params,)
        
    if sys.argv[1] ==  'plot':
        plot_name = sys.argv[3] 
        params = {}
        if 'global_emin' in dir(spec):
            params['global_emin']           = spec.global_emin
        if 'temp' in dir(spec):
            params['temp']                  = spec.temp
        sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/analysis_tool/ave_atom_pos")))
        import ave_atom
        sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/analysis_tool/nye_ddplot")))
        import ddplot_nye
        for i in range(len(spec.start_temperature)):
            exec(f'import info_{spec.element_struct}_{spec.start_temperature[i]}K as mi')
            ave_file_path = ave_atom.average(mi.latt_const, spec.cutoff_factor, spec.partial_dislocation, plot_name,
                    main_path, spec.pot_element, spec.pot_id, spec.element_struct, spec.config_style, spec.dislocation_type, 
                    spec.simulation_type, unit_cell_size, T=spec.start_temperature[i],
                    ave_file_num=spec.ave_file_num, ave_interval=spec.ave_interval, plotting_region=spec.plotting_region,
                    **params)
            ddplot_nye.main_plot(ave_file_path, spec.pot_element, spec.dislocation_type, spec.thetamax, spec.disl_line_pos,
                  mi.latt_const, spec.cutoff_factor, spec.ddmax_factor, spec.plotting_region,
                  spec.reference, spec.config_style, spec.element_struct, spec.arrowwidth, spec.arrowscale,
                  spec.nye_res, spec.nye_res, plot_name,
                  nye_colorscale=spec.nye_colorscale, test = spec.test)
    
    if sys.argv[1] ==  'post_sd':
        plot_name = sys.argv[3]
        params = {}
        if 'global_emin' in dir(spec):
            params['global_emin']           = spec.global_emin
        if 'temp' in dir(spec):
            params['temp']                  = spec.temp
        sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/analysis_tool/ave_atom_pos")))
        import ave_atom
        sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/analysis_tool/post_analysis")))
        import nye_post
        for i in range(len(spec.start_temperature)):
            exec(f'import info_{spec.element_struct}_{spec.start_temperature[i]}K as mi')
            ave_file_path = ave_atom.ave_file_path(spec.partial_dislocation,plot_name,
                    main_path, spec.pot_element, spec.pot_id, spec.element_struct, spec.config_style, spec.dislocation_type, 
                    spec.simulation_type, unit_cell_size, T=spec.start_temperature[i],
                    **params)
            path = os.path.abspath(os.path.join(ave_file_path, f'../dd_nye_fig_{plot_name}/analysis_raw.txt'))
            nye_post.nye_standard_deviation(spec.start_temperature[i], path, spec.dislocation_type, mi.latt_const, shift=15)
             
