#!/usr/bin/python3.8
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),".."))) #yes
import default_path as path
spec      = __import__(sys.argv[2].replace('.py',''))
pot_path  = path.pot_path
sys.path.append(os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id)))
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../dislocation_boy/msd")))
import msd
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../potential",spec.pot_element,spec.pot_id)))
import pot_mod
import importlib
import multiprocessing

# Read information in spec file and path
main_path = path.output_path

#####################
#   Main Function   #
#####################

# find potential file, if file is in package then use it, otherwise find outer file
for i in range(len(pot_mod.pot_path)):
    if os.path.exists(pot_mod.pot_path[i]) == False:
        __, pot_name = os.path.split(pot_mod.pot_path[i])
        pot_mod.pot_path[i] = os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id, pot_name))

# Customize
params = {}

if 'sbatch_job' in dir(spec):
    params['sbatch_job']            = spec.sbatch_job
if 'ovito' in dir(spec):
    params['ovito']                 = spec.ovito
if isinstance(spec.start_temperature, int) == True:
            spec.start_temperature = [spec.start_temperature]
unit_cell_size = [spec.num_unit_cell_x, spec.num_unit_cell_y, spec.num_unit_cell_z]

# import the module to constructure the final atomic structure
module_path = os.path.abspath(msd.__file__)

if __name__ == "__main__":
    pool = multiprocessing.Pool()
    if sys.argv[1] == 'calc':
        for i in range(len(spec.start_temperature)):
            exec(f'import info_{spec.element_struct}_{spec.start_temperature[i]}K as mi')
            pool.apply_async(msd.calc_msd,(main_path, spec.pot_element, spec.pot_id, unit_cell_size, spec.element_struct,
                    spec.start_temperature[i], mi.latt_const, mi.mass, pot_mod.pot_path, 
                    spec.warm_runningstep, spec.ave_interval, spec.ave_times, pot_mod.in_pot,
                    spec.partition, spec.mem, pot_mod.module_load, pot_mod.appexe, spec.ncpu, spec.dump,),
                    params,)
        pool.close()
        pool.join()

