import sys 
import os
sys.path.append(os.path.abspath(os.path.join(__file__,'../../disl_aniso')))
import input_dict
from perfect_cryst_constructor import perfect_cryst_constructor, perfect_hcp_constructor
sys.path.append(os.path.abspath(os.path.join(__file__,'../../write_output')))
from write_output import write_datafile, write_job
config_style     = 'normal'
from write_in import in_file

def calc_msd(main_path, pot_element, pot_name, unit_cell_size, 
            element_struct, start_temperature, latt_const, mass,  pot_path,
            warm_runningstep, ave_interval, ave_times, in_pot,
            partition, mem, module_load, appexe, ncpu, dump,
            ovito=False, sbatch_job=False):
    
    if element_struct =='fcc':
        dislocation_type = 'screw'
    if element_struct =='hcp':
        dislocation_type = 'screw_a_basal'
    if element_struct =='bcc':
        dislocation_type = 'normal'  
        
    local = locals()
    exec(f'config = input_dict.{element_struct}_{dislocation_type}')
    config = local['config']
    
    # Creat a new directory path for installing output file
    box  = str(unit_cell_size[0])+'_'+str(unit_cell_size[1])+'_'+ str(unit_cell_size[2])
    T = str(int(start_temperature))+ 'K'
    dir = os.path.join(main_path, pot_element, pot_name, 'msd', element_struct, T, box)
    if not os.path.isdir(dir):
        os.makedirs(dir)
    directory = dir
    dump_dir = os.path.join(dir, 'dump_file')
    if not os.path.isdir(dump_dir):
        os.makedirs(dump_dir)
    print("Atomic structure is prepared ......")
    print("work_path = " + dir)

    # construct a perfect crystal
    if element_struct == 'hcp':
        atom_coor, atom_type, box_boundary, new_box_lattice_const, frame_new, repeat_para = \
                perfect_hcp_constructor(config, config_style, latt_const, pot_element, 
                        dislocation_type, unit_cell_size)
    if element_struct == 'fcc' or element_struct == 'bcc':
        atom_coor, atom_type, box_boundary, repeat_para = perfect_cryst_constructor(config, config_style,
                element_struct, dislocation_type, latt_const, unit_cell_size,)

    # write perfect configuration for dislocation analysis
    write_datafile(pot_element, dislocation_type, mass, directory,
              atom_coor, atom_type, box_boundary, suffix='_perfect_ref.dat')
    
    # write in file to calculate msd
    in_file(directory, pot_element, dislocation_type, pot_path, pot_name, mass, in_pot,
             start_temperature, warm_runningstep, ave_interval, ave_times, dump)
    
    # write job_file
    write_job(pot_element, dislocation_type, partition, mem, module_load, appexe, directory, ncpu=ncpu)
    
    # ovito
    if ovito == True:
        init_path = os.getcwd()
        os.chdir(directory)
        os.system(f"ovito {pot_element}_{dislocation_type}_perfect_ref.dat")  
        os.chdir(init_path)
      
    # cd job.sh directory to submit job
    if sbatch_job == True:
        init_path = os.getcwd()
        os.chdir(directory)
        os.system('sbatch job.sh')
        os.chdir(init_path)
    
