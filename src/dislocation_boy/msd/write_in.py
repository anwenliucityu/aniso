import os

def in_file(directory, pot_element, dislocation_type, pot_path, pot_name, mass, in_pot,
              start_temp, warm_runningstep, ave_interval, ave_times, dump):
    case_name = pot_element + '_' + dislocation_type
    ini_path = os.getcwd()
    os.chdir(directory)
    for i in range(len(pot_path)):
        pot_name = os.path.basename(pot_path[i])
        os.system(f'ln -s {pot_path[i]} {pot_name}')
    os.chdir(ini_path)
    filename = 'in.lammps'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
    file = open(file_path_name, 'w')
    file.write("#----------------------------------------#\n"
               "#---   basic info & initialization    ---#\n"
               "#----------------------------------------#\n")
    file.write(f"clear\n"
                "units         metal\n"
                "dimension     3\n"
                "boundary      p p p\n"
                "atom_style    atomic\n"
                "atom_modify   map array\n"
               f"read_data     {case_name}_perfect_ref.dat\n"
               f"mass          1 {mass}\n\n")   

    file.write("#----------------------------------------#\n"
               "#---            potential             ---#\n"
               "#----------------------------------------#\n")

    for i in range(len(in_pot)):
        file.write(f"{in_pot[i]}\n")
    file.write('\n')
    
    start_temp = float(start_temp)
    file.write("#----------------------------------------#\n"
               "#              initialize                #\n"
               "#----------------------------------------#\n")
    file.write(f"fix       sys all npt temp {start_temp} {start_temp} 0.1 aniso 1.0 1.0 1.0\n"
               f"velocity  all create {start_temp} 32544 rot yes dist gaussian\n"
               f"run       {warm_runningstep}\n\n")

    total = int(ave_interval*ave_times)
    file.write("#----------------------------------------#\n"
               "#              compute                   #\n"
               "#----------------------------------------#\n")
    file.write("compute     frame all msd com yes\n"
               "variable    msd equal c_frame[4]\n"
              f"fix         average all ave/time {ave_interval} {ave_times} {total} v_msd file msd_info.txt\n\n")

    file.write("#---------------------------------------#\n"
               "#            thermo and dump            #\n"
               "#---------------------------------------#\n")
    total = int(total*3)
    file.write("thermo_style       custom step temp pe press lx ly lz\n"
               "thermo             100\n")
    if dump == True:
        file.write("shell              cd dump_file\n"
                   "dump               config all custom 2000 dump.relax_* type mass x y z\n")               
    file.write(f"run                {total}\n"
                "print"' "'"All done"'"')


