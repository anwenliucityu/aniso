import os
import time

def check_lmp(config, interval):
    '''
    check if lammps is finished every INTERVAL seconds.
    '''
    case_name = config["case_name"]   
    os.chdir('dump_file')
    ovito_dump_file = "ovito " + "dump." + case_name + "_0.cfg"
    if os.path.isfile("dump.final.cfg") == True:
        print("This simulation has done before, plase see results directly")
        os.system(ovito_dump_file)
    else:
        while os.path.isfile("dump.final.cfg") == False:
            time.sleep(interval)
        else:
            os.system(ovito_dump_file)
