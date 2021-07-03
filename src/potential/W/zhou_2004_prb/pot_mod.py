#!/usr/bin/python3.8
import sys
import os


#########################################
#       information of potential        #
#    DP potential for Ti by Wenqi.T     #
#########################################

# lammps in file data
pot_element    =    ''
pot_type       =    'eam'
in_pot         =    ['pair_style         eam/alloy',
                     'pair_coeff         * *  W_Zhou04.eam.alloy W', 
                     'neighbor           0.3 bin',
                     'neigh_modify       delay 10']
pot_path       =     [os.path.join(os.path.split(__file__)[0], 'W_Zhou04.eam.alloy')]

# potential information
pot_cutoff     =    8.0

# information for slurm
module_load    =    'module load lammps/20210210/gcc_openmpi_znver2_openkim'
appexe         =    'lmp'

