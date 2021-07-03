#!/usr/bin/python3.8
import sys
import os


#########################################
#       information of potential        #
#    DP potential for Ti by Wenqi.T     #
#########################################

# lammps in file data
pot_element    =    'Al'
pot_type       =    'eam'
in_pot         =    ['pair_style         eam/alloy',
                     'pair_coeff         * * Al99.eam.alloy Al', 
                     'neighbor           0.3 bin',
                     'neigh_modify       delay 10']
pot_path       =     [os.path.join(os.path.split(__file__)[0], 'Al99.eam.alloy')]

# potential information
pot_cutoff     =    8.0

# information for slurm
module_load    =    'module load lammps/gcc_openmpi_zen2'
appexe         =    'lmp'

