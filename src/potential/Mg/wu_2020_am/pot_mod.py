#!/usr/bin/python3.8
import sys
import os

#########################################
#       information of potential        #
#    DP potential for Ti by Wenqi.T     #
#########################################

# lammps in file data
pot_element    =    'Mg'
pot_type       =    'meam'
in_pot         =    ['pair_style         meam/c',
                     'pair_coeff         * *  library.meam Mg Mg.meam Mg', 
                     'neighbor           0.3 bin']
pot_path       =     [os.path.join(os.path.split(__file__)[0], 'library.meam'),
                      os.path.join(os.path.split(__file__)[0], 'Mg.meam')]

# potential information
pot_cutoff     =    8.4

# information for slurm
module_load    =    'lammps/20201029/aocc_openmpi_znver2'
appexe         =    'lmp'

