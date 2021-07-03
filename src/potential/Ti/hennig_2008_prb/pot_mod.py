#!/usr/bin/python3.8
import sys
import os

#########################################
#       information of potential        #
#    DP potential for Ti by Wenqi.T     #
#########################################

# lammps in file data
pot_element    =    'meam/spline'
pot_type       =    'meam'
in_pot         =    ['pair_style         meam/spline',
                     'pair_coeff         * * Ti.meam.spline Ti', 
                     'neighbor           0.3 bin']
pot_path       =     [os.path.join(os.path.split(__file__)[0], 'Ti.meam.spline')]

# potential information
pot_cutoff     =    8.4

# information for slurm
module_load    =    'module load lammps/20210210/gcc_openmpi_znver2_openkim'
appexe         =    'lmp'

