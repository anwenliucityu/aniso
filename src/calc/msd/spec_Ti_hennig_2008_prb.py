#!/usr/bin/python3.8 

##########################################################
#       calc meam square displacement at finite T        #
##########################################################

# information of potential
pot_id             =   'hennig_2008_prb'
pot_element        =   'Ti'
element_struct     =   'bcc'

# information of simulation box
num_unit_cell_x     =   5
num_unit_cell_y     =   5
num_unit_cell_z     =   5

# simulation details
start_temperature   =   [800]#[1200,1300,1400,1500,1600]#[100,200,300,400,500,600,700,800,900,1000,1100,1200]
warm_runningstep    =   10000
ave_interval        =   5
ave_times           =   1000
dump                =   True

# dump info

#--------------------------------------
# Define the job specs
#--------------------------------------
ncpu        =   8
partition   =   'xlong'
project     =   'default'
mem         =   '48G'

#---------------------------------------
# Additional function (turn TRUE if needed)
#---------------------------------------
'''
    sbatch_job:
        ---True:  job will be sbatched  
        ---False: job wont be sbatched  (Default)
    ovito:
        ---True: open ovito to see the anisotropic dislocation configuration before MD (Default)
'''
sbatch_job          =   True
ovito               =   False
