#!/usr/bin/python3.8 

##########################################################
#    introduce a dislocation with anisotropic solution   #
##########################################################

# information of potential
pot_id             =   'wen_2021_6'
pot_element        =   'Ti'

# information of material
'''
    All dislocation types for hcp metal:
        --- screw <a> basal
        ---  edge <a> basal
        --- screw <a> prismI
        ---  edge <a> prismI
        ---  edge <a> pyrI
        --- screw <c> prismI
        ---  edge <c> prismI
        ---  edge <c> prismII
        --- screw <c+a> pyrII 
        ---  edge <c+a> pyrII
        --- mixed <c+a> prismI
        --- mixed <c+a> pyrI
    All dislocation types for fcc metal:
'''
element_struct     =   'hcp'
dislocation_type    =   'screw_a_basal'
config_shape        =   ['cylinder', 'quadrupolar']
config_style        =   config_shape[1]

# information of simulation box
'''
    If config_style=='cylinder',    z is dislocation line direction, x is slip direction. 
    If config_style=='quadrupolar', x is dislocation line direction, y is slip direction.
'''
num_unit_cell_x     =   1
num_unit_cell_y     =   6
num_unit_cell_z     =   8

# information of dislocation
'''
    For cylinder,    disl_center is the shift according to cylinder center. UNIT: lattice constant
        --- boundary_freeze_width is the width of frozen region, default=3. UNIT: potential cutoff
    For quadrupolar, disl_center is the fraction of box size.
'''
disl_center_x         =   0.25
disl_center_y         =   0.25
boundary_freeze_width =   3.0

# simulation details
'''
    Calc_atomic_stress :   (Default=False)
    Simulation_type:
        Energy_minimization:
            If termal_assist == True: (Default=False)
                --- temp:         the temperature heated to help find global minimized structure. (Default=300K)
                --- cooling_rate: running_steps needed to cool every 1K  (Default=500) 
'''
start_temperature   =   0
simulations         =   ['energy_minimization', 'metastable', 'finite_T']
simulation_type     =   simulations[1]

if simulation_type =='energy_minimization':
    global_emin         =   True
    temp                =   77
    cooling_rate        =   20

if simulation_type =='metastable':
    temp                =   600
    running_steps       =   200000
    dump_interval       =   500

if simulation_type =='finite_T':
    running_steps       =   10000
    dump_interval       =   100
    min_step            =   10

#--------------------------------------
# Define the job specs
#--------------------------------------
ncpu        =   12
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
    calc_aniso_stress:
        ---True:  output a .cfg file with stress field calculated by anisotropic field
        ---False: turn this calculation (Default)
    calc_atomic_stress:
        ---True:  output MD calculated stress field with the dump file.
        ---False: no MD stress info will be recorded (Default)
    output_perfect:
        ---True: output perfect crystal atomic structure (Default)
    ovito:
        ---True: open ovito to see the anisotropic dislocation configuration before MD (Default)
'''
sbatch_job          =   False
calc_aniso_stress   =   False
calc_atomic_stress  =   False
output_perfect      =   True
ovito               =   False
