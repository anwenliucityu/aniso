#!/usr/bin/python3.8 

##########################################################
#  introduce  dislocation(s) with anisotropic solution   #
##########################################################

# information of potential
pot_id             =   'zhou_2004_prb'
pot_element        =   'W'

# information of material
'''
    All dislocation types for hcp metal:
        --- screw <a> basal                   [100,100,1], cutoff=1.3, thetamax=10, nye_colorscale=0.08
        ---  edge <a> basal
        --- screw <a> prismI
        ---  edge <a> prismI
        ---  edge <a> pyrI
        --- screw <c> prismI
        ---  edge <c> prismI
        ---  edge <c> prismII
        --- screw <c+a> pyrII (full+partial)  [90,180,1]
        --- screw <c+a> PyrI (in testing)
        ---  edge <c+a> pyrII
        --- mixed <c+a> prismI
        --- mixed <c+a> pyrI
    All dislocation types for fcc metal:
    
    Al dislocation type for bcc metal:
        --- screw <111>            [60,100,1], 7722 atoms
'''
element_struct      =   'bcc'
dislocation_type    =   'screw'
partial_dislocation =   False
partial_split_dis   =   6
config_shape        =   ['cylinder', 'quadrupolar']
config_style        =   config_shape[0]

# information of simulation box
'''
    If config_style=='cylinder',    z is dislocation line direction, x is slip direction.
    If config_style=='quadrupolar', x is dislocation line direction, y is slip direction.
'''
num_unit_cell_x     =   60
num_unit_cell_y     =   100
num_unit_cell_z     =   1

# information of dislocation
'''
    For cylinder,    disl_center is the shift according to cylinder center. UNIT: lattice constant
        --- boundary_freeze_width is the width of frozen region, default=3. UNIT: potential cutoff
    For quadrupolar, disl_center is the fraction of box size.
'''
disl_center_x         =   0.25
disl_center_y         =   0.25
boundary_freeze_width =   2.0

# simulation details
'''
    Calc_atomic_stress :   (Default=False)
    Simulation_type:
        Energy_minimization:
            If termal_assist == True: (Default=False)
                --- temp:         the temperature heated to help find global minimized structure. (Default=300K)
                --- cooling_rate: running_steps needed to cool every 1K  (Default=500)
'''
start_temperature   =   0#[100,200,300,400,500,600,700,800,900,1000,1100,1200]#[1200,1300,1400,1500,1600]#[500,600,700,800,900,1000,1100,1200]
simulations         =   ['energy_minimization', 'metastable', 'finite_T']
simulation_type     =   simulations[1]

if simulation_type =='energy_minimization':
    global_emin         =   False
    temp                =   300
    cooling_rate        =   500

if simulation_type =='metastable':
    temp                =   200
    running_steps       =   50000
    dump_interval       =   500

if simulation_type =='finite_T':
    # min step is used for visualize instant core structure with fewer average
    running_steps       =   50000
    dump_interval       =   500
    min_step            =   10

#--------------------------------------
# Define the job specs
#--------------------------------------
ncpu        =   32
partition   =   'xlong'
project     =   'default'
mem         =   '60G'

#--------------------------------------
# customize DDplot and Nye tensor
#--------------------------------------
'''
    core_region: the region near core used to plot
    reference:   0:perfect frame   1: dislocated frame
'''
ave_file_num        =   1
ave_interval        =   1
ddmax_factor        =   1./2
cutoff_factor       =   1.3
reference           =   0  

# plotting settings
arrowwidth          = 1/100
arrowscale          = 1.2
plotting_region     = [20,20]
nye_colorscale      = 0.08
nye_res             = 200
thetamax            = 10
test                = False
disl_line_pos       = True

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
        ---True: open ovito to see the initial configuration before MD (Default)
'''
sbatch_job          =   False
calc_aniso_stress   =   False
calc_atomic_stress  =   False
output_perfect      =   True
ovito               =   True
