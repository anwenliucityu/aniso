#!/usr/bin/python3.8 

#########################################
#      information of DP potential      #
#########################################

# information of material
'''
    units:
        temperature      : Kelvin
        lengh            : Angstrom
        mass             : g/cm^3
        elastic constant : GPa
        * msd            :
'''
latt_const     =    {'a':2.94918519344792, 'c': 4.72994758498517}
mass           =    47.867
c11            =    129.682472025899 
c22            =    133.777014512807
c33            =    144.645181954884
c12            =    91.1444311324859
c13            =    55.9095665601803
c23            =    55.3539997734204
c44            =    33.3776932523097
c55            =    33.5500643692152
c66            =    21.1876659735957
c14            =    -0.11013140707936
c15            =    0.030056536319783
c16            =    -0.417375766112135
c24            =    0.32705818203068
c25            =    0.394494207032879
c26            =    1.93469165705306
c34            =    1.09676737819415
c35            =    -0.224786620564975
c36            =    0.0447038751605501
c45            =    -0.365418717338909
c46            =    -0.915204946662133
c56            =    -0.214455751547776

# initialize elastic constant
elastic_const  = [c11, c12, c13, c14, c15, c16,
                       c22, c23, c24, c25, c26,
                            c33, c34, c35, c36,
                                 c44, c45, c46,
                                      c55, c56,
                                           c66]
                                           
###########################################
#            additional info              #
###########################################
msd            = 0.120911

