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
'''
latt_const     =    {'a': 2.93686668549779, 'c': 4.64745391981822}
mass           =    47.867
c11            =    172.230830688429
c22            =    172.097894399289
c33            =    196.457673987692
c12            =    80.105653142221
c13            =    79.2412242325933
c23            =    79.2074876622586
c44            =    42.6461878238122
c55            =    42.6461968896823
c66            =    46.0780912347327
c14            =    2.37818711711327e-12
c15            =    3.23216040118292e-13
c16            =    8.57858764185267e-13
c24            =    1.06208485061325e-13
c25            =    -7.10666412026032e-14
c26            =    -1.033869215599e-11
c34            =    1.2255448772547e-12
c35            =    -3.62349528343996e-13
c36            =    3.11590357690169e-11
c45            =    2.25188966084107e-14
c46            =    -3.8091578198713e-14
c56            =    7.26116492407895e-12

# initialize elastic constant
elastic_const  = [c11, c12, c13, c14, c15, c16,
                       c22, c23, c24, c25, c26,
                            c33, c34, c35, c36,
                                 c44, c45, c46,
                                      c55, c56,
                                           c66]
