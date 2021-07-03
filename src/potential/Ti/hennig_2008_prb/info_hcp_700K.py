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
latt_const     =    {'a': 2.95242356133086, 'c': 4.7380249642853}
mass           =    47.867
c11            =    131.530103185252
c22            =    129.116191454514 
c33            =    140.187208631292
c12            =    88.4659438281159
c13            =    54.513876340782
c23            =    52.2921164282582
c44            =    31.5074611350635
c55            =    31.2843815212149
c66            =    18.7445864958192
c14            =    0.336900500595694
c15            =    1.90809866798563
c16            =    0.301539017910008
c24            =    0.750458743018045
c25            =    -2.09960287331169
c26            =    1.69979750792887
c34            =    -0.145296073395725
c35            =    0.329497459469322
c36            =    0.822383880218219
c45            =    0.60960669181458
c46            =    -1.10842091718985
c56            =    -1.60945832952132

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
msd            = 0.136579

