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
latt_const     =    {'a': 2.93312780358881, 'c': 4.6546126928194}
mass           =    47.867
c11            =    165.107922174266
c22            =    164.955110093071
c33            =    193.346751509692
c12            =    86.5983538788439
c13            =    80.9057873357119
c23            =    80.8696795724422
c44            =    41.7911249463292
c55            =    41.7911309378191
c66            =    39.2742041217947
c14            =    5.89655618532365e-13
c15            =    -3.68414578128537e-13
c16            =    1.32713984956229e-12
c24            =    6.82711475796849e-13
c25            =    1.09113766382858e-13
c26            =    8.6586815890581e-13
c34            =    -1.72238810841595e-13
c35            =    -5.09450390751519e-13
c36            =    7.08697789502659e-13
c45            =    6.24352821451245e-14
c46            =    -2.44444877715413e-14
c56            =    -1.04910163983448e-12

# initialize elastic constant
elastic_const  = [c11, c12, c13, c14, c15, c16,
                       c22, c23, c24, c25, c26,
                            c33, c34, c35, c36,
                                 c44, c45, c46,
                                      c55, c56,
                                           c66]