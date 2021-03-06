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
latt_const     =    {'a': 2.93692663105696, 'c': 4.64740824186872}
mass           =    47.867
c11            =    172.177221462474
c22            =    172.039267936055
c33            =    196.545123690037
c12            =    80.9715999584296
c13            =    79.4918427697146
c23            =    79.4568462969497
c44            =    42.3984902856213
c55            =    42.3984992668365
c66            =    45.6209256677673
c14            =    -2.60333370656155e-12
c15            =    7.60986403322819e-13
c16            =    -1.35463522251162e-12
c24            =    2.41155025527145e-12
c25            =    7.49406282364205e-13
c26            =    1.00594499888498e-11
c34            =    1.36426089790656e-12
c35            =    7.7131411311487e-13
c36            =    -3.10554470419626e-11
c45            =    -1.42001292415135e-14
c46            =    3.782847982092e-14
c56            =    -1.57670340351805e-11

# initialize elastic constant
elastic_const  = [c11, c12, c13, c14, c15, c16,
                       c22, c23, c24, c25, c26,
                            c33, c34, c35, c36,
                                 c44, c45, c46,
                                      c55, c56,
                                           c66]
