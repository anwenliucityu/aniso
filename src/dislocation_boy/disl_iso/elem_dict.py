'''list of properties of elements'''

elem_dict = {
    "EAM_Mishin_Al" : {
        "elem_symbol"       : 'Al',
        "pot_file_name"     : 'Al99.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 4.05,      # lattice constant (Angstrom)
        "mass"              : 26.98,
        "pot_cutoff"        : 6.287,     # cut-off radius of potential (Angstrom)
        "C11"               : 1.14e11,   # (Pa)
        "C12"               : 0.619e11,
        "C44"               : 0.316e11,
        "shear_modulus"     : 0.316e11,  # C44
        "poisson_ratio"     : 0.331016,  # C12/(2 * (C44 + C12))

        "elastic_const_initia"  :  [0.619e11+2*0.316e11, 0.619e11, 0.619e11,   0,   0,   0,
                                               0.619e11+2*0.316e11, 0.619e11,   0,   0,   0,
                                                         0.619e11+2*0.316e11,   0,   0,   0,
                                                                         0.316e11,   0,   0,
                                                                              0.316e11,   0,
                                                                                   0.316e11],

        "elastic_const_initial"  :  [1.14e11, 0.619e11, 0.619e11,   0,   0,   0,
                                               1.14e11, 0.619e11,   0,   0,   0,
                                                         1.14e11,   0,   0,   0,
                                                             0.316e11,   0,   0,
                                                                  0.316e11,   0,
                                                                       0.316e11],
    },

    "EAM_Mishin_Al_300K" : {
        "elem_symbol"       : 'Al',
        "pot_file_name"     : 'Al99.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 4.065,      # lattice constant (Angstrom)
        "mass"              : 26.98,
        "pot_cutoff"        : 6.287,     # cut-off radius of potential (Angstrom)
        "C11"               : 1.10e11,   # (Pa)
        "C12"               : 0.643e11,
        "C44"               : 0.337e11,
        "shear_modulus"     : 0.337e11,  # C44
        "poisson_ratio"     : 0.328061,  # C12/(2 * (C44 + C12))

        "elastic_const_initial"  :  [1.10e11, 0.643e11, 0.643e11,   0,   0,   0,
                                               1.10e11, 0.643e11,   0,   0,   0,
                                                         1.10e11,   0,   0,   0,
                                                             0.337e11,   0,   0,
                                                                  0.337e11,   0,
                                                                       0.337e11],
    },

    "EAM1_Mishin_Cu" : {
        "elem_symbol"       : 'Cu',
        "pot_file_name"     : 'Cu01.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 3.615,     # lattice constant (Angstrom)
        "mass"              : 63.546,
        "pot_cutoff"        : 5.50679,   # cut-off radius of potential (Angstrom)
        "C11"               : 1.699e11,  # (Pa)
        "C12"               : 1.226e11,
        "C44"               : 0.762e11,
        "shear_modulus"     : 0.762e11,  # C44
        "poisson_ratio"     : 0.3083501, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [1.699e11, 1.226e11, 1.226e11,   0,   0,   0,
                                               1.699e11, 1.226e11,   0,   0,   0,
                                                         1.699e11,   0,   0,   0,
                                                              0.762e11,   0,   0,
                                                                   0.762e11,   0,
                                                                        0.762e11],

    },
        "EAM_Mishin_Ni" : {
        "elem_symbol"       : 'Ni',
        "pot_file_name"     : 'Ni99.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 3.52,     # lattice constant (Angstrom)
        "mass"              : 58.71,
        "pot_cutoff"        : 5.804,   # cut-off radius of potential (Angstrom)
        "C11"               : 2.47e11,  # (Pa)
        "C12"               : 5.19e11,
        "C44"               : 1.59e11,
        "shear_modulus"     : 1.25e11,  # C44
        "poisson_ratio"     : 0.27106227106227104, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [2.47e11, 1.48e11, 1.48e11,   0,   0,   0,
                                               2.47e11, 1.48e11,   0,   0,   0,
                                                         2.47e11,   0,   0,   0,
                                                              1.25e11,   0,   0,
                                                                   1.25e11,   0,
                                                                        1.25e11],

    },
        "EAM_Ma_W_300K" : {
        "elem_symbol"       : 'W',
        "pot_file_name"     : 'W_eam3.fs',
        "pair_style"        : 'eam/fs',
        "cryst_structure"   : 'bcc',     # crystal structure
        "latt_const"        : 3.202639,     # lattice constant (Angstrom)
        "mass"              : 183.84,
        "pot_cutoff"        : 5.5,   # cut-off radius of potential (Angstrom)
        "C11"               : 2.47e11,  # (Pa)
        "C12"               : 1.48e11,
        "C44"               : 1.25e11,
        "shear_modulus"     : 1.25e11,  # C44
        "poisson_ratio"     : 0.2748091603053435, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [2.47e11, 1.48e11, 1.48e11,   0,   0,   0,
                                               2.47e11, 1.48e11,   0,   0,   0,
                                                         2.47e11,   0,   0,   0,
                                                              1.25e11,   0,   0,
                                                                   1.25e11,   0,
                                                                        1.25e11],
    },
        "EAM_Mendelev_Fe_300K" : {
        "elem_symbol"       : 'Fe',
        "pot_file_name"     : 'Fe_2.eam.fs',
        "pair_style"        : 'eam/fs',
        "cryst_structure"   : 'bcc',     # crystal structure
        "latt_const"        : 2.85887562574632,     # lattice constant (Angstrom)
        "mass"              : 58.845,
        "pot_cutoff"        : 5.5,   # cut-off radius of potential (Angstrom)
        "C11"               : 2.43e11,  # (Pa)
        "C12"               : 1.45e11,
        "C44"               : 1.77e11,
        "shear_modulus"     : 1.77e11,  # C44
        "poisson_ratio"     : 0.2251552795031056, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [2.43e11, 1.45e11, 1.45e11,   0,   0,   0,
                                               2.43e11, 1.45e11,   0,   0,   0,
                                                         2.43e11,   0,   0,   0,
                                                              1.16e11,   0,   0,
                                                                   1.16e11,   0,
                                                                        1.16e11],

     },
        "Zhang_Nb_EAM" : {
        "elem_symbol"       : 'Nb',
        "pot_file_name"     : 'Nb.eam.fs',
        "pair_style"        : 'eam/fs',
        "cryst_structure"   : 'bcc',     # crystal structure
        "latt_const"        : 3.3023736761276,     # lattice constant (Angstrom)
        "mass"              : 90,
        "pot_cutoff"        : 5.5,   # cut-off radius of potential (Angstrom)  
        "C12"               : 1.32699709292376e11,  # (Pa)
        "C44"               : 0.923768303994206e11,
     },

        "Ravelo_Ta_EAM" : {
        "elem_symbol"       : 'Ta',
        "pot_file_name"     : 'Ta_Ravelo.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'bcc',     # crystal structure
        "latt_const"        : 3.30399999955882,# lattice constant (Angstrom)
        "mass"              : 180.95,
        "pot_cutoff"        : 5.5,   # cut-off radius of potential (Angstrom)  
        "C12"               : 1.60699974286487e11,  # (Pa)
        "C44"               : 0.818044529639637e11,
     },

        "Zhou_Pt_EAM" : {
        "elem_symbol"       : 'Zhou',
        "pot_file_name"     : 'Pt_Zhou04.eam.alloy',
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 3.92003421444833,# lattice constant (Angstrom)
        "mass"              : 195,
        "pot_cutoff"        : 5.5,   # cut-off radius of potential           (Angstrom)  
        "C12"               : 2.50462387622492e11,  # (Pa)
        "C44"               : 0.761074990718758e11,
     },
     
}
