from copy import deepcopy
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt
import atomman as am
import os
import atomman.unitconvert as uc
import multiprocessing as mp

def fileprocess(path):
    filenames = []
    filenum = 0
    for lists in os.listdir(path):
        sub_path = os.path.join(path, lists)
        if os.path.isfile(sub_path):
            filenum = filenum+1
            new_path = os.path.join(path, lists)
            filenames.append(new_path)
    # sort the list according to its dump sequence
    filenames.sort(key=lambda x: int(x.split('_')[-2][:-5]))
    return filenum, filenames

def plot_and_save(base_path, 
                  disl_path, all_path, dislocation_type, thetamax, disl_line_pos,
                  latt_const, cutoff_factor, b_len, ddmax_factor, plotting_region,
                  reference, config_style, element_struct, arrowwidth, arrowscale, 
                  xbins, ybins, plot_name,
                  nye_colorscale=None,):
    base_system = am.load('atom_data', base_path)
    disl_system = am.load('atom_data',disl_path)

    if config_style=='cylinder':
        disl_system.pbc = [False, False, True]
        base_system.pbc = [False, False, True]
    else:
        disl_system.pbc = [True, True, True]
        base_system.pbc = [True, True, True]

    if element_struct =='hcp':
        alat = latt_const['a']
    else:
        alat = latt_const
    if reference == 0:
        neighbors = base_system.neighborlist(cutoff = cutoff_factor*alat)
        dd = am.defect.DifferentialDisplacement(base_system, disl_system, neighbors=neighbors, reference=0)
    if reference == 1:
        neighbors = disl_system.neighborlist(cutoff = cutoff_factor*alat)
        dd = am.defect.DifferentialDisplacement(disl_system, base_system, neighbors=neighbors, reference=1)

    ddmax = b_len*ddmax_factor   # a/2<110> fcc dislocations use |b|/4
    
    # Set dict of keyword parameter values (just to make settings same for all plots below)
    params = {}
    params['xlim'] = (-plotting_region[0]/2,plotting_region[0]/2)
    params['ylim'] = (-plotting_region[1]/2,plotting_region[1]/2)
    params['figsize'] = 20         # Only one value as the other is chosen to make plots "regular"
    params['arrowwidth'] = arrowwidth    # Made bigger to make arrows easier to see
    params['arrowscale'] = arrowscale#2.4     # Typically chosen to make arrows of length ddmax touch the corresponding atom circles
    params['atomcmap'] = 'gray'
    params['b_length'] = b_len #to dismiss some ddvector if less than 0.1*burers

    if config_style=='quadrupolar':
        params['plotxaxis'] = 'y'
        params['plotyaxis'] = 'z'
        if 'screw' in dislocation_type:
            params['nye_index'] = [0,0]
            ddcomponent = 'x'
        if 'edge' in dislocation_type:
            params['nye_index'] = [0,1]
            ddcomponent = 'y'

    if config_style=='cylinder':
        params['plotxaxis'] = 'x'
        params['plotyaxis'] = 'y'
        if 'screw' in dislocation_type:
            params['nye_index'] = [2,2]
            ddcomponent = 'z'
        if 'edge' in dislocation_type:
            params['nye_index'] = [2,0]
            ddcomponent = 'x'

    ####params for nye###
    params['cmap'] = 'bwr'
    params['xbins'] = xbins
    params['ybins'] = ybins
    params['scale'] = nye_colorscale
    params['disl_line_pos'] = disl_line_pos

    if reference == 0:
        strain = am.defect.Strain(base_system, basesystem=disl_system, cutoff=cutoff_factor*alat, theta_max=thetamax)
    if reference == 1:
        strain = am.defect.Strain(disl_system, basesystem=base_system, cutoff=cutoff_factor*alat, theta_max=thetamax)
    strain.save_to_system()
    straindict = strain.asdict()

    #### plot seetings ##
    par={}
    par['fontsize'] =40
    nye = 'Nye'+'['+str(params['nye_index'][0]+1)+','+str(params['nye_index'][1]+1)+']'
    if reference == 0:
        fig = dd.plot(base_system, ddcomponent, ddmax, **params)
    if reference == 1:
        fig = dd.plot(disl_system, ddcomponent, ddmax, **params)
    #line_index = all_path.index(disl_path)+1
    #plt.title(nye+ f' with DD[{ddcomponent}]   ' ,**par)

    # plot time scale:
    plt.xlabel('x ('+ u'\u212B'+')', fontsize = 45)
    plt.ylabel('y ('+ u'\u212B'+')', fontsize = 45).set_rotation(0)

    path, name = os.path.split(disl_path)
    new_path = os.path.join(path, f'../dd_nye_fig_{plot_name}')
    new_name = name + str('.png')
    new_fig_path = os.path.join(new_path,new_name)
    plt.savefig(new_fig_path)

    line_index = all_path.index(disl_path)+1
    data_path = os.path.join(new_path,'analysis_raw.txt')
    file = open(data_path,'a')
    file.write(f"{line_index:>8d} {fig.x0:>10.8f} {fig.y0:>10.8f} {fig.x_sq:>10.8f} {fig.y_sq:>10.8f} {fig.xy:>10.8f}\n")
    file.close()

def main_plot(ave_file_path, pot_element, dislocation_type, thetamax, disl_line_pos,
                  latt_const, cutoff_factor, ddmax_factor, plotting_region,
                  reference, config_style, element_struct, arrowwidth, arrowscale,
                  xbins, ybins, plot_name,
                  nye_colorscale=None, test = True):
    base_path = os.path.abspath(os.path.join(ave_file_path, '../'+pot_element+'_'+dislocation_type+'_perfect_ref.dat'))
    filenum, disl_paths = fileprocess(ave_file_path)
    b_path = os.path.join(ave_file_path, '../b_length.dat')
    f = open(b_path)
    b_len = float(f.read())

    # open a txt for analysis writing
    path, name = os.path.split(disl_paths[0])
    data_path = os.path.abspath(os.path.join(path,f'../dd_nye_fig_{plot_name}/analysis_raw.txt'))
    if not os.path.exists(os.path.abspath(os.path.join(path,f'../dd_nye_fig_{plot_name}'))):
        os.makedirs(os.path.abspath(os.path.join(path,f'../dd_nye_fig_{plot_name}')))
    file = open(data_path,'w')
    file.write(f"id x0 y0 xx yy xy\n")
    file.close()
    
    if test == True:
        plot_and_save(base_path, disl_paths[0], disl_paths, dislocation_type, thetamax, disl_line_pos,
                  latt_const, cutoff_factor, b_len, ddmax_factor, plotting_region,
                  reference, config_style, element_struct, arrowwidth, arrowscale,
                  xbins, ybins, plot_name,
                  nye_colorscale=nye_colorscale,)
    
    else:
        pool = mp.Pool()
        par = {}
        par['nye_colorscale'] = nye_colorscale
        for i in range(filenum):
            pool.apply_async(plot_and_save, (base_path, disl_paths[i], disl_paths, 
                   dislocation_type, thetamax, disl_line_pos,
                   latt_const, cutoff_factor, b_len, ddmax_factor, plotting_region,
                   reference, config_style, element_struct, arrowwidth, arrowscale,
                   xbins, ybins, plot_name,),
                   par,)
        pool.close()
        pool.join()

if __name__ == '__main__':
    ave_file_path = "/gauss12/home/cityu/anwenliu/scratch/run/core/Ti/hennig_2008_prb/finite_T/hcp/cylinder/screw_a_basal/1100K/100_100_2/ave_file"
    pot_element = 'Ti'
    dislocation_type = 'screw_a_basal'
    thetamax = 10
    disl_line_pos = False
    latt_const = {'a':2.9651961942725 , 'c':4.76664357448851 }
    cutoff_factor =1.3
    ddmax_factor = 1./2
    plotting_region = [60,60]
    reference = 0
    config_style = 'cylinder'
    element_struct = 'hcp'
    arrowwidth = 1/500
    arrowscale = 1.2
    main_plot(ave_file_path, pot_element, dislocation_type, thetamax, disl_line_pos,
                  latt_const, cutoff_factor, ddmax_factor, plotting_region,
                  reference, config_style, element_struct, arrowwidth, arrowscale,
                  nye_colorscale=None, test = False)
