import numpy as np
import pdb
import matplotlib as plt

def populate_matrix(specdata, exc_lines, em_filters, laser_powers, exposure_times):
    # Input: a dictionary of spectral data
    # Input: 
    # Output: a numpy array with our spectral mixing model,

    # Generate empty matrix c_m,m',n
    # For example:
    # 2 exc lines
    # 2 em filters
    # 4 FPS
    # Empty 2 x 2 x 4 array
    FPs = specdata['FPs'] #get FP names
    c_3d = np.zeros((len(exc_lines), len(em_filters), len(FPs)))
    # Loop through the entries and populate
    for m, exc_line in enumerate(exc_lines):
        for mp, em_filter in enumerate(em_filters):
            for n, FP in enumerate(FPs):
                # Collect excitation terms
                laser_spec =  np.array(specdata['lasers'])[:,m]
                # dichroic_in = 1- np.array(specdata['dichroic'])[:,1]
                abs_spec = np.array(specdata['EX_EM'][0,n])
                pdb.set_trace()
                laser_power = laser_powers[m]
                exc_part = np.sum(laser_spec * laser_power *  abs_spec) # Collect excitation bits and sum/integrate
                # Collect emission terms
                em_spec = np.array(specdata['EX_EM'][0,n])
                beam_split_trans = np.array(specdata['beam_split'])[:,mp] #not sure which cameras gets which beam yet
                filter_em = np.array(specdata['filter_trans'])[:,mp]
                Quan_eff = np.array(specdata['QE_cameras'])[:,mp] #make sure this is correct camera order 
                #constants
                Quan_yield = np.array(specdata['QY'])[n]
                t_exp = np.array(exposure_times)[m]
                em_part = np.sum(em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp)# Collect emission bits and sum/integrate
                # # plt.pyplot.plot(specdata['Lambdas'],em_part)
                # # plt.pyplot.show()

                # # plt.pyplot.plot(specdata['Lambdas'],laser_spec)
                # # plt.pyplot.plot(specdata['Lambdas'],dichroic_in)
                # # plt.pyplot.plot(specdata['Lambdas'],abs_spec)
                # # plt.pyplot.show()

                # plt.pyplot.plot(specdata['Lambdas'],em_spec)
                # plt.pyplot.plot(specdata['Lambdas'],beam_split_trans)
                # plt.pyplot.plot(specdata['Lambdas'],filter_em)
                # plt.pyplot.show()       
                # pdb.set_trace()         
                c_3d[m,mp,n] = exc_part * em_part

    # Collapse over first two dimensions
    c_2d = np.reshape(c_3d, (c_3d.shape[0]*c_3d.shape[1], c_3d.shape[2]))
    
    return c_2d
