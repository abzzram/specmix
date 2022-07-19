import numpy as np
import pdb
import matplotlib as plt

def populate_matrix(specdata, exc_lines, laser_powers, exposure_times):
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
    c_3d = np.zeros((len(exc_lines), len(specdata['filters']), len(FPs)))

    #assembled paired laser lines
    paired_lasers = np.zeros((len(specdata['Lambdas']),len(exc_lines)))
    adjusted_lasers = specdata['lasers'] * laser_powers
    # plt.pyplot.plot(specdata['Lambdas'],specdata['filter_trans'][:,2])
    # plt.pyplot.show()
    # pdb.set_trace()
    #remove duplicate laser lines,
    adjusted_lasers = adjusted_lasers.loc[:,~adjusted_lasers.columns.duplicated()].copy()
    paired_lasers[:,0] = adjusted_lasers.loc[:,exc_lines[0]].sum(1) #add the frist laser pair into one array
    paired_lasers[:,1] = adjusted_lasers.loc[:,exc_lines[1]].sum(1) #add the second laser pair into one array

    #assemble paired filters (sum filter pairs)
    paired_filters =  np.zeros((len(specdata['Lambdas']),len(specdata['filters'])))
    paired_filters[:,0] = specdata['filter_trans'][:,[0,1]].sum(1) #sum first pair
    paired_filters[:,1] = specdata['filter_trans'][:,[2,3]].sum(1) # sum second pair
    cameras = ['BSI','Andor'] 
    # plt.pyplot.plot(specdata['Lambdas'],paired_filters[:,0])
    # plt.pyplot.show()
    # Loop through the entries and populate
    for m, exc_line in enumerate(exc_lines):
        for mp, em_filter in enumerate(specdata['filters']):
            for n, FP in enumerate(FPs):
                # Collect excitation terms
                # laser_spec =  np.array(specdata['lasers'])[:,m]
                laser_spec = np.array(paired_lasers[:,m])
                abs_spec = np.array(specdata['EX_EM'][0,n])
                ex_prod = (laser_spec  *  abs_spec)
                exc_part = np.sum(laser_spec  *  abs_spec) # Collect excitation bits and sum/integrate
                # Collect emission terms
                em_spec = np.array(specdata['EX_EM'][1,n])
                beam_split_trans = np.array(specdata['beam_split'])[:,mp] #not sure which cameras gets which beam yet
                filter_em = np.array(paired_filters)[:,mp]
                Quan_eff = np.array(specdata['QE_cameras'])[:,mp] #make sure this is correct camera order 
                dichroic_trans = np.array(specdata['dichroic'])[:,1]
                #constants
                Quan_yield = np.array(specdata['QY'])[n]
                t_exp = np.array(exposure_times)[m]
                em_part = np.sum(em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp * dichroic_trans)# Collect emission bits and sum/integrate
                em_prod = (em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp * dichroic_trans)
                #plot resulting spectra
            
                # em = (em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp * dichroic_trans)
                # plt.pyplot.plot(specdata['Lambdas'],em)
                # plt.pyplot.show()

                # plt.pyplot.plot(specdata['Lambdas'],laser_spec)
                # # plt.pyplot.plot(specdata['Lambdas'],dichroic_in)
                # plt.pyplot.plot(specdata['Lambdas'],abs_spec)
                # plt.pyplot.show()
                # if n ==3 and m ==1 and mp ==1:
                # plt.pyplot.plot(specdata['Lambdas'],laser_spec)
                # plt.pyplot.plot(specdata['Lambdas'],abs_spec)
                # ex_spec2 = np.array(specdata['EX_EM'][0,1])#
                # plt.pyplot.plot(specdata['Lambdas'],ex_spec2)
                # plt.pyplot.legend(('laser line','excitation spec'),loc='center left', bbox_to_anchor=(0, 1.1), ncol=3)
                # plt.pyplot.xlabel('Wavelength')
                # plt.pyplot.ylabel('Tramsmission') 
                # plt.pyplot.show()                   
                # em_spec2 = np.array(specdata['EX_EM'][1,1])#
                plt.pyplot.plot(specdata['Lambdas'],em_spec)
                plt.pyplot.plot(specdata['Lambdas'],beam_split_trans)
                plt.pyplot.plot(specdata['Lambdas'],filter_em)
                plt.pyplot.plot(specdata['Lambdas'],Quan_eff)
                plt.pyplot.plot(specdata['Lambdas'],dichroic_trans)
                # plt.pyplot.plot(specdata['Lambdas'],em_spec2)
                plt.pyplot.fill_between(specdata['Lambdas'],em_prod)
                # plt.pyplot.plot(specdata['Lambdas'],em_prod)
                plt.pyplot.legend(('em','beam_splitter','filter','QE','dichroic','product'),loc='center left', bbox_to_anchor=(0, 1.1),ncol=3)
                plt.pyplot.xlabel('Wavelength')
                plt.pyplot.ylabel('Tramsmission')
                print(exc_line, em_filter, FP,cameras[mp])      
                plt.pyplot.show() 
                #     pdb.set_trace()         
                c_3d[m,mp,n] = exc_part * em_part

    # Collapse over first two dimensions
    c_2d = np.reshape(c_3d, (c_3d.shape[0]*c_3d.shape[1], c_3d.shape[2]))
    
    return c_2d
