import numpy as np
import pdb
import matplotlib as plt

def populate_matrix(specdata, exc_lines, laser_powers, exposure_times,**kwargs):
    """ 
    # Input: a dictionary of spectral data
    # Output: a numpy array with our spectral mixing model channel-fluorophore sensitivity constant
        Rows = channels: in the order of Laser pair 1: cam1, cam2, then Laser pair 2: cam1, cam2
        Columns: FPs in the order provided
        Diagonals indicate sensitivity constant of desired FP for each channel. Off-diagonals indicate cross talk 
    # Generate empty matrix c_m,m',n 
     """

    FPs = specdata['FPs'] #get FP names
    cameras = specdata['cameras']
    # c_3d = np.zeros((len(exc_lines), len(specdata['filters']), len(FPs)))

    #assembled paired laser lines
    paired_lasers = np.zeros((len(specdata['Lambdas']),len(exc_lines)))
    # adjusted_lasers = specdata['lasers'] * laser_powers
    adjusted_lasers = specdata['lasers'] * laser_powers / specdata['laser_widths']
    #remove duplicate laser lines,
    adjusted_lasers = adjusted_lasers.loc[:,~adjusted_lasers.columns.duplicated()].copy()
    paired_lasers[:,0] = adjusted_lasers.loc[:,exc_lines[0]].sum(1) #add the frist laser pair into one array
    paired_lasers[:,1] = adjusted_lasers.loc[:,exc_lines[1]].sum(1) #add the second laser pair into one array

    #assemble paired filters 
    #initalize wavelength by filter pair by 2 (# of pairs) matrix to store fitler transmission, and assemble filter pairs
    count = [] #count the numbers of filters in each pair so we can initialize correct matrix
    for filt in specdata['filters']:
        count.append(len(filt)) 
    paired_filters =  np.zeros((len(specdata['Lambdas']),np.max(count),len(specdata['filters']))) #c
    #fill in paired filter matrix 
    j = -1
    for ipair, pairs in enumerate(specdata['filters']):
        for ifilt,filt in enumerate(specdata['filters'][ipair]):
            j=j+1
            paired_filters[:,[ifilt],[ipair]] = specdata['filter_trans'][:,[j]]
    #for plots 
    k=0
    nRows = nCols = len(FPs)
    dayRowCol = np.array([i + 1 for i in range(nRows * nCols)]).reshape(nRows, nCols)
    fig1, ax1 = plt.pyplot.subplots(nrows=nRows,ncols=nCols,figsize=(20,8))
    fig2, ax2 = plt.pyplot.subplots(nrows=nRows,ncols=nCols,figsize=(20,8))
    #initialize c_3d
    c_3d = np.empty((len(exc_lines),np.max(count),len(FPs)))
    c_3d[:] = np.NaN
    for m, exc_line in enumerate(exc_lines):

        for mp, em_filter in enumerate(specdata['filters'][m]):
            for nFP, FP in enumerate(FPs):
                k = k+1 #for subplots
                rowIdx, colIdx = np.argwhere(dayRowCol == k)[0]
                # Collect excitation terms
                laser_spec = np.array(paired_lasers[:,m])
                abs_spec = np.array(specdata['EX_EM'][0,nFP])
                dichroic_ref = np.array(1 - specdata['dichroic'])[:,1] 
                # dichroic_ref =1
                ex_prod = (dichroic_ref * laser_spec  *  abs_spec)
                exc_part = np.sum(dichroic_ref * laser_spec  *  abs_spec) # Collect excitation bits and sum/integrate
                
                #collect emission terms 
                #optional changing of beamsplitter settings, also picks correct camera to use 
                number_excitations = len(exc_line)
                if number_excitations == 1 and 'beamsplitter' in kwargs and len(FPs) > 2: #only urn on single excitation, where correct bs and camera need to be defined
                    # cur_beamsplitter = specdata['beam_split'][:,[0,1]] 
                    curr_bs = specdata['beam_split'][:,kwargs['beamsplitter'][m]] 
                    curr_cam = specdata['QE_cameras'][:,kwargs['beamsplitter'][m]] 
                    beam_split_trans = np.array(curr_bs[:,mp])
                    Quan_eff = np.array(curr_cam[:,mp]) 
                    camera_name = cameras[kwargs['beamsplitter'][m][0]] #will only run when on single excitation, so m will always be 0
                else: #define bs and QE when optional/custom bs is not provided 
                    beam_split_trans = np.array(specdata['beam_split'])[:,mp] 
                    Quan_eff = np.array(specdata['QE_cameras'])[:,mp] 
                    camera_name = cameras[mp]
                # Collect other emission terms
                em_spec = np.array(specdata['EX_EM'][1,nFP])
                filter_em = np.array(paired_filters)[:,mp,m]
                dichroic_trans = np.array(specdata['dichroic'])[:,1] #dichroic always the same 
                #constants
                Quan_yield = np.array(specdata['QY'])[nFP]
                t_exp = np.array(exposure_times)[m]
                em_part = np.sum(em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp * dichroic_trans)# Collect emission bits and sum/integrate
                em_prod = (em_spec * beam_split_trans * filter_em * Quan_eff *  dichroic_trans) #for visualization
                c_3d[m,mp,nFP] = exc_part * em_part
                #plot ex spectra
                axis2 = ax2[rowIdx, colIdx]
                axis2.plot(specdata['Lambdas'],laser_spec)
                axis2.plot(specdata['Lambdas'],abs_spec)
                axis2.plot(specdata['Lambdas'],dichroic_ref)
                axis2.fill_between(specdata['Lambdas'],ex_prod)

                #plot em spectra
                axis1 = ax1[rowIdx, colIdx]
                axis1.plot(specdata['Lambdas'],em_spec)
                axis1.plot(specdata['Lambdas'],beam_split_trans)
                axis1.plot(specdata['Lambdas'],filter_em)
                axis1.plot(specdata['Lambdas'],Quan_eff)
                axis1.plot(specdata['Lambdas'],dichroic_trans)
                axis1.fill_between(specdata['Lambdas'],em_prod)
                #legned 
                if k == len(FPs)**2:
                    axis1.legend(('emission','beam_splitter','filter','QE','dichroic','captured emission'),loc='center left', bbox_to_anchor=(1, 1.1),ncol=1)
                    axis2.legend(('laser pair','excitation','dichroic'),loc='center left', bbox_to_anchor=(1, 1.1),ncol=1)       
                axis1.set_xlabel('Wavelength')
                axis2.set_xlabel('Wavelength')
                axis1.set_ylabel('Tramsmission')
                axis2.set_ylabel('Tramsmission')
                print('Lasers: ',exc_line, em_filter, FP,camera_name, c_3d[m,mp,nFP],'Exposure time: ',t_exp) 
                axis1.set_title(FP + ' emission.' + '\n'  + 'Filter: ' + em_filter + ' Camera: ' + camera_name ,fontsize=8)      
                axis2.set_title(FP + ' excitation. Lasers: ' +  (' | '.join(exc_line) ),fontsize=8)      
                # pdb.set_trace()

    # Collapse over first two dimensions
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.suptitle('Theoretical emission captured by Dragonfly', y = 1) 
    fig2.suptitle('Theoretical excitation spectra/laser lines', y = 1) 
    fig1.savefig('./example_output_em.pdf', bbox_inches=None, dpi=300,facecolor='white')
    fig2.savefig('./example_output_ex.pdf', bbox_inches=None, dpi=300,facecolor='white')
    plt.pyplot.show()
    c_2d = np.reshape(c_3d, (c_3d.shape[0]*c_3d.shape[1], c_3d.shape[2]))
    np.set_printoptions(suppress=True)
    #get rid of nans in umixing matrix, this happends when each camera has different number of laser lines
    if np.isnan(c_2d).any():
        print('...\n...\n Warning: number of excitation lines are different for each camera. Removing NaN rows from unmixing matrix \n...\n...\n')
        c_2d = c_2d[~np.isnan(c_2d).any(axis=1)]
    print(' Channel-fluorophore cross talk matrix: \n',c_2d)
    return c_2d
