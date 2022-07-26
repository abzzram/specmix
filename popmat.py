import numpy as np
import pdb
import matplotlib as plt

def populate_matrix(specdata, exc_lines, laser_powers, exposure_times,**kwargs):
    # Input: a dictionary of spectral data
    # Output: a numpy array with our spectral mixing model,
    #optional input:
        #'beamsplitter = 'none'
            #this is for when no beamsplitter was used 
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
    # adjusted_lasers = specdata['lasers'] * laser_powers
    adjusted_lasers = specdata['lasers'] * laser_powers / specdata['laser_widths']
    #remove duplicate laser lines,
    adjusted_lasers = adjusted_lasers.loc[:,~adjusted_lasers.columns.duplicated()].copy()
    paired_lasers[:,0] = adjusted_lasers.loc[:,exc_lines[0]].sum(1) #add the frist laser pair into one array
    paired_lasers[:,1] = adjusted_lasers.loc[:,exc_lines[1]].sum(1) #add the second laser pair into one array

    #assemble paired filters 
    #initalize wavelength by filter pair by 2 (# of pairs) matrix to store fitler transmission, and assemble filter pairs
    paired_filters =  np.zeros((len(specdata['Lambdas']),len(specdata['filters']),len(specdata['filters'])))
    #is there a better way to do the following? 
    paired_filters[:,[0,1],[0]] = specdata['filter_trans'][:,[0,1]]
    paired_filters[:,[0,1],[1]] = specdata['filter_trans'][:,[2,3]]
    cameras = specdata['cameras']
    #if beam splitter is set to none, change values to 1 
    if kwargs['beamsplitter'] == 'none':
        specdata['beam_split'][:,:] = 1
    k=0
    nRows = nCols = len(FPs)
    dayRowCol = np.array([i + 1 for i in range(nRows * nCols)]).reshape(nRows, nCols)
    fig1, ax1 = plt.pyplot.subplots(nrows=nRows,ncols=nCols,figsize=(20,8))
    fig2, ax2 = plt.pyplot.subplots(nrows=nRows,ncols=nCols,figsize=(20,8))
    for m, exc_line in enumerate(exc_lines):
        for mp, em_filter in enumerate(specdata['filters']):
            for n, FP in enumerate(FPs):
                k = k+1 #for subplots
                rowIdx, colIdx = np.argwhere(dayRowCol == k)[0]
                # Collect excitation terms
                laser_spec = np.array(paired_lasers[:,m])
                abs_spec = np.array(specdata['EX_EM'][0,n])
                ex_prod = (laser_spec  *  abs_spec)
                exc_part = np.sum(laser_spec  *  abs_spec) # Collect excitation bits and sum/integrate
                # Collect emission terms
                em_spec = np.array(specdata['EX_EM'][1,n])
                beam_split_trans = np.array(specdata['beam_split'])[:,mp] #not sure which cameras gets which beam yet
                filter_em = np.array(paired_filters)[:,mp,m]
                Quan_eff = np.array(specdata['QE_cameras'])[:,mp] #make sure this is correct camera order 
                dichroic_trans = np.array(specdata['dichroic'])[:,1]
                #constants
                Quan_yield = np.array(specdata['QY'])[n]
                t_exp = np.array(exposure_times)[m]
                em_part = np.sum(em_spec * beam_split_trans * filter_em * Quan_eff * Quan_yield * t_exp * dichroic_trans)# Collect emission bits and sum/integrate
                em_prod = (em_spec * beam_split_trans * filter_em * Quan_eff *  dichroic_trans) #for visualization
                c_3d[m,mp,n] = exc_part * em_part

                #plot ex spectra
                axis2 = ax2[rowIdx, colIdx]
                axis2.plot(specdata['Lambdas'],laser_spec)
                axis2.plot(specdata['Lambdas'],abs_spec)
                axis2.fill_between(specdata['Lambdas'],ex_prod)

                #plot em spectra
                axis1 = ax1[rowIdx, colIdx]
                axis1.plot(specdata['Lambdas'],em_spec)
                axis1.plot(specdata['Lambdas'],beam_split_trans)
                axis1.plot(specdata['Lambdas'],filter_em)
                axis1.plot(specdata['Lambdas'],Quan_eff)
                axis1.plot(specdata['Lambdas'],dichroic_trans)
                axis1.fill_between(specdata['Lambdas'],em_prod)
                if k ==4:
                    axis1.legend(('em','beam_splitter','filter','QE','dichroic','captured emission'),loc='center left', bbox_to_anchor=(1, 1.1),ncol=3)
                    axis2.legend(('ex','laser pair'),loc='center left', bbox_to_anchor=(1, 1.1),ncol=3)                  
                if k>=12:
                    axis1.set_xlabel('Wavelength')
                    axis2.set_xlabel('Wavelength')
                if k%4 ==1:
                    axis1.set_ylabel('Tramsmission')
                    axis2.set_ylabel('Tramsmission')
                print(exc_line, em_filter, FP,cameras[mp], c_3d[m,mp,n])   
                axis1.set_title(FP + ' emission.' + '\n'  + 'Filter: ' + em_filter[m] + ' Camera: ' + cameras[mp] ,fontsize=8)      
                axis2.set_title(FP + ' excitation. Lasers: ' +  (' | '.join(exc_line) ),fontsize=8)      

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
    print('Unmixing matrix: \n',c_2d)
    return c_2d
