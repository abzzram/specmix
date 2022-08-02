#Functions for spectral modeling 
#import stuff
import json
import pandas as pd
import csv
import numpy as np
import ast
import wget
import json
import array 
from numpy import nan
# %matplotlib inline
import matplotlib.pyplot as plt
import pdb
import os
from os.path import exists

#make it all into a function for retreiving SP spectra 
#returns spectra in EX_EM (2 (ex/em) by n (#of FPs) by spectra (arranged from wavelength 300 to 800))
#also return quantum yield of FPs in QY (N by 1)
#also returns 'lambas', which is the wavelegnths that EX/EM corresponds to (300 to 800)
def get_FP_spectra(FPs):
#download all FP data from FP base, modify Talon's code 
    # filename = wget.download('https://www.fpbase.org/api/proteins/spectra/?name__iexact=mTagBFP2&default_state__qy__gte=0.7&format=json')
    FP_file = 'FP_spectra.wget'
    file_exists = os.path.exists(FP_file) #check if file exists 
    if not file_exists:
        print('downloading spectra...')
        filename = wget.download('https://www.fpbase.org/api/proteins/spectra/?name__iexact=mTagBFP2&default_state__qy__gte=0.7&format=json', out = 'FP_spectra.wget')
        # url = 'https://www.fpbase.org/api/proteins/spectra/?name__iexact=mTagBFP2&default_state__qy__gte=0.7&format=json'
        # !wget 
    else:
        print('spectra file found...')
        filename = './FP_spectra.wget'

    with open(filename, 'r') as f:
        fc = json.load(f)
        f.close()
#extract names of flourescent proteins 
    FP_names = []
    for i in fc:
        iname = (i['name'])
        FP_names.append(iname)       
#get indexes of desired proteins 
    FP_inds = []
    for fp in FPs:
        FP_inds.append([FP_names.index(name) for name in FP_names if name == fp])
    FP_inds = np.array(FP_inds) #convert to array
    
    #get spectra for each FP, 
    lambdas = np.arange(300,801,1) #range wavelengths to store
    EX_EM = np.zeros((2,len(FPs),len(lambdas)))
    QY = np.zeros(len(FPs)); #quantum yield
    for i,ifp in enumerate(FP_inds):
        lambda_ex = np.array(fc[int(ifp)]['spectra'][0]['data']) # n by 2 array of wavelengths and data
        for n, ilambda in enumerate(lambda_ex):
            if (lambdas == ilambda[0]).any(): #skip wavelengths that aren't shared 
                jj = np.argwhere(lambdas == ilambda[0])[0][0] #find index of ex wavelength containedin databse
                EX_EM[0,i,jj]  = lambda_ex[n,1]
        lambda_em = np.array(fc[int(ifp)]['spectra'][1]['data']) # n by 2 array of em wavelengths and data
        for n2, ilambda2 in enumerate(lambda_em):
            if (lambdas == ilambda2[0]).any(): #skip wavelengths that aren't shared 
                jj = np.argwhere(lambdas == ilambda2[0])[0][0] #find index of ex wavelength containedin databse
                EX_EM[1,i,jj]  = lambda_em[n2,1]
    #get quantum yield
        QY[i] = fc[int(ifp)]['spectra'][1]['qy']
    # os.remove('FPs.wget') #delete the downloaded FP file 
    return(EX_EM,QY, lambdas)


#function for retriving camera QEs
#put it all into a function
def get_QEs(Lambdas, cam1path, cam2path,cameras):
# def get_QEs(Lambdas, bsi_path, ixon_path):
    #Retreive quantum efficiency of the BSI Prime Express and iXon cameras
    #inputs: wavelengths corresponding to QE needed (should be the same as all other spectra data )
    #output, wavelength by 2 array of QEs for BSI (col 0) and iXon (col 1)
        #for intepolating nans ('https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array')
        #other inputs are the file locations for each camera's QE data 
        #Lets define first a simple helper function in order to make it more straightforward to handle indices and logical indices of NaNs:
    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]
    
    #retrive camera quantum efficiency data
    # BSI = pd.read_csv(bsi_path)
    # iXon = pd.read_csv(ixon_path)
    cam1 = pd.read_csv(cam1path)
    cam2 = pd.read_csv(cam2path)    
    #round the QE wavelengths for BSI data. (names of columns are as they were when downloaded)
    if cameras[0] == 'BSI_Prime_Express':
        cam1.BSI_New = cam1.BSI_New.round(0)
        #turn into percetages 
        cam1.QE = cam1.QE / 100
    if cameras[1] == 'BSI_Prime_Express':
        cam2.BSI_New = cam2.BSI_New.round(0)
        #turn into percetages 
        cam2.QE = cam2.QE / 100
    if cameras[0] == 'Andor_iXon':
        cam1['IXON-L-888 Sensor QE'] = cam1['IXON-L-888 Sensor QE'] / 100
    if cameras[1] == 'Andor_iXon':
        cam2['IXON-L-888 Sensor QE'] = cam2['IXON-L-888 Sensor QE'] / 100        

    #for each wavelength in lamdas, find the QE value for the camera 
    QE_cameras = np.empty((len(Lambdas),len(cameras))) #initialize QE for both cameraize, lambda by 2 (one col per camera)
    QE_cameras[:] = np.nan
    #loop through wavelength in Lambdas
    for i,wvl in enumerate(Lambdas):
        #chekc which cameras are requested
        if cameras[0] == 'BSI_Prime_Express':
            wvl_index = cam1.loc[cam1['BSI_New']== wvl].index.values #get index of wvl in QE curve 
        if cameras[1] == 'BSI_Prime_Express':
            wvl_index2 = cam2.loc[cam2['BSI_New']== wvl].index.values #get index of wvl in QE curve 
        if cameras[0] == 'Andor_iXon':
            wvl_index = cam1.loc[cam1['Wavelength (nm)']== wvl].index.values #get index of wvl in QE curve 
        if cameras[1] == 'Andor_iXon':
            wvl_index2 = cam2.loc[cam2['Wavelength (nm)']== wvl].index.values #get index of wvl in QE curve 
        #if wavelengt is found, save the QE values
        if len(wvl_index)> 0:
            QE_cameras[i,0] = cam1.iloc[int(wvl_index),1]#get QE for BSI. This will be in second col
        if len(wvl_index2)> 0:
            QE_cameras[i,1] = cam2.iloc[int(wvl_index2),1]#get QE for iXon. This will be in first col
    #interpolate nans
    for icam in range(2):
        y = QE_cameras[:,icam]
        nans, x= nan_helper(y)
        y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        QE_cameras[:,icam] = y
    return QE_cameras
   


#fucntion for getting filter data
def get_em_filters(filter_folder, filters, wavelengths):
    #assemble filter transmission spectra for input filters. For certain filters the data
    #is missing values, so this function will also fill those in 
    #Returns wavelength by n filters array of transmission data 
# example:
# filter_folder = './Drangonfly_transmission_spectra/Semrock_filters_bs/'
# wavelengths = np.arange(300, 801,1)#define wavelengths 
# filters = ['TR-DFLY-F521-038' , 'TR-DFLY-F698-077'] #input filter names  
# filter_trans = get_em_filters(filter_folder, filters, Lambdas)
# plt.plot(Lambdas, filter_trans[:,1]) 
    
    all_filters = filters[0] + filters[1] #convert to list
    filter_trans = np.zeros((len(wavelengths),len(all_filters))) #initialize filters 
    colnames = ['Wavelength','Transmission']
    #loop through filters to assemle filter transmission 
    for i, iname in enumerate(all_filters):
        fname = filter_folder + iname + '.txt'#assemble file name 
        cur_filt = pd.read_csv((fname),delimiter='\t', skiprows=0,names=colnames)
        for iwvl,wvl in enumerate(wavelengths):
            wvl_index = cur_filt.loc[cur_filt['Wavelength']== wvl].index.values #get index of wvl in filter data
        #if wavelengt is found, save the values
            if len(wvl_index)> 0: 
                filter_trans[iwvl,i] = cur_filt.iloc[wvl_index,1]#get values
        #521 filter has missing values, fill those in 
        if iname == 'TR-DFLY-F521-038':
            missing_wvl = np.arange(505,538,1)
        #get indeces of missing wavlengths 
            fill_in = ((wavelengths[:, None] == missing_wvl).argmax(axis=0))
            Max_val = (max(filter_trans[:,i]))
            filter_trans[fill_in,i] = Max_val #for 521, fill in max val(.97 )for missing vals 
        #698 filter has missing values, fill those in 
        elif iname == 'TR-DFLY-F698-077':
            missing_wvl = np.arange(663,734,1)
        #get indeces of missing wavlengths 
            fill_in = ((wavelengths[:, None] == missing_wvl).argmax(axis=0))
            Max_val = (max(filter_trans[:,i]))
            filter_trans[fill_in,i] = Max_val #for 698, fill in max valfor missing vals 
    return filter_trans

#beam splitter
def get_beam_spliiter(bs_folder, bs, wavelengths):
    #retrieve beam splitter data. returns wavelength by 2(1-Transmission, Transmission) array
    #no reflection data so it's calculated as 1-transmission. 
    #Transmission lets 565+ light through 
    #example:
        # bs_folder = './Drangonfly_transmission_spectra/Semrock_filters_bs/'
        # wavelengths = np.arange(300, 801,1)#define wavelengths 
        # bs = ['TR-DFLY-CMDM-565']
        # beam_split = get_beam_spliiter(bs_folder, bs, wavelengths)
        # plt.plot(wavelengths, beam_split[:,:]) 
        # plt.title('Beam Splitter??')
        # plt.legend(['1 - Transmission','Transmission',])
        # plt.xlabel('Wavelength')
        # plt.ylabel('Tramsmission')
        # plt.show()

    colnames = ['Wavelength','Transmission']
    beam_split = np.zeros((len(wavelengths),2)) #initialize beam spliter data, col 1 is blue, col 2 is red  

    fname = bs_folder + bs[0] + '.txt'#assemble file name 
    # print(fname)
    bs = pd.read_csv((fname),delimiter='\t', skiprows=0,names=colnames) #load data

    # print(cur_bs)
    for iwvl,wvl in enumerate(wavelengths):
        wvl_index = bs.loc[bs['Wavelength']== wvl].index.values #get index of wvl in  data
        if len(wvl_index)> 0: 
            beam_split[iwvl,1] = bs.iloc[wvl_index,1]#get values    
    beam_split[:,0] =  1 - beam_split[:,1]
    return(beam_split)

##function for retrieving default file paths
def get_filepaths(datafolder, **kwargs):
    #fix: doesn't work if only filters is provided as kwarg
    #assemble paths names of parts in Dragonfly_tramsission_spectra folder
    #input: location of folder where Dragonfly_transmission_spectra folder lives
    #optional inputs: non default filters (default: [['TR-DFLY-F450-050','TR-DFLY-F600-050'],['TR-DFLY-F521-038','TR-DFLY-F698-077']] 
        #filters must be defined as they are in filters folder in Dragonfly Transmission Spectra/Data/semrock_filters_bs folder
        #or non default cameras( default: )
    #assemble name of filters and beam spliter
    #returns dictionary of paths and names of default filters
    bsi_path =  datafolder + 'Drangonfly_transmission_spectra/BSI_Prime_Express/BSI_Prime_Express_QE.csv'
    ixon_path = datafolder + 'Drangonfly_transmission_spectra/iXonCamera/IXON-L-888 Sensor QE.csv'
    laser_file = datafolder + 'Drangonfly_transmission_spectra/Lasers/Laser_lines.csv'
    laser_widths = datafolder +  'Drangonfly_transmission_spectra/Lasers/Andor ILE lasers - Sheet1.csv'
    dichroic_file = datafolder + 'Drangonfly_transmission_spectra/Quad_pass_filter/Dichroic_transmission.csv'
    filter_folder = datafolder + 'Drangonfly_transmission_spectra/Semrock_filters_bs/'
    bs_folder = datafolder + 'Drangonfly_transmission_spectra/Semrock_filters_bs/'

    bs = ['TR-DFLY-CMDM-565'] #name of beam spliter
    if kwargs:
        if kwargs['filters']:
            print('using requested filters...',kwargs['filters'])
            filters = kwargs['filters']

        if kwargs['cameras']: 
            cameras = kwargs['cameras']
            print('using requested cameras...')
            default_cams = ['Andor_iXon','BSI_Prime_Express']
            cam1 = kwargs['cameras'][0]
            cam2 = kwargs['cameras'][1]
            # pdb.set_trace()
            if cam1 == default_cams[0]:
                cam1_path = ixon_path
            elif cam1 == default_cams[1]:
                cam1_path = bsi_path
            else:
                print('Specficied camera1 not found, try: Andor_iXon (or) BSI_Prime_Express')
            if cam2 == default_cams[0]:
                cam2_path = ixon_path
            elif cam2 == default_cams[1]:
                cam2_path = bsi_path
            else:
                print('Specficied camera2 not found, try: Andor_iXon (or) BSI_Prime_Express')
    else:
        print('using default fitlers...', [['TR-DFLY-F450-050','TR-DFLY-F600-050'],['TR-DFLY-F521-038','TR-DFLY-F698-077']])
        filters = [['TR-DFLY-F450-050','TR-DFLY-F600-050'],['TR-DFLY-F521-038','TR-DFLY-F698-077']] #input filter names  
        print('using default cameras (Andor iXon | BSI Prime Express ')
        cam1_path = ixon_path 
        cam2_path = bsi_path
        cameras = ['Andor_iXon','BSI_Prime_Express']
    # cam1_path =
    # cam2_path =  
    paths = {"cameras":cameras,"cam1_path":cam1_path, "cam2_path":cam2_path,"laser_file":laser_file, "laser_widths":laser_widths, "dichroic_file":dichroic_file,"filter_folder":filter_folder,"bs_folder":bs_folder,"filters":filters,"bs":bs}
    # paths = {"bsi_path":bsi_path, "ixon_path":ixon_path,"laser_file":laser_file, "laser_widths":laser_widths, "dichroic_file":dichroic_file,"filter_folder":filter_folder,"bs_folder":bs_folder,"filters":filters,"bs":bs}
    return(paths)

#function for retreiving all data together in a list
def get_spectra(FPs, paths, laser_lines,**kwargs):
    #function for retreving QE curves, filter spectra, FP spectra, dichroic mirror, laser data
    #inputs: 
    # FPs: list of flourescent proteins in experiment
    # paths: output of get_filepaths function
    #laser_lines, waveslengths of lasers used (ie 405)
        #optional input:
        #'beamsplitter = 'none'
            #this is for when no beamsplitter was used 
    #outputs Exicitation and Emmission for FPs (2 (ex/em) by n (#of FPs) by spectra (arranged from wavelength 300 to 800))
    #lasers will have a combined wavelength by lasers excitation data for the each pair inputted lasers

    EX_EM, QY, Lambdas = get_FP_spectra(FPs) #get spectra
    # QE_cameras = get_QEs(Lambdas,paths['bsi_path'],paths['ixon_path']) #get camera QE
    QE_cameras = get_QEs(Lambdas,paths['cam1_path'],paths['cam2_path'],paths['cameras']) #get camera QE
    #load saved laser lines
    all_lasers = pd.read_csv(paths['laser_file'])
    
    lasers = all_lasers.loc[:, laser_lines[0]+ laser_lines[1]] #should this be changed to not hard code adding the first two? what you we had three pairs?
    #load laser widths 
    #count the number of total lasers, sorry I didn't know how to easily do this :( 
    n_lasers = 0
    for laser_pair in laser_lines:
        for laser in laser_pair:
            n_lasers  = n_lasers + 1

    laser_widths_csv = pd.read_csv(paths['laser_widths'])
    laser_widths = np.zeros((1,n_lasers))
    ind = -1
    for laser_pair in laser_lines:
        for laser in laser_pair:
            ind = ind + 1
            laser_index = np.argwhere(laser_widths_csv['Wavelength'].values == int(laser)) #get index of laser
            laser_widths[:,ind] = laser_widths_csv['Width (+-)'].iloc[int(laser_index)]


    #load dichoric
    dichroic = pd.read_csv(paths['dichroic_file'])

    #load filters
    filter_trans = get_em_filters(paths['filter_folder'], paths['filters'], Lambdas)
    #load beam splitter 
    beam_split = get_beam_spliiter(paths['bs_folder'], paths['bs'], Lambdas)
        #if beam splitter is set to none, change values to 1 
    if kwargs:
        if kwargs['beamsplitter'] == 'none':
            beam_split[:,:] = 1
    #assemble dict
    specdata = {"Lambdas":Lambdas, "EX_EM":EX_EM,"cameras":paths['cameras'],"QE_cameras":QE_cameras,"lasers":lasers,"laser_widths":laser_widths,"dichroic":dichroic,"filters":paths['filters'],"filter_trans":filter_trans,"beam_split":beam_split,"QY":QY, "FPs":FPs}
    return(specdata)
