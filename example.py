##import stuff
import json
import pandas as pd
import csv
import numpy as np
import ast
import wget
import json
import array 
from numpy import nan
# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import glob
import os.path
import util
import pdb


#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 
# #the following part should be a function because it won't change much for e
# bsi_path =  datafolder + 'Drangonfly_transmission_spectra/BSI_Prime_Express/BSI_Prime_Express_QE.csv'
# ixon_path = datafolder + 'Drangonfly_transmission_spectra/iXonCamera/IXON-L-888 Sensor QE.csv'
# laser_file = datafolder + 'Drangonfly_transmission_spectra/Lasers/Laser_lines.csv'
# dichroic_file = datafolder + 'Drangonfly_transmission_spectra/Quad_pass_filter/Dichroic_transmission.csv'
# filter_folder = datafolder + 'Drangonfly_transmission_spectra/Semrock_filters_bs/'
# bs_folder = datafolder + 'Drangonfly_transmission_spectra/Semrock_filters_bs/'
# filters = ['TR-DFLY-F521-038' , 'TR-DFLY-F698-077'] #input filter names  
# bs = ['TR-DFLY-CMDM-565'] #name of beam spliter 

# paths = {"bsi_path":bsi_path, "ixon_path":ixon_path,"laser_file":laser_file, "dichroic_file":dichroic_file,"filter_folder":filter_folder,"bs_folder":bs_folder,"filters":filters,"bs":bs}#Define Experiment specfic info 
paths = util.get_filepaths(datafolder)

FPs = 'mTagBFP2','mNeonGreen','mScarlet','iRFP670' #these have to be exact matches as FPbase name
laser_powers = [10, 10, 10, 10]
laser_lines = ['405', '637']
exposure_times = [100, 100]
specdata = util.get_spectra(FPs, paths, laser_lines,exposure_times)
pdb.set_trace()
#use specdata to assemble unmixing matrix 

#define experimental data
#get unmixing matrix
#unmix data 