#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[18]:


#define functions 
# get_ipython().run_line_magic('run', "-i 'Functions for spectral modeling.py' #define functions")
import util

# In[82]:

print("Script starting...")
#load FP data and retrieve wavelengths 
FPs = 'mTagBFP2','mNeonGreen','mScarlet','iRFP670' #these have to be exact matches as FPbase name
datafolder = './Data/'
folder_bsi = datafolder + 'Drangonfly_transmission_spectra/BSI_Prime_Express/'
folder_ixon = datafolder + 'Drangonfly_transmission_spectra/iXonCamera/'
EX_EM, QY, Lambdas = util.get_FP_spectra(FPs) #get spectra
QE_cameras = util.get_QEs(Lambdas,folder_bsi,folder_ixon) #get camera QE

#load saved laser lines
laser_file = datafolder + './Drangonfly_transmission_spectra/Lasers/Laser_lines.csv'
lasers = pd.read_csv(laser_file)

#load dichoric
dichroic_file = datafolder + './Drangonfly_transmission_spectra/Quad_pass_filter/Dichroic_transmission.csv'
dichroic = pd.read_csv(dichroic_file)
#load filters
filter_folder = datafolder + './Drangonfly_transmission_spectra/Semrock_filters_bs/'
filters = ['TR-DFLY-F521-038' , 'TR-DFLY-F698-077'] #input filter names  
filter_trans = util.get_em_filters(filter_folder, filters, Lambdas)
#load beam splitter 
bs_folder = datafolder + './Drangonfly_transmission_spectra/Semrock_filters_bs/'
bs = ['TR-DFLY-CMDM-565']
beam_split = util.get_beam_spliiter(bs_folder, bs, Lambdas)

##one function that has a list of info for functions, generates cross talk matrix 
#one unmix function (returns csv that has rows/cols labeled)
#two line workflow 


# In[83]:

## save the figs,
print("Starting plots...")
FP_to_plot = 3
plt.plot(Lambdas,EX_EM[0,FP_to_plot,:])#plot excitation for 3rd protein in our FPs
plt.plot(Lambdas,EX_EM[1,FP_to_plot,:]) #plot emission 
plt.title(FPs[FP_to_plot])
plt.legend(('Ex','Em'))
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')
plt.show()


# In[84]:


#plot QE
plt.plot(Lambdas, QE_cameras[:,0])
plt.plot(Lambdas, QE_cameras[:,1])
plt.title('Camera Quantum Efficieny Curves')
plt.legend(('BSI Prime Express','Andor iXon'))
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')
plt.show


# In[85]:


#plot filters
plt.plot(Lambdas, filter_trans[:,:]) 
plt.title('Emission Fitlers')
plt.legend(filters)
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')
plt.show()


# In[86]:


#plot beam splitter
plt.plot(Lambdas, beam_split[:,:]) 
plt.title('Beam Splitter')
plt.legend(['1 - Transmission','Transmission',])
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')
plt.show()


# In[87]:


#plot laser lines
lasers = pd.read_csv(laser_file)
plt.plot(Lambdas,lasers.iloc[:,1:])
plt.title('Laser Lines')
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')


# In[88]:


#plot dichroic 
#load dichroic 
dichroic_file = datafolder + './Drangonfly_transmission_spectra/Quad_pass_filter/Dichroic_transmission.csv'
dichroic = pd.read_csv(dichroic_file)
#plot dichroic
plt.plot(Lambdas,dichroic.iloc[:,1:])
plt.title('Dichroic Mirror')
plt.xlabel('Wavelength')
plt.ylabel('Tramsmission')
plt.show()


# In[ ]:




