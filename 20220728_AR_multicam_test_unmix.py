from calendar import c
import numpy as np
import util
import popmat
import pdb

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

paths = util.get_filepaths(datafolder) #there's optional inputs for different filter pairs!

FPs = [['mTagBFP2','Texas Red'],['BODIPY FL goat anti-mouse IgG antibody in pH 7.2','Texas Red']] #these have to be exact matches as FPbase name
#alternatives to scarlet are mStrawberry, mTangerine, but QY is low
exc_lines = [['405','561'],['488','637']]# input paired laser lines #needs to be in quotes so can index into the csv
laser_powers = [1,1,1,1] #Percentages 
exposure_times = [100, 200]
specdata = util.get_spectra(FPs, paths, exc_lines)
# assemble unmixing matrix
c_2d = popmat.populate_matrix(specdata, exc_lines, laser_powers, exposure_times)



