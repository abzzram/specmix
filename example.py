import numpy as np
import util
import popmat
import pdb

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

paths = util.get_filepaths(datafolder) #add optional inputs to change filters 

FPs = ['mTagBFP2','mNeonGreen','mScarlet','iRFP670'] #these have to be exact matches as FPbase name
exc_lines = [['405','561'],['488','637']]# input paired laser lines #needs to be in quotes so can index into the csv
# exc_lines = [['405','561'],['488','561']]# input paired laser lines #needs to be in quotes so can index into the csv

##TO DO, fix get_spectra so that new exclines works, then fix popmat to add the paired laser lines 
em_filters = ['>530', '<530'] # TODO
laser_powers = [.5, .6, .7, .8] #Percentages 
exposure_times = [100, 100]
specdata = util.get_spectra(FPs, paths, exc_lines)



# assemble unmixing matrix
c_2d = popmat.populate_matrix(specdata, exc_lines, em_filters, laser_powers, exposure_times)
pdb.set_trace()
# simulate
abundances = np.array([1, 0.5, 0.2, 0.7])
intensities = c_2d @ abundances

# inverse
abundance_est = np.linalg.pinv(c_2d) @ intensities

