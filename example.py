import util
import pdb

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

paths = util.get_filepaths(datafolder)

FPs = 'mTagBFP2','mNeonGreen','mScarlet','iRFP670' #these have to be exact matches as FPbase name
laser_powers = [10, 10, 10, 10]
laser_lines = ['488','514','561','685'] #needs to be in quotes so can index into the csv
exposure_times = [100, 100]
specdata = util.get_spectra(FPs, paths, laser_lines,exposure_times)
pdb.set_trace()
print(specdata.keys())
#use specdata to assemble unmixing matrix 

#define experimental data
#get unmixing matrix
#unmix data 
