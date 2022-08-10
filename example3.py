#testing script for loading example zstack images and unmixing
import numpy as np
import util
import popmat
import pdb

#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 
filts = [['TR-DFLY-F600-050'],['TR-DFLY-F698-077']]
cams = ['BSI_Prime_Express','BSI_Prime_Express'] #images taken on same cameras 
FPs = [['mScarlet'],['iRFP670']] #these have to be exact matches as FPbase name

exc_lines = [['561'],['637']]# input paired laser lines #needs to be in quotes so can index into the csv

laser_powers = [ 1,1] #Percentages 
exposure_times = [100, 100]
c_2d = util.specmix_matrix(datafolder,FPs, exc_lines, laser_powers, exposure_times,cameras = cams, filters= filts,beamsplitter = 'false')

