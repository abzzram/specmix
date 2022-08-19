from calendar import c
import numpy as np
import util
import popmat
import pdb
import tifffile
import matplotlib.pyplot as plt
#define data paths for QE curves, should data live in the server? 
datafolder =  './Data/' #where the Dragonfly_transmission_spectra folder live 

paths = util.get_filepaths(datafolder) #there's optional inputs for different filter pairs!

FPs = [['DAPI_chroma','Texas Red'],['BODIPY FL goat anti-mouse IgG antibody in pH 7.2','mTurquoise2']] #these have to be exact matches as FPbase name
#alternatives to scarlet are mStrawberry, mTangerine, but QY is low
exc_lines = [['405','561'],['488','637']]# input paired laser lines #needs to be in quotes so can index into the csv
laser_powers = [1,1,1,1] #Percentages 
exposure_times = [100, 200]
cams = ['Andor_iXon','BSI_Prime_Express']
filts = [['TR-DFLY-F450-050','TR-DFLY-F600-050'],['TR-DFLY-F525-050','TR-DFLY-F698-077']]
c_2d = util.specmix_matrix(datafolder,FPs, exc_lines, laser_powers, exposure_times,cameras = cams,filters = filts)

# specdata = util.get_spectra(FPs, paths, exc_lines)
# assemble unmixing matrix
# c_2d = popmat.populate_matrix(specdata, exc_lines, laser_powers, exposure_times)
#load in data 
filepath = './Data/Example_datasets/20220728_AR_multicam_test/Invitrogen_slide_3/Invitrogen_slide_3_MMStack_Pos0.ome.tif'
with tifffile.TiffFile(filepath) as tiff:
    image_data = np.array([page.asarray() for page in tiff.pages if page.asarray() is not None])
# imgplot = plt.imshow(image_data[2])
num_stacks = 21
num_chans =4

zstacks = np.zeros((num_stacks,num_chans,1024,1024))
#first et chans 1 and 2 from first 42 slices of the uploaded image 
istack = 0
for islice in range(0,42): #get first two chanels 
    #get chan 1 from even slices below 42
    if islice % 2 == 0:
        zstacks[istack,0,:,:] = image_data[islice]
    #get chan 2 from odd slices below 42
    elif islice % 2 != 0:
        zstacks[istack,1,:,:] = image_data[islice]
        print('done getting odds for isclice ', islice)
        print('istack = ', istack)
    if islice % 2 != 0: #go to next stack on next odd number
        istack = istack +1
#now get chans 3 and 4 from the last 42 slices of the uploaded image 
istack = 0 #reset stack number so it startsback from 0
for islice in range(42,84): #get first two chanels 
    #get chan 1 from even slices below 42
    if islice % 2 == 0:
        zstacks[istack,2,:,:] = image_data[islice]
        print('done getting evens for islice ', islice)
        print('istack = ', istack)
    #get chan 2 from odd slices below 42
    elif islice % 2 != 0:
        zstacks[istack,3,:,:] = image_data[islice]
        print('done getting odds for isclice ', islice)
        print('istack = ', istack)
    if islice % 2 != 0: #go to next stack on next odd number
        istack = istack +1

         
pdb.set_trace()


#frist ag = all dimensions, ij are frist two dims
# unmixed = np.einsum('ij,xyj->xyi', np.linalg.pinv(c_2d), zstacks)
#
unmixed = np.einsum('ij,zjxy->zixy', np.linalg.pinv(c_2d), zstacks)
#might need to switch first ij
#save output as tiff 
plt.imshow(unmixed[10,0,...])
plt.show()


