import numpy as np

def populate_matrix(specdata, exc_lines, em_filters, FPs):
    # Input: a dictionary of spectral data
    # Input: 
    # Output: a numpy array with our spectral mixing model,

    # Generate empty matrix c_m,m',n
    # For example:
    # 2 exc lines
    # 2 em filters
    # 4 FPS
    # Empty 2 x 2 x 4 array
    c_3d = np.zeros((len(exc_lines), len(em_filters), len(FPs)))

    # Loop through the entries and populate
    for m, exc_line in enumerate(exc_lines):
        for mp, em_filter in enumerate(em_filters):
            for n, FP in enumerate(FPs):
                # Collect all terms
                laser_spec =  np.array(specdata['lasers'])[:,m]
                dichroic_spec = np.array(specdata['dichroic'])[:,1]
                abs_spec = np.array(specdata['EX_EM'][0,n])
                exc_part = np.sum(laser_spec * dichroic_spec * abs_spec) # Collect excitation bits and sum/integrate
                em_part = (mp+1)*(n+1) # Collect emission bits and sum/integrate

                c_3d[m,mp,n] = exc_part*em_part

    # Collapse over first two dimensions
    c_2d = np.reshape(c_3d, (c_3d.shape[0]*c_3d.shape[1], c_3d.shape[2]))
    
    return c_2d
