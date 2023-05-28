import numpy as np
from ligo.skymap.distance import conditional_pdf

def get_relative_probs(M, dist, distinfo):
    """
    Weight galaxies in single pixel by mass and distance CDF.
    
    Parameters
    ----------
    M : numpy array
        Masses of the galaxies
    dist : np array
        Distances to each galaxy (in Mpc)
    distinfo :
        contains (distmu, distsigma, distnorm) for the signal at each pixel.
        
    Returns
    ----------
    1d numpy array
        Relative probs of every galaxy (not normalized or corrected for pixel prob)
    """
    #all arrays should be the same length
    assert len(M) == len(dist)
    assert len(dist) == len(distinfo)
    
    distmu, distsigma, distnorm = distinfo.T
    ignore_idxs = (distsigma <= 0) | (distnorm <= 0)
    c_pdf = conditional_pdf(dist, *distinfo.T)
    c_pdf[ignore_idxs] = 0.

    return c_pdf * M
    

def distribute_pixel_prob(pixel_probs, pixel_mappings, cat, distinfo):
    """
    Distribute the probability of a single pixel across the galaxies
    in that pixel.
    
    Parameters
    ----------
    pixel_prob : 1d np array
        Probabilities associated with each pixel (NOT PROB DENSITY)
    gal_idxs : 2d np array
        Indices of galaxy catalog associated with each pixel. One row per pixel.
        Expressed as an "adjacency matrix", 1's if gal included in pixel, 0 otherwise
    cat : pandas df
        imported GLADE catalog
    distinfo : 2d numpy array.
        contains (distmu, distsigma, distnorm) for each pixel. One 3-tuple per row.
        
    Returns
    ----------
    1d numpy array
        The probabilities of each galaxy corrected for pixel probs
    """
    
    dist = np.array(cat['dist_Mpc']) # get distances in Mpc
    stellar_M = np.array(cat['M*']) # masses
    
    relative_probs = get_relative_probs(stellar_M, dist, distinfo)
    #probs_tiled = np.repeat(relative_probs[np.newaxis,:], len(pixel_probs), axis=0)
    uniq_idxs = np.unique(pixel_mappings)
    normed_probs = np.zeros(len(relative_probs))
    
    for i in range(len(uniq_idxs)):
        normed_c = np.sum(relative_probs[pixel_mappings == i])
        normed_probs[pixel_mappings == i] = relative_probs[pixel_mappings == i] * pixel_probs[i] / normed_c
        
    """
    rel_tiled = relative_probs[np.newaxis,:] * pixel_idxs
    normed_consts = np.sum(rel_tiled, axis=1)
    
    normed_consts[normed_consts <= 0] = np.inf

    normed_probs_tiled = rel_tiled * (pixel_probs / normed_consts)[:,np.newaxis] # distributes probability 1 among the galaxies within each pixel
    
    print(pixel_probs[:2], normed_probs_tiled[:2])
    # condense back down into 1d
    normed_probs = np.sum(normed_probs_tiled, axis=0)
    """
    return normed_probs # so all probs sum to pixel_prob
    
    
def get_gal_pixel_mapping(ipix):
    """
    From list of thetas and phis, gets galaxy indices associated with each pixel.
    
    Returned as a 2d binary matrix, where each row is a pixel, and each column is a galaxy.
    Each cell is 1 if that galaxy is in that pixel, and 0 otherwise.
    
    Parameters
    ----------
    theta : 1d numpy array (TODO: update this)
        theta values in radians for galaxies
    phi : 1d numpy array
        phi values in radians for galaxies
        
    Returns
    ----------
    binary_matrix : 2d numpy array
        1 if galaxy j is in pixel i, 0 otherwise
    """
    unique_i, inverse_idxs = np.unique(ipix, return_inverse=True)
    """
    gal_idxs = np.arange(len(ipix))
    
    binary_matrix = np.zeros( (len(unique_i), len(ipix)) )
    binary_matrix[inverse_idxs, gal_idxs] = 1
    
    assert ~np.any(np.sum(binary_matrix, axis=0) != 1)
    """
    print(inverse_idxs)
    
    
    return inverse_idxs, unique_i

    
    
    
    
    
