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
    
    c_pdf = conditional_pdf(dist, *distinfo)
    
    return c_pdf * M
    

def distribute_pixel_prob(pixel_probs, pixel_idxs, cat, distinfo):
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
    probs_tiled = np.tile(relative_probs, (1, len(pixel_probs)))
    
    normed_probs_tiled = probs_tiled * gal_idxs / np.sum(probs_tiled * gal_idxs, axis=1) # distributes probability 1 among the galaxies within each pixel
    
    # condense back down into 1d
    normed_probs = np.sum(normed_probs_tiled, axis=0)
    
    return pixel_probs * normed_probs # so all probs sum to pixel_prob
    
    
def get_gal_idxs_flattened(theta, phi):
    """
    From list of thetas and phis, gets galaxy indices associated with each pixel.
    
    Parameters
    ----------
    theta : 1d numpy array
        theta values in radians for galaxies
    phi : 1d numpy array
        phi values in radians for galaxies
    """
    pass
    
    
    
