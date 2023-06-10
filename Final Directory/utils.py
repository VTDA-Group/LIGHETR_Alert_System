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

    
def glade_completeness(dists):
    """
    Return GLADE's completeness at array of specified distances.
    Assumes more conservative completeness obtained by comparison with
    expected constant luminosities in each cell of volume space.
    
    Parameters
    ----------
    dist : numpy array, float
        galaxy distances, in Mpc
    """
    # in intervals of 50 Mpc, then 100 Mpc
    dist_refs = np.array([0, 50., 100., 150., 200., 300., 400., 500., 600., 700., 800.])
    dist_to_comp = np.array([100., 80, 65., 50., 45., 42., 40., 35., 30., 25., 20.])
    
    closest_dist_idxs = np.array([np.argmin(d - dist_refs) for d in dists)])
    return dist_to_comp[closest_dist_idxs]
    
    

def calculate_integrated_completeness(dist_mu, dist_sigma):
    """
    Assuming we're using a 3-sigma distance cut for each pixel, determine
    the completeness we expect for that pixel integrated across that distance range.
    
    Approximates integral as 20-bin Riemann sum.
    """
    dist_min = dist_mu - 3. * dist_sigma
    dist_max = dist_mu + 3. * dist_sigma
    
    dist_bins = np.array([np.linspace(dist_min[i], dist_max[i], num=20) for i in range(len(dist_mu))])
    
    pix_comp_per_dist = np.array([glade_completeness(d_bin) for d_bin in dist_bins])
    pix_comp_summed = np.sum(3. * dist_bins**2 * pix_comp_per_dist, axis=1)
    pix_comp_avged = pix_comp_summed * (dist_bins[:,1] - dist_bins[:,0]) / (dist_max**3 - dist_min**3)
    return pix_comped_avged

    
    
def adjust_probs_for_completeness(normed_probs, dists):
    """
    Reduce effective probabilities to account for GLADE+ completeness at different
    distances. Assumes that, for each galaxy at distance d, there are (1 / completeness - 1)
    similar galaxies not visible at that distance, so the adjusted probability is:
    
    prob_adj = prob_norm / (1 / completeness) = prob_norm * completeness
    
    Note that
    
    Parameters
    ----------
    normed_probs : numpy array
        normed probabilities of each galaxy, accounting for GW signal pixel probs
    dists : numpy array
        galaxy distances, in Mpc
    """
    
