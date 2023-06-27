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
    in that pixel. Also accounts for < 100% completeness.
    
    STEPS TO GET PROPER MISSING PROBS:
    - For each pixel:
        - Get completeness for each distance bin
        - Calculate conditional_pdf for each distance bin
        - Treat unseen galaxies per bin as an effective single galaxy
        - For that effective galaxy, ranked prob = completeness * M_bin * conditional_pdf
        - Distribute pixel prob among "effective" galaxies and real galaxies
        - Sum probs of effective galaxies as "missing prob"
    
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
    missing_prob : float
        Total probability across sky occupied by galaxies we cannot see
    """
    
    dist = np.array(cat['dist_Mpc']) # get distances in Mpc
    stellar_M = np.array(cat['M*']) # masses
    
    #uniq_idxs = np.unique(pixel_mappings)
    relative_probs = get_relative_probs(stellar_M, dist, distinfo)
    normed_probs = np.zeros(len(relative_probs))

    completeness_per_bin, dist_bins = calculate_completeness_per_bin(distinfo[:,0], distinfo[:,1])
    total_visible_mass = np.sum(stellar_M)
    vis_mass_per_pixel = total_visible_mass / len(pixel_probs)
    missing_prob = 0.
    for i in range(len(pixel_probs)):
        gals_in_pixel = np.where(pixel_mappings == i)[0]
        #if len(gals_in_pixel) == 0:
        #    #missing_prob += (1. - completeness_per_bin[i]) * pixel_probs[i]
        #    continue
        missing_M_per_bin = (1. / completeness_per_bin[i] - 1.) * distribute_val_among_volume_space(vis_mass_per_pixel, dist_bins[i])
        distinfo_single = distinfo[np.newaxis,i] * np.ones((len(dist_bins[i]),3))
        relative_miss_probs_sum = np.sum(get_relative_probs(missing_M_per_bin, dist_bins[i], distinfo_single))
        # incorporate missing galaxies when normalizing probabilities
        normed_c = np.sum(relative_probs[gals_in_pixel]) + relative_miss_probs_sum
        if normed_c == 0.:
            normed_c = 1.
        normed_probs[gals_in_pixel] = relative_probs[gals_in_pixel] * pixel_probs[i] / normed_c
        missing_prob += relative_miss_probs_sum * pixel_probs[i] / normed_c
        
    assert np.all(normed_probs >= 0)
    return normed_probs, missing_prob # so all probs sum to pixel_prob

    
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
    dist_to_comp = np.array([80., 80., 65., 50., 45., 42., 40., 35., 30., 25., 20.]) / 100.
    
    closest_dist_idxs = np.array([np.argmin((d - dist_refs)**2) for d in dists])
    return dist_to_comp[closest_dist_idxs]
    

def calculate_completeness_per_bin(dist_mu, dist_sigma):
    """
    Assuming we're using a 3-sigma distance cut for each pixel, determine
    the completeness we expect for that pixel integrated across that distance range.
    
    Approximates integral as 20-bin Riemann sum.
    """
    dist_min = dist_mu - 3. * dist_sigma
    dist_min[(dist_min <= 0.) | np.isinf(dist_min) | np.isnan(dist_min)] = 0.0001
    dist_max = dist_mu + 3. * dist_sigma
    dist_max[(dist_max <= 0.) | np.isinf(dist_max) | np.isnan(dist_max)] = 1e10
            
    dist_bins = np.array([np.linspace(dist_min[i], dist_max[i], num=20) for i in range(len(dist_mu))])
    
    pix_comp_per_dist = np.array([glade_completeness(d_bin) for d_bin in dist_bins])
    
    return pix_comp_per_dist, dist_bins
    #pix_comp_summed = np.sum(3. * dist_bins**2 * pix_comp_per_dist, axis=1)
    #pix_comp_avged = pix_comp_summed * (dist_bins[:,1] - dist_bins[:,0]) / (dist_max**3 - dist_min**3)
    #return pix_comped_avged
    

def distribute_val_among_volume_space(val, dists):
    """
    Used to distribute value accordng to weights and spherical shell size at each distance.
    """
    return dists**2 * val / np.sum(dists**2) # more volume space occupied by further out distances

    
