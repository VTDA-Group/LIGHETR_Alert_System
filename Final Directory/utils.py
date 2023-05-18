import numpy as np

def get_relative_probs(M, dist, distmu, distsigma, distnorm):
    """
    Weight galaxies in single pixel by mass and distance CDF.
    """
    pass
    

def distribute_pixel_prob(pixel_prob, pixel_idxs, cat, distinfo):
    """
    Distribute the probability of a single pixel across the galaxies
    in that pixel.
    
    Parameters
    ----------
    pixel_prob : int
        Probability associated with that pixel
    pixel_idxs : np array
        Indices of distinfo list associated with that pixel
    cat : pandas df
        imported GLADE catalog
    distinfo : 3-element list or np array
        contains (distmu, distsigma, distnorm) for the signal
    """
    pass
    
