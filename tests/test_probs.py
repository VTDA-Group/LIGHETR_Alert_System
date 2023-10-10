import pytest
import os
from astropy.time import Time
import numpy as np

from lighetr_alert_system.observatory import Observatory
from lighetr_alert_system.prob_obs import prob_observable

def test_prob_single_pixel(dummy_alert, midnight_time, het_loc):
    """Test that prob_obs is doing what's expected."""
    dummy_alert.read_in_skymap()
    dummy_alert.flatten_skymap()
    skymap = dummy_alert.unload_skymap()
    
    # put all probability in one cell
    skymap[:] = 0.
    skymap[0] = 1.
    
    full_sky_obs = Observatory(
        "GEN",
        (-90., 90.),
        het_loc,
    )
    _, frac_visible_gen = prob_observable(
        skymap, full_sky_obs, midnight_time, dummy_alert.directory
    )
    assert frac_visible_gen == 1
    assert skymap[0] == 1.
    assert np.all(skymap[1:] == 0.)
    assert os.path.exists(
        os.path.join(
            dummy_alert.directory,
            "GEN_Full_Visibility.pdf"
        )
    )
    
   
def test_prob_even_dist(dummy_alert, midnight_time, het_loc):
    """Test that prob_obs is doing what's expected."""
    dummy_alert.read_in_skymap()
    dummy_alert.flatten_skymap()
    skymap = dummy_alert.unload_skymap()
    
    # put all probability in one cell
    skymap[:] = 1. / len(skymap)
    
    full_sky_obs = Observatory(
        "GEN",
        (-90., 90.),
        het_loc,
    )
    _, frac_visible_gen = prob_observable(
        skymap, full_sky_obs, midnight_time, dummy_alert.directory
    )
    assert frac_visible_gen < 1
    assert frac_visible_gen > 0
    assert np.any(skymap == 0.) # alters non-visible areas
    assert np.any(skymap == 1. / len(skymap))
    assert os.path.exists(
        os.path.join(
            dummy_alert.directory,
            "GEN_Full_Visibility.pdf"
        )
    )
    
    
    