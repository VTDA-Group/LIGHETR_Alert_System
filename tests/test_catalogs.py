import pytest
import os
from astropy.time import Time
import numpy as np

from lighetr_alert_system.observatory import Observatory
from lighetr_alert_system.prob_obs import prob_observable
from lighetr_alert_system.get_galaxies import write_catalogs

def test_write_catalog_works(dummy_alert, midnight_time, het_loc):
    """Just checks that it works."""
    dummy_alert.read_in_skymap()
    dummy_alert.flatten_skymap()
    skymap = dummy_alert.unload_skymap()
    
    full_sky_obs = Observatory(
        "GEN",
        (-90., 90.),
        het_loc,
        alert=dummy_alert,
    )
    prob_observable(
        skymap, full_sky_obs, midnight_time, dummy_alert.directory
    )
    full_sky_obs.skymap = skymap

    cattop_all, logptop_all = write_catalogs(
        observatories = [full_sky_obs,],
        savedir = dummy_alert.directory,
    )
    assert len(cattop_all) == 1
    assert len(logptop_all) == 1
    
    assert len(cattop_all[0]) > 0
    assert len(logptop_all[0]) > 0
    
    assert np.sum(np.exp(logptop_all[0])) < 1.
    
    assert os.path.exists(
        os.path.join(
            dummy_alert.directory,
            "GEN_Visible_Galaxies_prob_list.dat"
        )
    )