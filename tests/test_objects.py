import os
import astropy.units as u
import astropy
import numpy as np


from lighetr_alert_system.observatory import Observatory


def test_alert(dummy_alert):
    """Test alert generation."""
    alert = dummy_alert
    assert alert.max_dist == 300.
    
    assert alert.read_in_skymap()
    
    alert.flatten_skymap()
    assert os.path.exists(
        alert.singleorder_file
    )
    
    alert.plot_skymap()
    assert os.path.exists(
        os.path.join(
            alert.directory,
            "skymap_plot.pdf"
        )
    )
    
    alert.plot_piechart()
    assert os.path.exists(
        os.path.join(
            alert.directory,
            "piechart.pdf"
        )
    )
    skymap = alert.unload_skymap()
    assert skymap is not None
    
        
def test_observatory(dummy_alert):
    """Test observatory functions."""
    loc = (-104.01472,30.6814,2025)

    observatory = Observatory(
        name='test',
        dec_range=(-60, 60), 
        loc=loc,
        alert=dummy_alert,
    )
    
    assert np.isclose(observatory.max_dec_rad, 5.*np.pi/6.)
    assert np.isclose(observatory.min_dec_rad, np.pi/6.)
    