"""Stores all alert info + helper functions."""
import dataclasses
import os
from dataclasses import dataclass, field
from astropy.time import Time
import requests
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

@dataclass
class Alert:
    """Holder for alert-specific information."""

    event_id: str = ""
    time: Time = field(default_factory=Time)
    alert_type: str = ""
    save_path: str = "."
    event_dict: Dict[str, float] = field(default_factory=dict)
    far: float = 1.0
    significance: bool = False
    
    def __post_init__(self):
        # make all related directories
        self.directory = os.path.join(
            self.save_path,
            f"Events/{self.event_id}/{self.time.mjd}/"
        )
        os.makedirs(self.directory, exist_ok=True)
        
        self.overview_file = os.path.join(self.directory, 'Basic_Info_About_Event.txt')
        self.multiorder_file = os.path.join(
            self.directory,
            f'multiorder_fits_{self.event_id}.fits'
        )
        self.singleorder_file = os.path.join(
            self.directory,
            f'flattened_multiorder_fits_{self.event_id}.fits'
        )
        self.fits_url = f'https://gracedb.ligo.org/api/superevents/{self.event_id}/files/bayestar.multiorder.fits'
        self.gracedb_site =  f'https://gracedb.ligo.org/superevents/{self.event_id}/'
        self.max_dist = 300. # cutoff based on what could be reasonably followed up
        
        self.test_event = ( self.event_id[0] in ['M', 'T'] )
        
    
    def read_in_skymap(self):
        try:
            print("Saving multiorder file")
            #download the multi-order fits file from the fits_url file location
            req = requests.get(self.fits_url)
            file = open(self.multiorder_file, 'wb')
            for chunk in req.iter_content(100000):
                file.write(chunk)
            file.close()
            return True
        except:
            print("URL request error")
            return False
        
        
    def flatten_skymap(self):
        os.system(
            f'ligo-skymap-flatten {self.multiorder_file} {self.singleorder_file} --nside 256'
        )
        
    
    def plot_skymap(self):
        """Plot the sky-localization from the flattened, single-order fits file.
        """
        plot_fn = os.path.join(self.directory, "skymap_plot.pdf")
        os.system(
            f"ligo-skymap-plot {self.singleorder_file} -o {plot_fn}"
        )
        
        
    def plot_piechart(self):
        try:
            event_classes, event_probs = np.asarray(list(self.event_dict.items())).T
        except:
            print("Event dict not assigned to alert.")
            return
        
        event_probs = event_probs.astype(float)
        event_classes = event_classes[event_probs > 0]
        event_probs = event_probs[event_probs > 0]
        
        # Making a pie chart of the type of event for the email
        fig1, ax1 = plt.subplots()
        patches,texts = ax1.pie(event_probs, startangle=90)
        ax1.legend(patches, event_classes, loc="best")
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig(
            os.path.join(self.directory,'piechart.pdf')
        )
        
        
    def unload_skymap(self):
        try:
            skymap, header = hp.read_map(
                self.singleorder_file, h=True, verbose=False
            )
            header = dict(header)
            obs_time = Time(header['DATE-OBS'],format='isot',scale='utc')
            time = Time.now()

            self.dist_mu = float(header['DISTMEAN'])
            self.dist_std = float(header['DISTSTD'])
            return skymap
        except:
            print("Unable to read map")
            return None

    