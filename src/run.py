import matplotlib.pyplot as plt
from astropy.time import Time
import json
import os
import pdb
from hop import Stream
from hop.io import StartPosition
import healpy as hp
from astropy.io import fits
import requests
import prob_obs_gen
import prob_obs_HET
import get_galaxies
import get_LST

import time as TIme
import numpy as np

from contact import *
from alert import *


def process_fits(fits_file, alert_message = None, skip_test_alerts = True):
        '''
        The format of these alerts is given in this website:
        https://emfollow.docs.ligo.org/userguide/content.html
        
        So, the probabilities of BBH, BNS, and so on is located in:
        alert_message['event']['classification']['BNS']
        alert_message['event']['classification']['BBH']
        alert_message['event']['classification']['NSBH']
        alert_message['event']['classification']['Noise']
        
        it's unclear whether the key for noise will be either 'Noise' or 'Terrestrial', so I'll make the code flexible enough to use both
        
        If you look at alert_message['superevent_id'], it will be [{T,M}]SYYMMDDabc. The first character can be either T or M. T means test, M means mock. Then, the S just means it's a superevent. Then you have the name of the event.
        
        '''
        contact_list_all = get_contact_lists(
            file_loc = 'contact_all_LIGO.json',
            people_to_contact = people_to_contact
        )
        contact_list_HET = get_contact_lists(
            file_loc = 'contact_all_LIGO.json',
            people_to_contact = people_to_contact
        )
        
        alert_type = alert_message['alert_type']
        
        superevent_id = alert_message['superevent_id']
        gracedb_site = 'https://gracedb.ligo.org/superevents/'+str(superevent_id)+'/'
        test_event = ( superevent_id[0] in ['M', 'T'] )
            
        if skip_test_alerts and test_event: #if we want to ignore test events, and this is a test event, ignore it.
            return
        
        print("\n\n=============================\nFound a real LIGO event: "+str(superevent_id)+" | "+str(alert_type)+" Alert\n")
        
        # create alert
        alert = Alert(
            alert_type = alert_type,
            time = Time(alert_message['time_created']),
            event_id = superevent_id,
        )

        with open(self.overview_file, 'w+') as data_out: 
            data_out.write(
                f'Found an event: {alert.event_id} | {alert.alert_type} Alert\n' +
                f'Found at time: {alert.time.mjd} MJD\n'
            )
        
        try:
            print("Saving multiorder file")
            #download the multi-order fits file from the fits_url file location
            req = requests.get(alert.fits_url)
            file = open(alert.multiorder_file, 'wb')
            for chunk in req.iter_content(100000):
                file.write(chunk)
            file.close()
            
        except:
            print("URL request error")
            return
        
        
        #flatten the multi-order fits file into single-order
        print("Flattening...")
        alert.flatten_skymap()
        
        #get the skymap and header of the flattened, single-order fits file
        print("Processing FITS...")
        try:
            alert.unload_skymap()
            
        except:
            print("Unable to read map")
            save_burst_info(alert, contact_list_all)
            return

        alert.plot_skymap()
        
        # get classification probs
        alert.event_dict = alert_message['event']['classification']
        alert.plot_piechart()

        alert.far = alert_message['event']['far']
        alert.significance = alert_message['event']['significant']

        
        if float(far) > 3.9E-7 and significance == 'False':
            #sending emails out to only people on the contact_list_file_loc_all_events file about the alert. Because there is likely no remnant.
            
            send_low_prob_info(alert, contact_list_all, reason="significance")
            return
        
        #if it's not at least 30% a (BNS or NSBH) signal, ignore it.
        if (alert.event_dict['NSBH']+alert.event_dict['BNS'])/np.sum(alert.event_dict.values()) < 0.3:
            if sizes[3] < 0.9: #send the email if terrestrial signal < 0.9
                send_low_prob_info(alert, contact_list_all)
            else:
                print("Probably terrestrial")
            return
        
        if alert.dist_mu - alert.dist_std > alert.max_dist:
            send_low_prob_info(alert, contact_list_all, reason="distance")     
            return
        
        print("Calculating probabilities")
        #find pixels in the 90% confidence region that will be visible to HET in the next 24 hours
        skymap1 = np.copy(skymap)
        m_gen, frac_visible_gen = prob_obs_gen.prob_observable(skymap, header, time, savedir = obs_time_dir, plot=True)
        timetill90_HET, m_HET, frac_visible_HET = prob_obs_HET.prob_observable(skymap1, header, time, savedir = obs_time_dir, plot=True)
        
        print("Two-dimensionally, percentage of pixels visible to HET: "+str(frac_visible_HET*100)+"%")

        alert_message['skymap_fits'] = singleorder_file_name
        alert_message['skymap_array_HET'] = m_HET
        alert_message['skymap_array'] = m_gen
        
        cattop, logptop = get_galaxies.write_catalog(alert_message, savedir = obs_time_dir, HET_specific_constraints = False)
        
        if timetill90_HET ==-99 or frac_visible_HET <= 0.1:
            print("HET can't observe the source.")
            send_not_visible_info(alert, contact_list_all)
            return
        else:
            cattop, logptop = get_galaxies.write_catalog(
                alert_message, savedir = obs_time_dir, HET_specific_constraints = True
            )
            if len(cattop) == 0:
                return
                
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90_HET)+"\nPercentage of visible pixels to HET: "+str(round(frac_visible_HET*100, 3))+"%")
            
            get_LST.get_LST(
                savedir = alert.directory,
                targf = os.path.join(alert.directory,'HET_Visible_Galaxies_prob_list.dat')
            )
            
            send_mapped_alert_info(alert, contact_list_HET)
            

            
if __name__ == "__main__":
    ###########Things start here####################
    contact_list_file_loc = 'contact_only_HET_BNS.json'
    contact_list_file_loc_all_events = 'contact_all_LIGO.json'
    #people_to_contact = ["Karthik"]
    people_to_contact = []

    #stream_start_pos = 1600
    stream_start_pos = StartPosition.EARLIEST
    #print("Starting stream at "+str(stream_start_pos))
    #stream = Stream(start_at=stream_start_pos)

    stream = Stream()

    num_messages = 0

    print("Listening for Alerts from kafka")

    with stream.open("kafka://kafka.scimma.org/igwn.gwalert", "r") as s:
        for message in s:
            event = message.content[0]['event']
            message_content = message.content[0]
            alert_type = message_content['alert_type']
            alert_time = message_content['time_created']

            alert_time = Time(alert_time)
            print(alert_time)
            num_messages+=1
            #if num_messages > 2:
            #    sys.exit()

            skymap = None
            if 'skymap' in message_content.keys():
                skymap = message_content['skymap']
            if event is not None and 'skymap' in event.keys():
                skymap = event['skymap']

            '''send that fits file into process_fits'''
            if skymap is not None and event is not None:
                print('Calling process_fits')
                process_fits(fits_file = skymap, alert_message = message_content,  skip_test_alerts = True)

