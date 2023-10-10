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
import numpy as np

from lighetr_alert_system.contact import *
from lighetr_alert_system.alert import *
from lighetr_alert_system.observatory import *
from lighetr_alert_system.prob_obs import prob_observable
from lighetr_alert_system import get_galaxies


def process_fits(alert_message, save_path='.', people_to_contact = None, skip_test_alerts = True):
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
        dir_path = os.path.dirname(os.path.realpath(__file__))
        contact_lists_followup = get_contact_lists(
            file_loc = os.path.join(
                dir_path,
               "../data/contact_only_events.json"
            ),
            people_to_contact = people_to_contact
        )
        contact_lists_all = get_contact_lists(
            file_loc = os.path.join(
                dir_path,
                "../data/contact_all_alerts.json"
            ),
            people_to_contact = people_to_contact
        )
        
        alert_type = alert_message['alert_type']
        
        superevent_id = alert_message['superevent_id']
        gracedb_site = 'https://gracedb.ligo.org/superevents/'+str(superevent_id)+'/'
        test_event = ( superevent_id[0] in ['M', 'T'] )
            
        if skip_test_alerts and test_event: #if we want to ignore test events, and this is a test event, ignore it.
            print("TEST EVENT - SKIPPING")
            return False
        
        print("\n\n=============================\nFound a real LIGO event: "+str(superevent_id)+" | "+str(alert_type)+" Alert\n")
        
        # create alert
        alert = Alert(
            alert_type = alert_type,
            time = Time(alert_message['time_created']),
            event_id = superevent_id,
            event_dict = alert_message['event']['classification'],
            far = alert_message['event']['far'],
            significance = alert_message['event']['significant'],
            save_path = save_path
        )
        
        with open(alert.overview_file, 'w+') as data_out: 
            data_out.write(
                f'Found an event: {alert.event_id} | {alert.alert_type} Alert\n' +
                f'Found at time: {alert.time.mjd} MJD\n'
            )
        
        if not alert.read_in_skymap():
            return False
        
        #flatten the multi-order fits file into single-order
        print("Flattening...")
        alert.flatten_skymap()
        
        #get the skymap and header of the flattened, single-order fits file
        print("Processing FITS...")
        skymap = alert.unload_skymap()
        if skymap is None:  
            save_burst_info(alert, contact_lists_all['email'])
            return False
        
        alert.plot_skymap()
        alert.plot_piechart()

        
        if float(alert.far) > 3.9E-7 and alert.significance == 'False':
            #sending emails out to only people on the contact_list_file_loc_all_events file about the alert. Because there is likely no remnant.
            send_low_prob_info(alert, contact_lists_all['email'], reason="significance")
            return False
        
        #if it's not at least 30% a (BNS or NSBH) signal, ignore it.
        if (alert.event_dict['NSBH']+alert.event_dict['BNS'])/np.sum(list(alert.event_dict.values())) < 0.3:
            noise_key = 'terrestrial'
            if noise_key not in alert.event_dict:
                noise_key = 'noise'
            if alert.event_dict[noise_key] < 0.9: #send the email if terrestrial signal < 0.9
                send_low_prob_info(alert, contact_lists_all['email'])
            else:
                print("Probably terrestrial")
            return False
        
        if alert.dist_mu - alert.dist_std > alert.max_dist:
            send_low_prob_info(alert, contact_lists_all['email'], reason="distance")     
            return False
        
        current_time = Time.now()
        print("Calculating probabilities")
        
        # also send to FLWO (same observatory)
        MMT_observatory = Observatory(
            name='MMT',
            dec_range=(-30.0, 90.0),
            loc=(-110.8852,31.6889,2616),
            alert=alert
        )
        
        magellan_observatory = Observatory(
            name='Magellan',
            dec_range=(-90.0, 30.0),
            loc=(-70.683575,-29.048217,2516),
            alert=alert,
        )

        #find pixels in the 90% confidence region that will be visible to HET in the next 24 hours
        #NOTE: THIS MUTATES THE SKYMAPS
        skymap1 = np.copy(skymap)
        _, frac_visible_MMT = prob_observable(
            skymap, MMT_observatory, current_time, alert.directory
        )
        _, frac_visible_magellan = prob_observable(
            skymap1, magellan_observatory, current_time, alert.directory
        )
        
        MMT_observatory.skymap = skymap
        magellan_observatory.skymap = skymap1
        
        magellan_observatory.frac_visible = frac_visible_magellan
        MMT_observatory.frac_visible = frac_visible_MMT
        
        print("Two-dimensionally, percentage of pixels visible to MMT: "+str(frac_visible_MMT*100)+"%")
        print("Two-dimensionally, percentage of pixels visible to Magellan: "+str(frac_visible_magellan*100)+"%")
        
        cattop_all, logptop_all = get_galaxies.write_catalogs(
            savedir = alert.directory,
            observatories = [MMT_observatory, magellan_observatory]
        )
        
        # will be visible to at least one telescope!
        """
        if HET_observatory.frac_visible <= 0.1:
            print("HET can't observe the source.")
            send_not_visible_info(HET_observatory, contact_lists_all['email'])
            return False
        """
        if len(cattop_all[0]) > 0:
            send_mapped_alert_info(
                MMT_observatory, contact_lists_all,
                contact_lists_followup, phase_II=False
            )
            
        if len(cattop_all[1]) > 0:
            send_mapped_alert_info(
                magellan_observatory, contact_lists_all,
                contact_lists_followup, phase_II=False
            )

        """
        get_LST.get_LST(
            savedir = alert.directory,
            targf = os.path.join(alert.directory,'HET_Visible_Galaxies_prob_list.dat')
        )
        """
        return True
            

            
if __name__ == "__main__":
    ###########Things start here####################
   
    SAVE_PATH = '../../'
    stream_start_pos = StartPosition.EARLIEST
    #print("Starting stream at "+str(stream_start_pos))
    stream = Stream(start_at=stream_start_pos)
    #stream = Stream()

    num_messages = 0

    print("Listening for Alerts from kafka")

    with stream.open("kafka://kafka.scimma.org/igwn.gwalert", "r") as s:
        for message in s:
            event = message.content[0]['event']
            message_content = message.content[0]
            alert_time = message_content['time_created']
            alert_time = Time(alert_time)
            
            # check alert within last 24 hours
            if alert_time.mjd < Time.now().mjd - 1.:
                continue
                
            print(alert_time)
            
            num_messages+=1
            #if num_messages > 2:
            #    sys.exit()

            '''send that fits file into process_fits'''
            if event is not None:
                print('Calling process_fits')
                process_fits(alert_message = message_content, save_path = SAVE_PATH, skip_test_alerts = True)

