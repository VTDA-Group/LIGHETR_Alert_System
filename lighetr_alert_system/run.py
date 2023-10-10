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
from lighetr_alert_system import get_galaxies, get_LST



def process_fits(alert_message, people_to_contact = None, skip_test_alerts = True):
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
        contact_lists_all = get_contact_lists(
            file_loc = os.path.join(
                dir_path,
               "../data/contact_all_LIGO.json"
            ),
            people_to_contact = people_to_contact
        )
        contact_lists_HET = get_contact_lists(
            file_loc = os.path.join(
                dir_path,
                "../data/contact_only_HET_BNS.json"
            ),
            people_to_contact = people_to_contact
        )
        
        alert_type = alert_message['alert_type']
        
        superevent_id = alert_message['superevent_id']
        gracedb_site = 'https://gracedb.ligo.org/superevents/'+str(superevent_id)+'/'
        test_event = ( superevent_id[0] in ['M', 'T'] )
            
        if skip_test_alerts and test_event: #if we want to ignore test events, and this is a test event, ignore it.
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
            if sizes[3] < 0.9: #send the email if terrestrial signal < 0.9
                send_low_prob_info(alert, contact_lists_all['email'])
            else:
                print("Probably terrestrial")
            return False
        
        if alert.dist_mu - alert.dist_std > alert.max_dist:
            send_low_prob_info(alert, contact_lists_all['email'], reason="distance")     
            return False
        
        current_time = Time.now()
        print("Calculating probabilities")
        
        HET_observatory = Observatory(
            name='HET',
            dec_range=(-12.0, 74.0),
            loc=(-104.01472,30.6814,2025),
            pupil=np.loadtxt(
                os.path.join(
                    dir_path,
                    '../data/hetpix.dat'
                )
            ),
            alert=alert
        )
        
        gen_observatory = Observatory(
            name='GEN',
            dec_range=(-90.0, 90.0),
            loc=(-104.01472,30.6814,2025),
            alert=alert,
        )
        
        #find pixels in the 90% confidence region that will be visible to HET in the next 24 hours
        #NOTE: THIS MUTATES THE SKYMAPS
        skymap1 = np.copy(skymap)
        _, frac_visible_gen = prob_observable(
            skymap, gen_observatory, current_time, alert.directory
        )
        _, frac_visible_HET = prob_observable(
            skymap1, HET_observatory, current_time, alert.directory
        )
        
        gen_observatory.skymap = skymap
        HET_observatory.skymap = skymap1
        
        gen_observatory.frac_visible = frac_visible_gen
        HET_observatory.frac_visible = frac_visible_HET
        
        print("Two-dimensionally, percentage of pixels visible to HET: "+str(frac_visible_HET*100)+"%")
        
        cattop_all, logptop_all = get_galaxies.write_catalogs(
            savedir = alert.directory,
            observatories = [gen_observatory, HET_observatory]
        )
        
        if HET_observatory.frac_visible <= 0.1:
            print("HET can't observe the source.")
            send_not_visible_info(HET_observatory, contact_lists_all['email'])
            return False
        
        else:
            cattop, logptop = cattop_all[1], logptop_all[1]
            
            if len(cattop) == 0:
                return False
                
            print("Percentage of visible pixels to HET: "+str(round(frac_visible_HET*100, 3))+"%")
            
            get_LST.get_LST(
                savedir = alert.directory,
                targf = os.path.join(alert.directory,'HET_Visible_Galaxies_prob_list.dat')
            )
            
            send_mapped_alert_info(HET_observatory, contact_lists_all, contact_lists_HET)
            return True
            

            
if __name__ == "__main__":
    ###########Things start here####################
   
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
                process_fits(alert_message = message_content, skip_test_alerts = True)
