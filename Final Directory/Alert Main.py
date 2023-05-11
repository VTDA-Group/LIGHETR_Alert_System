import matplotlib.pyplot as plt
from astropy.time import Time
#import diagnose
#from diagnose import *
import json
import os
import pdb
from hop import Stream
from hop.io import StartPosition
import healpy as hp
from astropy.io import fits
import requests
from prob_obs import *
from get_galaxies import *
import get_galaxies
import get_LST
from twilio_caller import *
from twilio_texter import *
from testing_emailer import *
import time as TIme
import make_phaseii

recent_April = Time('2023-02-16T00:00:00.00')

def write_to_file(file_loc, data_out, separator = ' ', headers=None, append=False):
    '''inputs: file_loc-location to which to write to, data_to_write-this is a 2d array that will be written to a file
       outputs: none
       This function writes a 2d array to a file'''
       
    if not append:
        out = open(file_loc ,'w')
    else:
        out = open(file_loc ,'a+')
    
    if type(data_out) is str:
        out.write(data_out+"\n")
        
    else:
        for i  in range(len(data_out)):
            if isinstance(data_out[i], float) or type(data_out[i]) is str:
                out.write(str(data_out[i])+separator)
            else:
                for j in range(len(data_out[i])):
                    out.write(str(data_out[i][j])+separator)
                
            out.write("\n")
    out.close()


def get_email_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['email_list']
    
def get_caller_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['caller_list']
    
def get_texter_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['texter_list']

def process_fits(fits_file, alert_message = None):
        
        alert_type = alert_message['alert_type']
        print("alert type: "+str(alert_type))
        alert_time = alert_message['time_created']
        
        alert_time = Time(alert_time)
        
        if alert_time.jd < recent_April.jd:
            return
        
        #so we've found a time that we want to look at. I'll make a directory for this time.
        obs_time_dir = str(alert_time.mjd)+"/"
        #obs_time_dir = str(alert_time)+"/"
        
        if not os.path.exists(obs_time_dir):
            os.mkdir(obs_time_dir)
        superevent_id = alert_message['superevent_id']
        fits_url = 'https://gracedb.ligo.org/api/superevents/'+str(superevent_id)+'/files/bayestar.multiorder.fits'
        
        
        #download the multi-order fits file from the fits_url filelocation
        url = fits_url
        multiorder_file_name = obs_time_dir+'multiorder_fits_'+str(superevent_id)+'.fits'
        req = requests.get(url)
        file = open(multiorder_file_name, 'wb')
        for chunk in req.iter_content(100000):
            file.write(chunk)
        file.close()
        
        #flatten the multi-order fits file into single-order
        singleorder_file_name = obs_time_dir+'flattened_multiorder_fits_'+superevent_id+'.fits'
        os.system('ligo-skymap-flatten '+str(multiorder_file_name)+' '+singleorder_file_name+' --nside 256')
        
        #open the flattened fits file
        fits_file = fits.open(singleorder_file_name)
        
        #get the skymap and header of the flattened, single-order fits file
        skymap,header = hp.read_map(singleorder_file_name, h=True, verbose=False)
        
        
        #plot the sky-localization from the flattened, single-order fits file
        m1 = hp.read_map(singleorder_file_name)
        hp.mollview(m1, rot = (180, -10, 0))
        plt.savefig(obs_time_dir+'fits_plotted.png')
        plt.close()


        # Print and save some values from the FITS header.
        header = dict(header)
        
        
        obs_time = Time(header['DATE-OBS'],format='isot',scale='utc')
        time = Time.now()
        # This is our testing time
        #time = Time('2023-04-17T07:00:00.00')
        dist = str(header['DISTMEAN']) + ' +/- ' + str(header['DISTSTD'])
        header['id'] = superevent_id

        
        
        
        
        
        
        '''
        The format of these alerts is given in this website:
        https://emfollow.docs.ligo.org/userguide/content.html
        
        So, the probabilities of BBH, BNS, and so on is located in:
        alert_message['event']['classification']['BNS']
        alert_message['event']['classification']['BBH']
        alert_message['event']['classification']['NSBH']
        alert_message['event']['classification']['Noise']
        
        it's unclear whether the key for noise will be either 'Noise' or 'Terrestrial', so I'll make the code flexible enough to use both
        
        '''
        
        # Making a pie chart of the type of event for the email
        noise_or_terrestrial = 'Terrestrial'
        if noise_or_terrestrial not in  alert_message['event']['classification'].keys():
            noise_or_terrestrial = 'Noise'
        event_poss = ['BBH', 'BNS', 'NSBH', noise_or_terrestrial]
        labels = []
        sizes = []
        
        for label in event_poss:
            val = float(alert_message['event']['classification'][label])
            if val >= 0:
                sizes.append(val)
                labels.append(label)
        
        
        #labels = ['%s (%.1f %%)'%(lab,pct) for lab,pct in zip(labels,sizes)]
        fig1, ax1 = plt.subplots()
        patches,texts = ax1.pie(sizes, startangle=90)
        ax1.legend(patches, labels, loc="best")
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig(obs_time_dir+'piechart.png')

        
        #find probabilities and list of galaxies visible to HET
        timetill90, m, visible_mask, never_visible_mask = prob_observable(skymap, header, time, savedir = obs_time_dir, plot=True)

        alert_message['skymap_fits'] = singleorder_file_name
        alert_message['skymap_array'] = m
        alert_message['visible_mask'] = visible_mask
        alert_message['never_visible_mask'] = never_visible_mask
        
        #print("In alert main, m at never visible: "+str(m[visible_mask]))
        #print("In alert main, m at never visible: "+str(m[never_visible_mask]))
        
        
        #if False:
        if timetill90 ==-99:
            print("HET can't observe the source.")
            
            return
        else:
            cattop, logptop = get_galaxies.write_catalog(alert_message, savedir = obs_time_dir)
            
            if len(cattop) == 0:
                print("Only visible during the day, I think")
                write_to_file(obs_time_dir+" observability.txt", "This object is only visible during the day, I think.")
                return
                
            #print("Source has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))))
            #print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(int(round(100 * probfull))))
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90))
            #write_to_file(obs_time_dir+" observability.txt", "Source has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))), append = False)
            #write_to_file(obs_time_dir+" observability.txt", "Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(int(round(100 * probfull))), append = True)
            write_to_file(obs_time_dir+" observability.txt", '{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90), append = True)
            
            mincontour = get_LST.get_LST(savedir = obs_time_dir,targf = obs_time_dir+'HET_Visible_Galaxies_prob_list.dat')
            
            
            
            print("people to contact: "+str(people_to_contact))
            
            
            
            email_body = 'A Neutron Star Merger has been detected by LIGO.\n{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nSource has a {:.1f}% chance of being observable now.\n\nPlease join this zoom call: https://us06web.zoom.us/j/87536495694".format(int(round(100 * prob)))
            #email(contact_list_file_loc = contact_list_file_loc, subject='LIGHETR Alert: NS Merger Detected', body = email_body, files_to_attach = [], people_to_contact = people_to_contact)
            
            
            calling_dict = get_caller_list(contact_list_file_loc)
            texting_dict = get_texter_list(contact_list_file_loc)
            
            make_phaseii.make_phaseii(lstfile = obs_time_dir+'LSTs_Visible.out', savedir = obs_time_dir)
            
            #call_people(calling_dict = calling_dict, people_to_contact = people_to_contact, message_to_say = 'NS Event Detected. Check email for information.')
            
            
            
###########Things start here####################
contact_list_file_loc = 'contact_only_HET_BNS.json'
people_to_contact = ['Karthik', 'Srisurya']


#stream_start_pos = 1600
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
        alert_type = message_content['alert_type']
        alert_time = message_content['time_created']
        
        alert_time = Time(alert_time)
        
        num_messages+=1
        #if num_messages > 2:
        #    sys.exit()
        
        # Is this from testing?
        if alert_time.jd < recent_April.jd:
            num_messages+=1
            continue
        
                
        skymap = None
        if 'skymap' in message_content.keys():
            skymap = message_content['skymap']
        if event is not None and 'skymap' in event.keys():
            skymap = event['skymap']

        '''send that fits file into process_fits'''
        if skymap is not None and event is not None:
            print('Calling process_fits')
            process_fits(fits_file = skymap, alert_message = message_content)
