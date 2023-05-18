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
from prob_obs import *
from get_galaxies import *
import get_galaxies
import get_LST
from twilio_caller import *
from twilio_texter import *
from testing_emailer import *
import time as TIme
import make_phaseii
import urllib, json
import urllib.request, json
from bs4 import BeautifulSoup
import requests

def find_files(url):

    soup = BeautifulSoup(requests.get(url).text)

    for a in soup.find_all('a'):
        yield a['href']


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

def file_on_website_exists(file):

    return False
    '''
    import requests
    r = requests.head(file, allow_redirects = True)
    if r.status_code == 200:
        return True
    return False'''

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
    
    
def process_fits_1(superevent_id, skip_test_alerts = False, test_event = False, fits_url = 'bayestar.multiorder.fits'):
        
        obs_time_dir = superevent_id+"/"
        
        print("inside process fits 2")
        
        
        
        
        #print(data)
        #print("data: " +str(data))
        
        
        #get the alert json file
        #get the real website where the json file may exist
        
        print("file existing: "+str(file_on_website_exists("https://gracedb.ligo.org/api/superevents/S230518h/files/S230518h-preliminary.json")))
        
        site = "https://gracedb.ligo.org/api/superevents/"+superevent_id+"/files/"
        
        #check if the prelim json exists
        found_file = False
        possible_filenames = ["preliminary", "update"]
        for filename in possible_filenames:
            file_to_try = site+superevent_id+"-"+filename+".json"
            print("trying: "+str(file_to_try))
            if file_on_website_exists(file_to_try):
                found_file = True
                url = "https://gracedb.ligo.org/api/superevents/"+superevent_id+"/files/"+superevent_id+"-"+str(filename)+".json"

        if not found_file:
            url = "https://gracedb.ligo.org/api/superevents/"+superevent_id+"/files/"+superevent_id+"-preliminary.json"
            print("Cannot automatically find the json file on the website. Will attempt the default: "+str(url))
        
        with urllib.request.urlopen(url) as url:
            alert_message = json.load(url)
                
        
        print("alert message: "+str(alert_message))
        

        fits_url = 'https://gracedb.ligo.org/api/superevents/'+str(superevent_id)+'/files/'+fits_url
        
        if not os.path.exists(obs_time_dir):
            os.mkdir(obs_time_dir)
        

        
        #download the multi-order fits file from the fits_url file location
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
        
        #get the skymap and header of the flattened, single-order fits file
        skymap,header = hp.read_map(singleorder_file_name, h=True, verbose=False)
        
        
        #plot the sky-localization from the flattened, single-order fits file
        os.system("ligo-skymap-plot %s -o %s" % (singleorder_file_name, obs_time_dir+"fits_plotted.png"))


        # Print and save some values from the FITS header.
        header = dict(header)
        obs_time = Time(header['DATE-OBS'],format='isot',scale='utc')
        time = Time.now()
        dist = str(header['DISTMEAN']) + ' +/- ' + str(header['DISTSTD'])
        header['id'] = superevent_id
        
        
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
        
        #if it's not at least 30% a (BNS or NSBH) signal, ignore it.
        if (sizes[1]+sizes[2])/(sizes[0]+sizes[1]+sizes[2]+sizes[3]) < 0.3:
            return

        
        
        
        #find pixels in the 90% confidence region that will be visible to HET in the next 24 hours
        timetill90, m, frac_visible = prob_observable(skymap, header, time, savedir = obs_time_dir, plot=True)

        print("Two-dimensionally, percentage of pixels visible to HET: "+str(frac_visible*100)+"%")

        #alert_message['skymap_fits'] = singleorder_file_name
        #alert_message['skymap_array'] = m
        
        #if False:
        if timetill90 ==-99 or frac_visible <= 0.0:
            print("HET can't observe the source.")
            write_to_file(obs_time_dir+" observability.txt", "HET can't observe this source.")
            return
        else:
            cattop, logptop = get_galaxies.write_catalog(alert_message, savedir = obs_time_dir)
            
            if len(cattop) == 0:
                return
                
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nPercentage of visible pixels to HET: "+str(round(frac_visible*100, 3))+"%")
            write_to_file(obs_time_dir+" observability.txt", '{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nPercentage of visible pixels to HET: "+str(round(frac_visible*100, 3))+"%", append = True)
            mincontour = get_LST.get_LST(savedir = obs_time_dir,targf = obs_time_dir+'HET_Visible_Galaxies_prob_list.dat')
            
            
            
            #sending emails out to everybody about the alert.
            email_subject = 'LIGHETR Alert: NS Merger Detected'
            email_body = 'A Neutron Star Merger has been detected by LIGO.\n{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nI have attached a figure here, showing the 90% contour of the sky localization where LIGO found a merger. The portion in bright green is not visible to HET because of declination limitations or because of sun constraints. The portion in the dimmer blue-green is visible to HET tonight. The percentage of pixels that are visible to HET is "+str(round(frac_visible*100, 3))+"% \n\nPlease join this zoom call: https://us06web.zoom.us/j/87536495694"
            if test_event:
                email_subject = '[TEST, Can Safely Disregard!] '+email_subject
                email_body = '[TEST EVENT!]' + email_body
                
                
            email(contact_list_file_loc = contact_list_file_loc, subject=email_subject, body = email_body, files_to_attach = [obs_time_dir+"HET_Full_Visibility.pdf"], people_to_contact = people_to_contact)
            
            
            #calling people
            calling_dict = get_caller_list(contact_list_file_loc)
            message_to_say = 'Lie Go Gravitational Wave Neutron Star Event Detected. Check email for information. '
            if test_event:
                message_to_say = 'This is a Test Alert. ' + message_to_say
            message_to_say = message_to_say * 10
            call_people(calling_dict = calling_dict, people_to_contact = people_to_contact, message_to_say = message_to_say)
            
            
            
            #sending phaseII file out to astronomer at HET
            make_phaseii.make_phaseii(lstfile = obs_time_dir+'LSTs_Visible.out', savedir = obs_time_dir)
            
            email_body = 'Attached is a phase ii submission file to HET with the list of galaxies we wish to observe, along with a file with LSTs of when they are visible.\nThanks for your help! \n\nZoom call: https://us06web.zoom.us/j/87536495694'
            subject='LIGHETR Alert: NS Merger Detected'
            
            if test_event:
                email_body = '[TEST, Can Safely Disregard!] ' + email_body
                subject = '[TEST, Can Safely Disregard!] ' + subject
            if "HET" in people_to_contact:

    
                email(contact_list_file_loc = contact_list_file_loc, subject = subject, body = email_body, files_to_attach = [obs_time_dir+"submission_to_HET.tsl", obs_time_dir+"LSTs_Visible.pdf"], people_to_contact = ['HET'])
            
            
###########Things start here####################
contact_list_file_loc = 'contact_only_HET_BNS.json'
#people_to_contact = ["Karthik", "Srisurya", "HET"]
people_to_contact = ["Karthik"]

print("Running")

process_fits_1(superevent_id = "S230518h")

print("did process fits 2")
