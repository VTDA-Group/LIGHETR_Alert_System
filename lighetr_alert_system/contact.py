import json
import dataclasses
import os
from dataclasses import dataclass, field
from typing import Dict, List
from email.message import EmailMessage
import ssl
import smtplib
import mimetypes
import healpy as hp
import os.path
import datetime
import numpy as np
from lighetr_alert_system import make_phaseii
from lighetr_alert_system.twilio_caller import *


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

            
@dataclass
class Email:
    """Helper class to send out LIGHETR emails.
    """
    
    subject: str = "LIGO HET Followup Test Email"
    body: str = """This is just a test."""
    recipients: List[str] = field(default_factory=list)
    sender_email: str = os.environ['LIGO_EMAIL']
    sender_password: str = os.environ['EMAIL_PASSWORD']
    files: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        pass
            
    def send_email(self):
        everyone = " , ".join(self.recipients)
        body = f"{self.body}\nList of all Recipients of this email: {everyone}"
        now = datetime.datetime.now()
        body += f"\nAttempted send at time: {now}"

        for recipient in self.recipients:
            em = EmailMessage()
            em['From'] = self.sender_email
            em['To'] = recipient
            em['Subject'] = self.subject
            em.set_content(body)

            #attaching files
            for path in self.files:
                if not os.path.isfile(path):
                    print("File: "+str(path)+" doesn't exist. Cannot attach it.")
                    continue
                ctype, encoding = mimetypes.guess_type(path)
                if ctype is None or encoding is not None:
                    # No guess could be made, or the file is encoded (compressed), so
                    # use a generic bag-of-bits type.
                    ctype = 'application/octet-stream'
                maintype, subtype = ctype.split('/', 1)

                with open(path, 'rb') as fp:
                    em.add_attachment(
                        fp.read(),
                        maintype = maintype,
                        subtype = subtype,
                        filename = path
                    )

            #sending email
            context = ssl.create_default_context()
            try:
                with smtplib.SMTP_SSL("smtp.gmail.com", 465, context = context) as smtp:
                    smtp.login(self.sender_email, self.sender_password)
                    smtp.sendmail(self.sender_email, recipient, em.as_string())
                return True
            except:
                print(f"Wasn't able to send email to: {recipient}")
                return False

    
    
def get_contact_lists(file_loc = 'contact_all_BNS.json', people_to_contact = None):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    jsonDict = {
        'email': jsonObject['email_list'],
        'call': jsonObject['caller_list'],
    }
    
    for k in jsonDict:
        # convert to list
        dict_k = jsonDict[k]
        if people_to_contact is None:
            jsonDict[k] = [dict_k[x] for x in dict_k]
        else:
            jsonDict[k] = [dict_k[x] for x in dict_k if x in people_to_contact]
    
    return jsonDict


def send_burst_info(alert, contact_list):
    subject = f'LIGHETR Alert: GW Burst Event Detected (No Optical Counterpart) Event: {alert.event_id}'
    body = 'A gravitational wave burst event was detected. This means there is no sky localization from LIGO. We should ignore this event.\n Happy days!'
    
    if alert.test_event:
        subject = '[TEST, Can Safely Disregard!] '+ subject
        body = '[TEST EVENT!]\n' + body

    try:
        #Let's NOT send emails for burst events. --> let's send out burst info emails >:) 
        email = Email(subject=subject, body=body, recipients=contact_list)
        email.send_email()

        with open(alert.overview_file, "a+") as data_out:
            data_out.write(f'This is a burst event: {alert.event_id} | {alert.alert_type} Alert\n')
    except:
        return False
        
    return True
                              
                              
def send_low_prob_info(alert, contact_list, reason="low_prob"):
    
    subject = f'({alert.event_id} {alert.alert_type}) LIGHETR Alert: GW Event Detected (No Optical Counterpart)'

    body = f'A gravitational wave event was detected.\n\n'

    for k in alert.event_dict:
        body += f'Probability of {k}: {alert.event_dict[k]}\n'

    body += f'Distance to object: {alert.dist_mu} +/- {alert.dist_std} Mpc\n\n'

    if reason == 'significance':
        body += 'We will ignore this event because of FAR and significance cuts.\n\n'
    elif reason == 'distance':
        body += f'We will ignore this event because of distance cuts. We\'re using a distance cut of mu-1sigma < {alert.max_dist} Mpc.\n\n'
    else:
        body += 'We will ignore this event because it is unlikely to have a meaningful optical counterpart.\n\n'

    body += f'This had a FAR of {alert.far}\nSignificance of event: {alert.significance}\n'
    body += f'Gracedb Site: {alert.gracedb_site}\n\nHappy days!'

    if alert.test_event:
        subject = '[TEST, Can Safely Disregard!] '+ subject
        body = '[TEST EVENT!]\n' + body
        
    try:
        email = Email(subject=subject, body=body, recipients=contact_list)
        email.send_email()

        with open(alert.overview_file, "a+") as data_out:
            data_out.write(body)
        return True
    except:
        return False
        
        
def send_not_visible_info(observatory, contact_list):
    alert = observatory.alert
    subject = f'({alert.event_id} {alert.alert_type}) LIGHETR Alert: NS Merger Detected (NOT VISIBLE TO {observatory.name}).'                             
    body = 'A Neutron Star Merger has been detected by LIGO. This event is not visible to HET.'
    body += f'Percentage of the 90% localization region of this event visible to {observatory.name}: {observatory.frac_visible}.'
    body += 'This is too small to be worth following up. This is a courtesy email stating that an event was detected by LIGO.\n\n'
                       
    for k in alert.event_dict:
        body += f'Probability of {k}: {alert.event_dict[k]}\n'
    
    body += f'Distance to object: {alert.dist_mu} +/- {alert.dist_std} Mpc\n\n'  
    body += f'This had a FAR of {alert.far}\nSignificance of event: {alert.significance}\n'
    body += f'Gracedb Site: {alert.gracedb_site}\n\nHappy days!'
                              
    if alert.test_event:
        subject = '[TEST, Can Safely Disregard!] '+ subject
        body = '[TEST EVENT!]\n' + body

    try:
        email = Email(
            subject=subject,
            body=body,
            recipients=contact_list,
            files=[
                os.path.join(
                    alert.directory, f"{observatory.name}_Full_Visibility.pdf"
                )
            ]
        )
        email.send_email()

        write_to_file(
            os.path.join(
                alert.directory, 'observability.txt',
            ),
            f"{observatory.name} can't observe this source."
        )
        with open(alert.overview_file, "a+") as data_out:
            data_out.write(body)
        return True
    except:
        return False
        
        

def send_mapped_alert_info(observatory, contact_dir, contact_dir_followup, phase_II = True):
    alert = observatory.alert
    
    #sending emails out to everybody about the alert.
    subject = f'({alert.event_id} {alert.alert_type}) LIGHETR Alert: NS Merger Detected.'
    body = 'A Neutron Star Merger has been detected by LIGO.\n'
    body += '\n\nI have attached a figure here, showing the 90% contour of the sky localization where LIGO found a merger.'
    body += f'The portion in bright green is not visible to {observatory.name} because of declination limitations or because of sun constraints.'
    body += f'The portion in the dimmer blue-green is visible to {observatory.name} tonight. The percentage of pixels that are visible to {observatory.name} is {round(observatory.frac_visible*100, 3)}'
    body += '\n\n=================\n'
    
    for k in alert.event_dict:
        body += f'Probability of {k}: {alert.event_dict[k]}\n'
        
    body += f'Distance to object: {alert.dist_mu} +/- {alert.dist_std} Mpc\n\n'  
    body += f'This had a FAR of {alert.far}\nSignificance of event: {alert.significance}\n'
    body += f'Gracedb Site: {alert.gracedb_site}\n\nHappy days!\n'
    body += 'Please join this zoom call: https://us06web.zoom.us/j/87536495694\n'
    
    body_II = f'Attached is a phase ii submission file to {observatory.name} with the list of galaxies we wish to observe, '
    body_II += f'along with a file with LSTs of when they are visible\nGracedb Site: {alert.gracedb_site}'
    body_II += '\n Happy days!'+'.\nThanks for your help! \n\nZoom call: https://us06web.zoom.us/j/87536495694'
    
    if alert.test_event:
        subject = '[TEST, Can Safely Disregard!] '+ subject
        body = '[TEST EVENT!]\n' + body
        body_II = '[TEST EVENT!]\n' + body_II

    email = Email(
        subject=subject,
        body=body,
        recipients=np.append(contact_dir['email'], contact_dir_followup['email']),
        files=[
            os.path.join(
                alert.directory, f"{observatory.name}_Full_Visibility.pdf"
            )
        ]
    )
    email.send_email()

    #write_str = '{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90_HET)
    write_str = f'\nPercentage of visible pixels to {observatory.name}: {round(observatory.frac_visible*100, 3)}%'
    write_to_file(
        os.path.join(
            alert.directory, 'observability.txt',
        ),
        write_str,
        append = True
    )
    
    with open(alert.overview_file, "a+") as data_out:
        data_out.write(body)


    #calling people
    calling_list = contact_dir['call']
    message_to_say = 'Lie Go Gravitational Wave Neutron Star Event Detected. Check email for information. '
    if alert.test_event:
        message_to_say = 'This is a Test Alert. ' + message_to_say
    message_to_say = message_to_say * 10
    call_people(calling_list = calling_list, message_to_say = message_to_say)

    if phase_II:
        #sending phaseII file out to astronomer at observatory
        make_phaseii.make_phaseii(
            lstfile = os.path.join(alert.directory,'LSTs_Visible.out'),
            savedir = alert.directory
        )

        email2 = Email(
            subject=subject,
            body=body_II,
            recipients=contact_dir_followup['email'],
            files=[
                os.path.join(
                    alert.directory, "submission_to_HET.tsl"
                ),
                os.path.join(
                    alert.directory, "LSTs_Visible.pdf"
                ),
            ]
        )
        email.send_email()