from important_functions import *
import json
import os

ligo_list_email = os.environ['LIGO_EMAIL']
ligo_list_password = os.environ['EMAIL_PASSWORD']


def get_email_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['email_list']

def email(contact_list_file_loc = 'contact_all_BNS.json', subject = None, body = None, files_to_attach = [], people_to_contact = []):
    
    email_list = get_email_list(file_loc = contact_list_file_loc)
    print(email_list)
    
    psu_email_receivers = ''
    texas_email_recievers = ''
    non_academic_email_receivers = ''
    
    if len(people_to_contact) == 0:
        people_to_contact = email_list.keys()
    
    
    print("people to contact: "+str(people_to_contact))
    all_email_recipients = []
    for person in people_to_contact:
        if person not in email_list.keys():
            print("No email for person: "+str(person))
            continue
        address = email_list[person]
        all_email_recipients.append(address)
    
    print('all email recipients: '+str(all_email_recipients))
        
    if subject is None:
        subject = "LIGO HET Followup Test Email"
    if body is None:
        body = """This is just a test."""

    email_sender = ligo_list_email
    email_password = ligo_list_password
    send_email(email_sender = email_sender, email_password = email_password, all_email_recipients = all_email_recipients, files = files_to_attach, subject = subject, body = body)

#email(contact_list_file_loc = 'contact_only_HET_BNS.json', subject = "IGNORE RECENT LIGO ALERT", body = 'Hello Everyone!\nThe recent alert that just got sent out was an update for the event S230529ay, which was detected back in May. We will not be following up on this event. Please disregard the alert regarding the event S230529ay.\nHave a wonderful day!\n~LIGHETR Alert Team')
