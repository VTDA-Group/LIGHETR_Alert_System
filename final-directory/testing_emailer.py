from important_functions import *
import json

def get_email_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['email_list']

def email(contact_list_file_loc = 'contact_all_BNS.json', subject = None, body = None, files_to_attach = [], people_to_contact = []):
            
    print("testing email")
    
    email_list = get_email_list(file_loc = contact_list_file_loc)
    print(email_list)
    
    psu_email_receivers = ''
    non_psu_email_receivers = ''
    
    if len(people_to_contact) == 0:
        people_to_contact = email_list.keys
    
    print("emails of people to contact: "+str(people_to_contact))
    
    
    for i in range(len(people_to_contact)):
        person = people_to_contact[i]
        address = email_list[person]
        print("looking at address: "+str(address))
        if 'psu.edu' in address:
            
            print("adding email address: "+str(email_list[person]))
            psu_email_receivers+=str(email_list[person])+", "
        else:
            non_psu_email_receivers+=str(email_list[person])+", "
            
    if psu_email_receivers[-2:] == ', ':
        psu_email_receivers = psu_email_receivers[:-2]
        
    if non_psu_email_receivers[-2:] == ', ':
        non_psu_email_receivers = non_psu_email_receivers[:-2]

    print("psu email receivers: "+str(psu_email_receivers))

    all_email_recipients = [psu_email_receivers, non_psu_email_receivers]
    
    print('all email recipients: '+str(all_email_recipients))
        
    if subject is None:
        subject = "LIGO HET Followup Test Email"
    if body is None:
        body = """I think this works?
        Let's see
        I like cheese"""

    email_sender = 'hetligo@gmail.com'
    email_password = 'xadlutfmejxfyqeq'

    print("sending email")
    send_email(email_sender = email_sender, email_password = email_password, all_email_recipients = all_email_recipients, files = files_to_attach, subject = subject, body = body)





#main_text()
#email(files_to_attach = ['60027.74060185185/fits_plotted.png'])
