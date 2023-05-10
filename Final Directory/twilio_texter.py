import os
import twilio
import time
from twilio.rest import Client
import pdb

# Set environment variables for your credentials
# Read more at http://twil.io/secure

account_sid = os.environ['TWILIO_ACCOUNT_SID']
auth_token = os.environ['TWILIO_AUTH_TOKEN']
client = Client(account_sid, auth_token)

    

def send_text_messages(reciever_dict, people_to_contact = [], message_to_send = 'LIGHETR Alert Message.\nThis is a test', from_ = "+18333749011"):
    client = Client(account_sid, auth_token)
    
    if len(people_to_contact) == 0:
        people_to_contact = reciever_dict.keys()
    
    print("People to contact: "+str(people_to_contact))
    print("message_to_send: "+str(message_to_send))
    
    
    for person in people_to_contact:
        if person not in reciever_dict.keys():
            print('Reciever '+str(person)+' doesn\'t have a known phone number.')
            continue
            
        
        
        relevant_info = reciever_dict[person]
        print("relevantinfo: "+str(relevant_info))
        message = client.messages.create(
            to=relevant_info,
            from_=from_,
            body=message_to_send
        )
        

        print(message.sid)
        time.sleep(1)
        
