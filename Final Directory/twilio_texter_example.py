import os
import twilio
import time
from twilio.rest import Client
import pdb

# Set environment variables for your credentials
# Read more at http://twil.io/secure

calling_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}
texter_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}

    

def send_text_messages(reciever_dict, people_to_contact = [], message_to_send = 'LIGHETR Alert Message.\nStay Fresh', from_ = "+16073886023"):
    
    account_sid = 'ACc430265c246c76afe3f2c2bc52fd7c8a'
    auth_token = "6a01135684c7da0946063d1d99f2a7f0"
    client = Client(account_sid, auth_token)
    
    if len(people_to_contact) == 0:
        people_to_contact = reciever_dict.keys()
    
    print("People to contact: "+str(people_to_contact))
    print("message_to_send: "+str(message_to_send))
    
    
    for homie in people_to_contact:
        if homie not in reciever_dict.keys():
            print('Reciever '+str(homie)+' doesn\'t have a known phone number.')
            continue
            
            
            
        diggies = reciever_dict[homie]
        print("diggies: "+str(diggies))
        message = client.messages.create(
            to=diggies,
            from_=from_,
            body=message_to_send
        )

        print(message.sid)
        time.sleep(1)
        
#send_text_messages(reciever_dict = texter_dict, people_to_contact = ['Kaylee'])

