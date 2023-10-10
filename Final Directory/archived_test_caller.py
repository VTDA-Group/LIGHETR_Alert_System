import os
import twilio
import time
from twilio.rest import Client
from twilio.twiml.voice_response import VoiceResponse, Say
import json
from twilio_caller import *

# Your Account Sid and Auth Token from twilio.com/user/account
# To set up environmental variables, see http://twil.io/secure
account_sid = os.environ['TWILIO_ACCOUNT_SID']
auth_token = os.environ['TWILIO_AUTH_TOKEN']
client = Client(account_sid, auth_token)


message_to_say = 'Neutron Star Merger Event Detected by LIGO. Check email for information.'

def get_caller_list(file_loc = 'contact_tester.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['caller_list']
    
def test_calling(file_loc = 'contact_tester.json', people_to_contact = []):
    
    calling_dict = get_caller_list(file_loc = file_loc)
    call_people(calling_dict = calling_dict, people_to_contact = people_to_contact, message_to_say = message_to_say)

test_calling(people_to_contact = ['Ashley'])
