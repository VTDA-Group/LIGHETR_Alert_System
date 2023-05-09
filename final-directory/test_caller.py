import os
import twilio
import time
from twilio.rest import Client
from twilio.twiml.voice_response import VoiceResponse, Say
import json
from twilio_caller_example import *

# Set environment variables for your credentials
# Read more at http://twil.io/secure

calling_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}
texter_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}

message_to_say = 'Neutron Star Merger Event Detected by LIGO. Check email for information.'

def get_caller_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['caller_list']
    
def test_calling(file_loc = 'contact_all_BNS.json', people_to_contact = []):
    
    calling_dict = get_caller_list(file_loc = file_loc)
    call_people(calling_dict = calling_dict, people_to_contact = people_to_contact, message_to_say = message_to_say)

test_calling(people_to_contact = ['Karthik'])
