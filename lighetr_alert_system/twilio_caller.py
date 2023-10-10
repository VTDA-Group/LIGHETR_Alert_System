import os
import twilio
import time
from twilio.rest import Client
from twilio.twiml.voice_response import VoiceResponse, Say
import json

account_sid = os.environ['TWILIO_ACCOUNT_SID']
auth_token = os.environ['TWILIO_AUTH_TOKEN']
client = Client(account_sid, auth_token)

#calling_list = ['+17726438132'} # Ashley
#texter_dict = {'+17726438132'} # Ashley

message_to_say = 'There is a BNS!.'

voice = 'Polly.Geraint'

def build_message_to_say(voice = voice, message_to_say = message_to_say):
    xml = ''
    if len(voice) == 0:
        xml = '<Response><Say>'
    else:
        xml = '<Response><Say voice = \"'+str(voice)+'\">'
    xml += message_to_say
    xml += '</Say>'
    xml += '</Response>'
    
    print("going to say: "+str(xml))
    
    return xml

def call_people(calling_list, from_ = "+18333749011", message_to_say = message_to_say):
    client = Client(account_sid, auth_token)

    for diggies in calling_list:
        
        call = client.calls.create(
          twiml=build_message_to_say(message_to_say = message_to_say),
          #twiml='message_to_speak.xml',
          #url = 'http://demo.twilio.com/docs/voice.xml',
          to=diggies,
          from_=from_,
        )

        print(call.sid)
        time.sleep(1)

# This is just a test
# call_people(people_to_contact = ['Ashley'])
