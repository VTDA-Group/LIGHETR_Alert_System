import os
import twilio
import time
from twilio.rest import Client
from twilio.twiml.voice_response import VoiceResponse, Say
import json



account_sid = os.environ['TWILIO_ACCOUNT_SID']
auth_token = os.environ['TWILIO_AUTH_TOKEN']
client = Client(account_sid, auth_token)



calling_dict = {'Ashley':'+17726438132'}
texter_dict = {'Ashley':'+17726438132'}

message_to_say = 'We got one boys. There be a neutron star merger. Check your email for the deets.'



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

def call_people(people_to_contact = [], from_ = "+18333749011", message_to_say = message_to_say, calling_dict = calling_dict):
    client = Client(account_sid, auth_token)
    
    if len(people_to_contact) == 0:
        people_to_contact = calling_dict.keys()
        
    for homie in people_to_contact:

        diggies = calling_dict[homie]
        call = client.calls.create(
          twiml=build_message_to_say(message_to_say = message_to_say),
          #twiml='message_to_speak.xml',
          
          
          #url = 'http://demo.twilio.com/docs/voice.xml',
          to=diggies,
          from_=from_,
        )


        print(call.sid)
        time.sleep(1)

call_people(people_to_contact = ['Ashley'])
