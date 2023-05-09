import os
import twilio
import time
from twilio.rest import Client
from twilio.twiml.voice_response import VoiceResponse, Say

# Set environment variables for your credentials
# Read more at http://twil.io/secure

calling_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}
texter_dict = {'Zhenyuan':'+18147772603', 'Kaylee':'+17863973538', 'OG':'+16789006318', 'Laura':'+14403615990', 'Mary':'+17034241176'}

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

def call_people(people_to_contact = [], from_ = "+16073886023", message_to_say = message_to_say, calling_dict = calling_dict):
    
    account_sid = "ACc430265c246c76afe3f2c2bc52fd7c8a"
    auth_token = "eb0888a1ffc6f472a26723b2edd2da39"


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

#call_people(reciever_dict = calling_dict, people_to_contact = ['Kaylee'])
