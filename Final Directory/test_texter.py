from twilio_texter import *
import json

message = 'Neutron Star Merger Event Detected by Ligo. Check email for information.'

def get_texter_list(file_loc = 'contact_tester.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['texter_list']
    
def test_texting(file_loc = 'contact_tester.json', people_to_contact = [], message_to_send = message):
    
    texting_dict = get_texter_list(file_loc = file_loc)
    send_text_messages(reciever_dict = texting_dict, people_to_contact = people_to_contact, message_to_send = message_to_send)

test_texting(people_to_contact = [], message_to_send = 'LIGO Alert System Test.')
