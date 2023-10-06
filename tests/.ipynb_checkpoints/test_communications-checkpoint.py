import pytest

from lighetr_alert_system.contact import (
    get_contact_lists,
    Email
)
from lighetr_alert_system.twilio_caller import call_people

def test_get_contact_lists(test_contact_dir):
    contact_lists = get_contact_lists(
        file_loc = test_contact_dir,
        people_to_contact = []
    )
    
    assert 'email' in contact_lists
    assert 'call' in contact_lists
    assert len(contact_lists['email']) == 0
    
    contact_lists2 = get_contact_lists(
        file_loc = test_contact_dir,
    )
    
    assert 'email' in contact_lists2
    assert 'call' in contact_lists2
    assert len(contact_lists2['email']) == 1
    
    
def test_send_email(test_contact_dir):
    """Send test email."""
    contact_lists = get_contact_lists(
        file_loc = test_contact_dir,
    )
    test_email = Email(
        subject = "LIGHETR PYTEST EMAIL",
        body = """This is just a test.""",
        recipients = contact_lists['email']
    )
    
    assert test_email.send_email()
    
    
def test_call(test_contact_dir):
    """Send test email."""
    contact_lists = get_contact_lists(
        file_loc = test_contact_dir,
    )
    calling_list = contact_lists['call']
    message_to_say = 'This is a pie test call. '
    message_to_say = message_to_say * 3
    call_people(calling_list = calling_list, message_to_say = message_to_say)
    assert True
    