import pytest

from lighetr_alert_system.contact import (
    get_contact_lists,
    Email,
    send_burst_info,
    send_low_prob_info,
    send_not_visible_info,
)
from lighetr_alert_system.twilio_caller import call_people
from lighetr_alert_system.observatory import Observatory

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
    """Send test phone call."""
    contact_lists = get_contact_lists(
        file_loc = test_contact_dir,
    )
    calling_list = contact_lists['call']
    message_to_say = 'This is a pie test call. '
    message_to_say = message_to_say * 3
    call_people(calling_list = calling_list, message_to_say = message_to_say)
    assert True
    
    
def test_send_burst_info(dummy_alert, test_contact_dir):
    """Send burst info email."""
    contact_list = get_contact_lists(
        file_loc = test_contact_dir,
    )['email']
    
    assert send_burst_info(dummy_alert, test_contact_dir)
    
    
def test_send_low_prob_info(dummy_alert, test_contact_dir):
    """Send emails from low BNS/NSBH probability or distance constraints."""
    contact_list = get_contact_lists(
        file_loc = test_contact_dir,
    )['email']
    dummy_alert.read_in_skymap()
    dummy_alert.flatten_skymap()
    _ = dummy_alert.unload_skymap()
    
    assert send_low_prob_info(dummy_alert, contact_list, reason="significance")
    assert send_low_prob_info(dummy_alert, contact_list, reason="low_prob")
    assert send_low_prob_info(dummy_alert, contact_list, reason="distance")
    
    
def test_send_not_visible_info(dummy_alert, het_loc, test_contact_dir):
    """Send emails for sources which are not visible by an observatory."""
    contact_list = get_contact_lists(
        file_loc = test_contact_dir,
    )['email']
    
    dummy_alert.read_in_skymap()
    dummy_alert.flatten_skymap()
    skymap = dummy_alert.unload_skymap()
    
    full_sky_obs = Observatory(
        "GEN",
        (-90., 90.),
        het_loc,
        alert=dummy_alert,
    )
    full_sky_obs.frac_visible = 0.01
    assert send_not_visible_info(full_sky_obs, contact_list)
    
    
    
    