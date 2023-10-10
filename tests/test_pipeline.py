import pytest

from lighetr_alert_system.run import process_fits

def test_full_pipeline(dummy_alert_message):
    success = process_fits(
        dummy_alert_message,
        people_to_contact = ['Kaylee',],
    )
    assert success
    