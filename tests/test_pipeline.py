import pytest

from lighetr_alert_system.run_HET import process_fits
import lighetr_alert_system.run_new_telescopes as rnt


def test_het_pipeline(dummy_alert_message):
    success = process_fits(
        dummy_alert_message,
        people_to_contact = ['Kaylee',],
    )
    assert success

    
def test_mmt_pipeline(dummy_alert_message):
    success = rnt.process_fits(
        dummy_alert_message,
        people_to_contact = ['Kaylee',],
    )
    assert success
    