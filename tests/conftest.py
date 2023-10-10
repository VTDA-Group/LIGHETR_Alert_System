import os

import numpy as np
import pandas as pd
import pytest

from lighetr_alert_system.alert import Alert
from astropy.time import Time

TEST_DIR = os.path.dirname(__file__)
DATA_DIR_NAME = "data"

# pylint: disable=missing-function-docstring, redefined-outer-name


@pytest.fixture
def test_data_dir():
    return os.path.join(TEST_DIR, DATA_DIR_NAME)

@pytest.fixture
def dummy_event_dict():
    return {
        "BBH": 0.1,
        "BNS": 0.5,
        "NSBH": 0.1,
        "terrestrial": 0.3
    }

@pytest.fixture
def dummy_alert(tmp_path, dummy_event_dict):
    return Alert(
        event_id='S230810af',
        time=Time.now(),
        alert_type='S',
        save_path=tmp_path,
        event_dict=dummy_event_dict,
    )


@pytest.fixture
def midnight_time():
    # midnight in Texas time
    time = '2023-01-01T08:00:00.0'
    return Time(time, format='isot', scale='utc')


@pytest.fixture
def het_loc():
    return (-104.01472,30.6814,2025)


@pytest.fixture
def test_contact_dir(test_data_dir):
    return os.path.join(test_data_dir, 'contact_tester.json')


@pytest.fixture
def dummy_alert_message(dummy_event_dict):
    alert_message = {}
    alert_message['alert_type'] = 'PRELIMINARY'
    alert_message['time_created'] = Time.now()
    alert_message['superevent_id'] = 'S230810af'
    alert_message['event'] = {
        'classification': dummy_event_dict,
        'far': 0.0,
        'significant': True
    }
    return alert_message
    