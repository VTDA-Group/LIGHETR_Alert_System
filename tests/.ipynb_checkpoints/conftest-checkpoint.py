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
        "BBH": 0.5,
        "BNS": 0.1,
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