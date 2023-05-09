# -*- coding: utf-8 -*-
# Copyright 2021-2022 Cardiff University
# Distributed under the terms of the BSD-3-Clause license

"""Tests for :mod:`igwn_auth_utils.requests`.
"""

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__credits__ = "Leo Singer <leo.singer@ligo.org>"

import os
import stat
from unittest import mock
from urllib.parse import urlencode

import pytest

from requests import RequestException

from .. import requests as igwn_requests
from ..error import IgwnAuthError
from .test_scitokens import rtoken  # noqa: F401


def _empty(*args, **kwargs):
    return []


def _igwnerror(*args, **kwargs):
    raise IgwnAuthError("error")


def mock_no_scitoken():
    return mock.patch(
        "igwn_auth_utils.scitokens._find_tokens",
        _empty,
    )


def mock_no_x509():
    return mock.patch(
        "igwn_auth_utils.requests.find_x509_credentials",
        _igwnerror,
    )


@pytest.fixture
def netrc(tmp_path):
    netrc = tmp_path / "netrc"
    netrc.write_text(
        "machine example.org login albert.einstein password super-secret",
    )
    netrc.chmod(stat.S_IRWXU)
    return netrc


def has_auth(session):
    return bool(
        session.auth
        or session.cert
        or "Authorization" in session.headers
    )


class TestSession:
    Session = igwn_requests.Session

    # -- SessionErrorMixin

    def test_raise_for_status_hook(self, requests_mock):
        # define a request that returns 404 (not found)
        requests_mock.get(
            "https://test.org",
            status_code=404,
            reason="not found",
        )

        # with the kwarg a RequestException is raised
        with pytest.raises(
            RequestException,
            match=r"404 Client Error: not found for url: https://test.org/",
        ):
            igwn_requests.get("https://test.org")

    # -- SessionAuthMixin

    def test_noauth_args(self):
        """Test that `Session(force_noauth=True, fail_if_noauth=True)`
        is invalid.
        """
        with pytest.raises(ValueError):
            self.Session(force_noauth=True, fail_if_noauth=True)

    def test_fail_if_noauth(self):
        """Test that `Session(fail_if_noauth=True)` raises an error
        """
        with pytest.raises(IgwnAuthError):
            self.Session(
                token=False,
                cert=False,
                url=None,
                fail_if_noauth=True,
            )

    def test_force_noauth(self):
        """Test that `Session(force_noauth=True)` overrides auth kwargs
        """
        sess = self.Session(cert="cert.pem", force_noauth=True)
        assert not has_auth(sess)

    @mock_no_scitoken()
    @mock_no_x509()
    def test_defaults(self):
        """Test that the `Session()` defaults work in a noauth environment
        """
        sess = self.Session()
        assert not has_auth(sess)

    # -- tokens

    def test_token_explicit(self, rtoken):  # noqa: F811
        """Test that tokens are handled properly
        """
        sess = self.Session(token=rtoken)
        assert sess.token is rtoken
        assert sess.headers["Authorization"] == (
            igwn_requests.scitoken_authorization_header(rtoken)
        )

    def test_token_serialized(self, rtoken):  # noqa: F811
        """Test that serialized tokens are handled properly
        """
        serialized = rtoken.serialize()
        sess = self.Session(token=serialized)
        assert sess.headers["Authorization"] == f"Bearer {serialized}"
        # will not deserialise a token for storage
        assert sess.token is None

    @mock.patch("igwn_auth_utils.requests.find_scitoken")
    def test_token_discovery(self, find_token, rtoken):  # noqa: F811
        find_token.return_value = rtoken
        sess = self.Session()
        assert sess.token is rtoken
        assert sess.headers["Authorization"] == (
            igwn_requests.scitoken_authorization_header(rtoken)
        )

    @mock.patch(
        "igwn_auth_utils.requests.find_scitoken",
        side_effect=IgwnAuthError,
    )
    def test_token_required_failure(self, _):
        with pytest.raises(IgwnAuthError):
            self.Session(token=True)

    @pytest.mark.parametrize(("url", "aud"), (
        ("https://secret.example.com:8008", [
            "https://secret.example.com",
            "ANY",
        ]),
        (None, None)
    ))
    @mock.patch("igwn_auth_utils.requests.find_scitoken")
    def test_token_audience_default(self, find_scitoken, url, aud):
        """Check that the default `token_audience` is set correctly.
        """
        self.Session(
            url=url,
            token=True,
            cert=False,
            auth=False,
        )
        find_scitoken.assert_called_once_with(aud, None)

    # -- X.509

    def test_cert_explicit(self):
        """Test that cert credentials are stored properly
        """
        sess = self.Session(token=False, cert="cert.pem")
        assert sess.cert == "cert.pem"
        assert has_auth(sess)

    @mock.patch(
        "igwn_auth_utils.requests.find_x509_credentials",
        return_value="test.pem",
    )
    def test_cert_discovery(self, _):
        """Test that automatic certificate discovery works
        """
        assert self.Session(token=False).cert == "test.pem"

    @mock.patch(
        "igwn_auth_utils.requests.find_x509_credentials",
        side_effect=IgwnAuthError,
    )
    def test_cert_required_failure(self, _):
        with pytest.raises(IgwnAuthError):
            self.Session(token=False, cert=True)

    # -- basic auth

    @mock.patch.dict(os.environ)
    @pytest.mark.parametrize(("url", "auth"), (
        ("https://example.org", ("albert.einstein", "super-secret")),
        ("https://bad.org", None),
    ))
    def test_basic_auth(self, netrc, url, auth):
        os.environ["NETRC"] = str(netrc)
        sess = self.Session(cert=False, token=False, url=url)
        assert sess.auth == auth

    # -- all

    @mock.patch("igwn_auth_utils.requests.find_scitoken")
    @pytest.mark.parametrize(("cert", "token", "auth"), (
        ("A", True, ("C", "D")),  # all
        (False, True, ("C", "D")),  # no cert
        ("A", False, ("C", "D")),  # no token
        ("A", True, False),  # no auth
        (False, False, ("C", "D")),  # no cert or token
    ))
    def test_multi_auth(
        self,
        find_token,
        rtoken,  # noqa: F811
        cert,
        token,
        auth,
    ):
        """Check that Session._init_auth records all auth options

        In case a remote host accepts X.509 but not tokens, but the user
        has a valid ANY token (for example).
        """
        find_token.return_value = rtoken
        sess = self.Session(cert=cert, token=token, auth=auth)
        assert sess.cert == (cert or None)
        if token:
            sess.headers["Authorization"].startswith("Bearer")
        assert sess.auth == (auth or None)


def test_get(requests_mock):
    """Test that `igwn_auth_utils.requests.get` can perform a simple request
    """
    requests_mock.get(
        "https://test.org",
        text="TEST",
    )
    assert igwn_requests.get("https://test.org").text == "TEST"


@mock.patch("igwn_auth_utils.requests.Session")
def test_get_session(mock_session):
    """Test that ``session`` for `igwn_auth_utils.requests.get` works
    """
    session = mock.MagicMock()
    assert igwn_requests.get("https://test.org", session=session)
    session.request.assert_called_once_with("get", "https://test.org")
    mock_session.assert_not_called()


def test_post(requests_mock):
    """Test that `igwn_auth_utils.requests.post` can perform a simple request.
    """
    data = {"a": 1, "b": 2}
    requests_mock.post(
        "https://example.com",
        text="THANKS",
    )
    # check that the correct response got passed through
    assert igwn_requests.post(
        "https://example.com",
        data=data,
    ).text == "THANKS"
    # check that the data was encoded into the request properly
    req = requests_mock.request_history[0]
    assert req.body == urlencode(data)
