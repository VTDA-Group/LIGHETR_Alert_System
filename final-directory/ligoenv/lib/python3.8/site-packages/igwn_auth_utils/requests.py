# -*- coding: utf-8 -*-
# Copyright 2021-2022 Cardiff University
# Distributed under the terms of the BSD-3-Clause license

"""Python Requests interface with IGWN authentication

This is heavily inspired by Leo Singer's excellent
:mod:`requests_gracedb` package.
"""

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__credits__ = "Leo Singer <leo.singer@ligo.org>"

from urllib.parse import urlparse

from safe_netrc import netrc

import requests

from scitokens import SciToken

from .error import IgwnAuthError
from .scitokens import (
    find_token as find_scitoken,
    target_audience as scitoken_audience,
    token_authorization_header as scitoken_authorization_header,
)
from .x509 import (
    find_credentials as find_x509_credentials,
)

_auth_session_parameters = """
    Discovery/configuration of authorisation/authentication methods
    is attempted in the following order:

    1.  if ``force_noauth=True`` is given, no auth is configured;

    2.  for SciTokens:

        1.  if a bearer token is provided via the ``token`` keyword argument,
            then use that, or

        2.  look for a bearer token by passing the ``token_audience``
            and ``token_scope`` keyword parameters to
            :func:`igwn_auth_utils.find_scitokens`;

    3.  for X.509 credentials:

        1.  if an X.509 credential path is provided via the ``cert`` keyword
            argument, then use that, or

        2.  look for an X.509 credential using
            :func:`igwn_auth_utils.find_x509_credential`

    4.  for basic auth (username/password):

        1.  if ``auth`` keyword is provided, then use that, or

        2.  read the netrc file located at :file:`~/.netrc`, or at the path
            stored in the :envvar:`$NETRC` environment variable, and look
            for a username and password matching the hostname given in the
            ``url`` keyword argument;

    5.  if none of the above yield a credential, and ``fail_if_noauth=True``
        was provided, raise a `ValueError`.

    Steps 2-4 are all tried independently, with all valid credentials
    (one per type) configured for the session.
    It is up to the request receiver to handle the multiple credential
    types and prioritise between them.

    Parameters
    ----------
    token : `scitokens.SciToken`, `str`, `bool`, optional
        Bearer token (scitoken) input, one of

        - a bearer token (`scitokens.SciToken`),
        - a serialised token (`str`, `bytes`),
        - `False`: disable using tokens completely
        - `True`: discover a valid token via
          :func:`igwn_auth_utils.find_scitoken` and
          error if something goes wrong
        - `None`: try and discover a valid token, but
          try something else if that fails

    token_audience : `str`, list` of `str`
        The value(s) of the audience (``aud``) claim to pass to
        :func:`igwn_auth_utils.find_scitoken` when discovering
        available tokens.

    token_scope : `str`
        The value(s) of the ``scope``  ``audience`` and ``scope`` to pass to
        :func:`igwn_auth_utils.find_scitoken` when discovering
        available tokens.

    cert : `str`, `tuple`, `bool`, optional
        X.509 credential input, one of

        - path to a PEM-format certificate file,
        - a ``(cert, key)`` `tuple`,
        - `False`: disable using X.509 completely
        - `True`: discover a valid cert via
          :func:`igwn_auth_utils.find_x509_credentials` and
          error if something goes wrong
        - `None`: try and discover a valid cert, but
          try something else if that fails

    auth :  `tuple`, `object`, optional
        ``(username, password)`` `tuple` or other authentication/authorization
        object to attach to a `~requests.Request`

    url : `str`, optional
        the URL/host that will be queried within this session; this is used
        to set the default ``token_audience`` and to access credentials
        via :mod:`safe_netrc`.

    force_noauth : `bool`, optional
        Disable the use of any authorisation credentials (mainly for testing).

    fail_if_noauth : `bool`, optional
        Raise a `~igwn_auth_utils.IgwnAuthError` if no authorisation
        credentials are presented or discovered.

    Raises
    ------
    ~igwn_auth_utils.IgwnAuthError
        If ``cert=True`` or ``token=True`` is given and the relevant
        credential was not actually discovered, or
        if ``fail_if_noauth=True`` is given and no authorisation
        token/credentials of any valid type are presented or discovered.

    See also
    --------
    requests.Session
        for details of the standard options

    igwn_auth_utils.find_scitoken
        for details of the SciToken discovery

    igwn_auth_utils.find_x509_credentials
        for details of the X.509 credential discovery
    """.strip()


def _find_cred(func, *args, error=True, **kwargs):
    """Find a credential and maybe ignore an `~igwn_auth_utils.IgwnAuthError`

    This is an internal utility for the `SessionAuthMixin._init_auth`
    method which shouldn't necessary fail if it doesn't
    find a credential of any one type, but should just move on to the
    next option.
    """
    try:
        return func(*args, **kwargs)
    except IgwnAuthError:
        if error:
            raise
        return


def _hook_raise_for_status(response, *args, **kwargs):
    """Response hook to raise exception for any HTTP error (status >= 400)

    Reproduced (with permission) from :mod:`requests_gracedb.errors`,
    authored by Leo Singer.
    """
    return response.raise_for_status()


class SessionErrorMixin:
    """A mixin for :class:`requests.Session` to raise exceptions for HTTP
    errors.

    Reproduced (with permission) from :mod:`requests_gracedb.errors`,
    authored by Leo Singer.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hooks.setdefault("response", []).append(
            _hook_raise_for_status,
        )


class SessionAuthMixin:
    """Mixin for :class:`requests.Session` to add support for IGWN auth.

    By default this mixin will automatically attempt to discover/configure
    a bearer token (scitoken) or an X.509 credential, with options to
    require/disable either of those, or all authentication entirely.

    {parameters}
    """
    def __init__(
        self,
        token=None,
        token_audience=None,
        token_scope=None,
        cert=None,
        auth=None,
        url=None,
        force_noauth=False,
        fail_if_noauth=False,
        **kwargs,
    ):
        # initialise session and new attributes
        super().__init__(**kwargs)

        #: The SciToken object to be serialised and sent with all requests,
        #: this is only populated if a `~scitokens.SciToken` object is
        #: passed directly, or discovered automatically (i.e. a serialised
        #: token will not be deserialized and stored).
        self.token = None

        # handle options
        if force_noauth and fail_if_noauth:
            raise ValueError(
                "cannot select both force_noauth and fail_if_noauth",
            )
        if force_noauth:
            return

        # find creds if we can
        _auth = self._init_auth(
            auth=auth,
            cert=cert,
            token=token,
            token_audience=token_audience,
            token_scope=token_scope,
            url=url,
        )

        # if no auth was found, and we need it, fail here
        if not _auth and fail_if_noauth:
            raise IgwnAuthError("no valid authorisation credentials found")

    def _init_auth(
        self,
        auth=None,
        cert=None,
        token=None,
        token_audience=None,
        token_scope=None,
        url=None,
    ):
        # bearer token (scitoken)
        if token in (None, True):
            token = self._find_token(
                token_audience,
                token_scope,
                url=url,
                error=bool(token),
            )
        if isinstance(token, SciToken):
            self.token = token
        if token:
            token = self._set_token_header(token)

        # cert auth
        if cert in (None, True):  # not disabled and not given explicitly
            cert = self._find_x509_credentials(error=bool(cert))
        if cert:
            self.cert = cert

        # basic auth (netrc)
        if auth in (None, True):  # not disabled and not given explicitly
            auth = self._find_username_password(url)
        if auth:
            self.auth = auth

        return token or self.cert or self.auth

    def _set_token_header(self, token):
        """Serialise a `scitokens.SciToken` and format and store as an
        Authorization header for this session.

        Parameters
        ----------
        token : `scitokens.SciToken`, `str`, `bytes`
            the token to serialize, or an already serialized representation
        """
        if isinstance(token, (str, bytes)):  # load a valid token
            header = f"Bearer {token}"
        else:
            header = scitoken_authorization_header(token)
        self.headers["Authorization"] = header
        return header

    @staticmethod
    def _find_x509_credentials(error=True):
        """Find an X.509 certificate for authorization
        """
        return _find_cred(find_x509_credentials, error=error)

    @staticmethod
    def _find_token(audience, scope, url=None, error=True):
        """Find a bearer token for authorization
        """
        if audience is None and url is not None:
            audience = scitoken_audience(url, include_any=True)
        return _find_cred(find_scitoken, audience, scope, error=error)

    @staticmethod
    def _find_username_password(url):
        """Use `safe_netrc.netrc` to find the username/password for basic auth
        """
        host = urlparse(url).hostname

        try:
            result = netrc().authenticators(host)
        except IOError:
            return

        try:
            username, _, password = result
        except TypeError:
            return
        return username, password


class Session(
    SessionAuthMixin,
    SessionErrorMixin,
    requests.Session,
):
    """`requests.Session` class with default IGWN authorization handling

    {parameters}

    Examples
    --------
    To use the default authorisation discovery:

    >>> from igwn_auth_utils.requests import Session
    >>> with Session() as sess:
    ...     sess.get("https://science.example.com/api/important/data")

    To explicitly pass a specific :class:`~scitokens.SciToken` as the token:

    >>> with Session(token=mytoken) as sess:
    ...     sess.get("https://science.example.com/api/important/data")

    To explicitly *require* that a token is discovered, and *disable*
    any X.509 discovery:

    >>> with Session(token=True, x509=False) as sess:
    ...     sess.get("https://science.example.com/api/important/data")

    To use default authorisation discovery, but fail if no credentials
    are discovered:

    >>> with Session(fail_if_noauth=True) as sess:
    ...     sess.get("https://science.example.com/api/important/data")

    To disable all authorisation discovery:

    >>> with Session(force_noauth=True) as sess:
    ...     sess.get("https://science.example.com/api/important/data")
    """
    __attrs__ = requests.Session.__attrs__ = [
        "token",
    ]


# update the docstrings to include the same parameter info
for _obj in (Session, SessionAuthMixin):
    _obj.__doc__ = _obj.__doc__.format(parameters=_auth_session_parameters)


def request(method, url, *args, session=None, **kwargs):
    """Send a request of the specific method to the specified URL.

    Parameters
    ----------
    method : `str`
        The method to use.

    url : `str`,
        The URL to request.

    session : `requests.Session`, optional
        The connection session to use, if not given one will be
        created on-the-fly.

    args, kwargs
        All other keyword arguments are passed directly to
        `requests.Session.request`

    Returns
    -------
    resp : `requests.Response`
        the response object

    See also
    --------
    requests.Session.request
        for information on how the request is performed
    """
    # user's session
    if session:
        return session.request(method, url, *args, **kwargs)

    # new session
    sess_kwargs = {k: kwargs.pop(k) for k in (
        "cert",
        "token",
        "token_audience",
        "token_scope",
    ) if k in kwargs}
    with Session(url=url, **sess_kwargs) as session:
        return session.request(method, url, *args, **kwargs)


_request_wrapper_doc = """
    Send an HTTP {METHOD} request to the specified URL with IGWN Auth attached.

    Parameters
    ----------
    url : `str`
        The URL to request.

    session : `requests.Session`, optional
        The connection session to use, if not given one will be
        created on-the-fly.

    args, kwargs
        All other keyword arguments are passed directly to
        :meth:`requests.Session.{method}`

    Returns
    -------
    resp : `requests.Response`
        the response object

    See also
    --------
    requests.Session.{method}
        for information on how the request is performed
""".strip()


def _request_wrapper_factory(method):
    def _request_wrapper(url, *args, session=None, **kwargs):
        return request(method, url, *args, session=session, **kwargs)

    _request_wrapper.__doc__ = _request_wrapper_doc.format(
        method=method,
        METHOD=method.upper(),
    )
    return _request_wrapper


# request methods
delete = _request_wrapper_factory("delete")
get = _request_wrapper_factory("get")
head = _request_wrapper_factory("head")
patch = _request_wrapper_factory("patch")
post = _request_wrapper_factory("post")
put = _request_wrapper_factory("put")
