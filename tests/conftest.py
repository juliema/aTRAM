"""A module that sets up the tests."""

import pytest
import tests.mock as mock


@pytest.fixture(scope='function', autouse=True)
def mocker(monkeypatch):
    """Testing."""
    mock.history = []
    mock.monkeypatch = monkeypatch
