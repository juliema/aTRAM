"""Utility to mock function calls."""

import inspect

HISTORY = []


def mock(monkeypatch, module, func_name):
    """Append the function call to the history.

    Save all of the arguments of each function call in the order they
    were called.
    """
    func = module.__dict__.get(func_name)
    sig = inspect.signature(func)
    arg_names = [p for p in sig.parameters]

    def mocked(*args, **kwargs):
        history = {arg_names[i]: a for i, a in enumerate(args)}
        history['module'] = module.__name__
        history['func'] = func_name
        HISTORY.append(history)

    monkeypatch.setattr(module, func_name, mocked)


def filter(module, func):
    """Filter the history to get only specific function calls."""
    return [h for h in HISTORY if h['module'] == module and h['func'] == func]
