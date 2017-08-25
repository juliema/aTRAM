"""Utility to mock function calls."""

import inspect
from itertools import cycle
from contextlib import contextmanager

history = []
monkeypatch = None


def it(module, func_name, returns=None):
    """
    Append the function call to the history.

    Save all of the arguments of each function call in the order they
    were called. You can pass in a set of return values in the returns
    argument. If it is a list the return values will cycle thru the list
    returning each item in turn. For simple data types the cycle is one item.
    """
    global history

    func = module.__dict__.get(func_name)
    sig = inspect.signature(func)
    arg_names = [p for p in sig.parameters]

    if returns is not None:
        returns = returns if isinstance(returns, list) else [returns]
        returns = cycle(returns)

    def mocked(*args, **kwargs):
        hist = {arg_names[i]: a for i, a in enumerate(args)}
        hist['module'] = module.__name__
        hist['func'] = func_name
        history.append(hist)
        if returns:
            return next(returns)

    monkeypatch.setattr(module, func_name, mocked)


def context(module, func_name, returns):
    """
    Append the function call to the history.

    This version or "mock.it" is used when you're mocking a function call
    inside a with statement.
    """
    global history

    func = module.__dict__.get(func_name)
    sig = inspect.signature(func)
    arg_names = [p for p in sig.parameters]

    returns = returns if isinstance(returns, list) else [returns]
    returns = cycle(returns)

    @contextmanager
    def mocked(*args, **kwargs):
        hist = {arg_names[i]: a for i, a in enumerate(args)}
        hist['module'] = module.__name__
        hist['func'] = func_name
        history.append(hist)
        yield next(returns)

    monkeypatch.setattr(module, func_name, mocked)


def filter(module, func=None):
    """
    Filter the history to get only specific function calls.

    Remove the given module and (optionally) functions names.
    """
    hist = []

    filtered = [h for h in history if h['module'] == module]
    if func:
        filtered = [h for h in filtered if h['func'] == func]

    for h in filtered:
        func_call = {k: v for k, v in h.items() if k != 'module'}
        if func:
            del func_call['func']
        hist.append(func_call)

    return hist
