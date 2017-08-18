"""Utility to mock function calls."""

import inspect
from itertools import cycle

history = []
monkeypatch = None


def it(module, func_name, returns=None):
    """Append the function call to the history.

    Save all of the arguments of each function call in the order they
    were called. You can pass in a set of return values in the returns
    function. If it is a list the return values will cycle thru the list
    returning each item in turn. For simple data types the cycle is one item.
    """
    func = module.__dict__.get(func_name)
    sig = inspect.signature(func)
    arg_names = [p for p in sig.parameters]

    if returns:
        returns = [returns] if not isinstance(returns, list) else returns
        returns = cycle(returns)

    def mocked(*args, **kwargs):
        hist = {arg_names[i]: a for i, a in enumerate(args)}
        hist['module'] = module.__name__
        hist['func'] = func_name
        history.append(hist)
        if returns:
            return next(returns)

    monkeypatch.setattr(module, func_name, mocked)


def filter(module, func):
    """Filter the history to get only specific function calls.

    Remove the given module and functions names. They are known.
    """
    hist = []
    filt = [h for h in history if h['module'] == module and h['func'] == func]

    for h in filt:
        func_call = {k: v for k, v in h.items() if k not in ['module', 'func']}
        hist.append(func_call)

    return hist
