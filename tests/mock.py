"""Utility to mock function calls."""

import inspect

history = []
returns = {}


def mock(monkeypatch, module, func_name, returns=None):
    """Append the function call to the history.

    Save all of the arguments of each function call in the order they
    were called.
    """
    func = module.__dict__.get(func_name)
    sig = inspect.signature(func)
    arg_names = [p for p in sig.parameters]

    def mocked(*args, **kwargs):
        hist = {arg_names[i]: a for i, a in enumerate(args)}
        hist['module'] = module.__name__
        hist['func'] = func_name
        history.append(hist)
        if returns:
            return returns

    monkeypatch.setattr(module, func_name, mocked)


def filter(module, func):
    """Filter the history to get only specific function calls."""
    hist = []
    filt = [h for h in history if h['module'] == module and h['func'] == func]
    for h in filt:
        call = {k: v for k, v in h.items() if k not in ['module', 'func']}
        hist.append(call)
    return hist
