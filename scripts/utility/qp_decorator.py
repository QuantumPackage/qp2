from functools import wraps

def cache(func):
    """
    A decorator for lazy evaluation off true function
    """
    saved = {}

    @wraps(func)
    def newfunc(*args):
        if args in saved:
            return saved[args]

        result = func(*args)
        saved[args] = result
        return result
    return newfunc
