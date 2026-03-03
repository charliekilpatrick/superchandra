"""
Input validation helpers for the superchandra pipeline.

Provides numeric string checks and similar validators used by the CLI
and pipeline (e.g. validating ra/dec arguments).
"""


def is_number(num):
    """
    Check whether a value is numeric (convertible to float).

    Parameters
    ----------
    num : str or any
        Value to test (typically a string from the command line).

    Returns
    -------
    bool
        True if the value can be converted to float, False otherwise
        (e.g. None, empty string, or non-numeric string).
    """
    try:
        float(num)
    except (ValueError, TypeError):
        return False
    return True
