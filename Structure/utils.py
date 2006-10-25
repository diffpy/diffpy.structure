"""Small shared functions.
"""

__id__ = "$Id$"

def isfloat(s):
    """True if argument can be converted to float"""
    try:
        x = float(s)
        return True
    except ValueError:
        pass
    return False
