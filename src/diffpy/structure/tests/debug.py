#!/usr/bin/env python3
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################


"""\
Convenience module for debugging the unit tests using

python -m diffpy.structure.tests.debug

Exceptions raised by failed tests or other errors are not caught.
"""


if __name__ == '__main__':
    from diffpy.structure.tests import testsuite
    suite = testsuite()
    suite.debug()
