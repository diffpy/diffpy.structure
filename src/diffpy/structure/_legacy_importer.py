#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2017 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""
Support import of old camel-case module names with DeprecationWarning.

The imported camel-case modules are aliases for the current module
instances.  Their `__name__` attributes are thus all in lower-case.

This module is deprecated and will be removed in the future.

NOTE: this module must be only imported from `diffpy.Structure`.
"""


import sys
from warnings import warn
import six

if six.PY2:
    import importlib
    class mock_importlib_abc(object):
        MetaPathFinder = object
        Loader = object
    importlib.abc = mock_importlib_abc
    sys.modules.setdefault('importlib.abc', mock_importlib_abc)
    del mock_importlib_abc

import importlib.abc

WMSG = "Module {!r} is deprecated.  Use {!r} instead."

# ----------------------------------------------------------------------------

class FindRenamedStructureModule(importlib.abc.MetaPathFinder):

    prefix = 'diffpy.Structure.'

    def find_spec(self, fullname, path=None, target=None):
        # only handle submodules of diffpy.Structure
        if not fullname.startswith(self.prefix):
            return None
        lcname = fullname.lower()
        spec = importlib.util.find_spec(lcname)
        if spec is not None:
            spec.name = fullname
            spec.loader = MapRenamedStructureModule()
        return spec


    if six.PY2:
        def find_module(self, fullname, path):
            # only handle submodules of diffpy.Structure
            loader = None
            if fullname.startswith(self.prefix):
                loader = MapRenamedStructureModule()
            return loader

# end of class FindRenamedStructureModule

# ----------------------------------------------------------------------------

class MapRenamedStructureModule(importlib.abc.Loader):
    """
    Loader for old camel-case module names.

    Import the current module and alias it under the old name.
    """

    def create_module(self, spec):
        lcname = spec.name.lower()
        mod = importlib.import_module(lcname)
        sys.modules[spec.name] = mod
        warn(WMSG.format(spec.name, lcname), DeprecationWarning, stacklevel=2)
        return mod


    def exec_module(self, module):
        return


    if six.PY2:
        from collections import namedtuple
        ModuleSpec = namedtuple('ModuleSpec', 'name')

        def load_module(self, fullname):
            spec = self.ModuleSpec(fullname)
            return self.create_module(spec)

# end of class MapRenamedStructureModule

# ----------------------------------------------------------------------------

# show deprecation warning for diffpy.Structure
warn(WMSG.format('diffpy.Structure', 'diffpy.structure'),
     DeprecationWarning, stacklevel=2)

# install meta path finder for diffpy.Structure submodules
sys.meta_path.append(FindRenamedStructureModule())

# End of file
