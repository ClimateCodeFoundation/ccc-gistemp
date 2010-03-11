"""This module extends sys.path to include this module's parent directory.

This module is imported for its deliberate side effect of extending
``sys.path`` to include the parent directory of this module. The idea
is that any script in the ``tool`` directory can import code in the
``code`` package, for example, as follows::

    import extend_path
    from code import step5

"""
__docformat__ = "restructuredtext"

import sys
import os

my_path = os.path.abspath(__file__)
parent = os.path.dirname(os.path.dirname(my_path))
if parent not in sys.path:
    sys.path[0:0] = [parent]

del my_path, parent
