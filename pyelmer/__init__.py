"""A python interface to Elmer."""
import pyelmer.elmer
import pyelmer.execute
import pyelmer.gmsh_utils
import pyelmer.post


from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
