# ngslib not supported on OS X so this global variable is used to  determine if library is installed at runtime
try:
    import ngslib
    NGSLIB_INSTALLED = True
except ImportError:
    NGSLIB_INSTALLED = False