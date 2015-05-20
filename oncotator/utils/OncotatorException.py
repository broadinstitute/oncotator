class OncotatorException(RuntimeError):
    """
    Generic runtime exception for ReCapSeg (incl. AllelicCapSeg) and auxiliary functionality.
    """

    def __init__(self, value, cause_exception_err=None, cause_exception_traceback=None):
        """Given an exception as well."""
        self.value = value

        if cause_exception_err is None:
            self.root_err_msg = ""
        else:
            self.root_err_msg = str(cause_exception_err)

        if cause_exception_traceback is None:
            self.root_traceback = ""
        else:
            self.root_traceback = str(cause_exception_traceback)

    def __str__(self):
        if self.root_err_msg is None:
            return repr(self.value)
        else:
            return str(self.value + "\nOriginal Exception:\n" + self.root_err_msg + "\n" +  self.root_traceback)