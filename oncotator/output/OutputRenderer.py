# LICENSE_GOES_HERE


from abc import abstractmethod


class OutputRenderer(object):
    """
    This is the base class for rendering MutationData instances.
    """

    @abstractmethod
    def __init__(self, filename, configFile="", otherOptions=None):
        """
        Constructor
        """
        raise NotImplementedError
    
    @abstractmethod
    def renderMutations(self, mutations, metadata=None, comments=None):

        """
        Render the given mutations.  Subclasses should also handle a list of comments.  However, the rendering of
            comments may be optional.

        :param mutations:
        :param metadata:
        :param comments:
        :raise: NotImplementedError
        """
        raise NotImplementedError