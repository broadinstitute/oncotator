from abc import abstractmethod


class SegmentDatasource():
    """Abstract class for datasources that can annotate segments (or regions).

    """

    @abstractmethod
    def annotate_segment(self, m):
        """Annotate the mutation as if it was a segment.

        :returns MutationData
        """
        # The default implementation raise a NotImplementedError
        # to ensure that any subclasses must override this method.
        raise NotImplementedError("Abstract class " + __name__ + " instantiated.")
