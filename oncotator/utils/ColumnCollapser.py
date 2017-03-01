from abc import abstractmethod


class ColumnCollapser(object):
    """
    Interface for column collapsing strategies.
    """
    def __init__(self, config_file=None):
        pass

    @abstractmethod
    def update_mutation(self, mut, new_annotation_source=None, copy_old_suffix=None):
        """Given a mutation, update the relevant annotations with new values.  Please note that this updates in place.

        :param copy_old_suffix: Create another annotation with the old value in it.  Use this suffix  If None, do not create the annotation.
        :param new_annotation_source: string giving the annotation source that should be used for this updating of the mutation
            If set to None, leave the old annotation source in place (use sparingly).
        :param mut: MutationData instance
        :return: None
        """
        raise NotImplementedError(__name__ + " cannot be instantiated so this method cannot be called..")

    @abstractmethod
    def retrieve_new_annotations_added(self, mut, copy_old_suffix=None):
        """
        This method simply returns a list of annotations that would be created if run through update_mutation.

        :param mut:annotated mutation with numeric fields to collapse (presumably)
        :type mut MutationData
        :param copy_old_suffix: suffix for the new annotations (that will hold the old values).  Please note that
        this method will return an empty list if None (default) is specified
        :return: list of str with annotation names.
        """
        raise NotImplementedError(__name__ + " cannot be instantiated so this method cannot be called.")
