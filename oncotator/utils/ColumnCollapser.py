import logging
from numpy.core.fromnumeric import mean
from oncotator.utils.ConfigUtils import ConfigUtils

class ColumnCollapser(object):
    """This class takes multiple values in an annotation (separated by "|") and collapses into a single value based on given rules.

        """

    MIN = "min"
    MEAN = "mean"

    # Any methods specified above should be added to this list.
    VALID_VALUES = [MIN, MEAN]

    RENDER_AS_INT_IF_POSSIBLE = [MIN]

    def __init__(self, config_file="column_collapse.config"):
        self._config_parser = ConfigUtils.createConfigParser(config_file, ignoreCase=False)
        self._columns_to_collapse = self._config_parser.options("columns_to_collapse")

        # Create a lookup table to get the method for each column
        self._method_dict = dict()
        if len(self._columns_to_collapse) != len(set(self._columns_to_collapse)):
            logging.getLogger(__name__).warn("Duplicate keys in " + config_file + " seen.  Some collapsing may produce unexpected values.")

        for c in self._columns_to_collapse:
            self._method_dict[c] = self._config_parser.get("columns_to_collapse", c).strip()

            # Default to mean if not specified
            if self._method_dict[c] is None or self._method_dict[c] == "":
                self._method_dict[c] = ColumnCollapser.MEAN

        # Basic validation
        problematic_method_assignments = dict()
        for c in self._columns_to_collapse:
            if self._method_dict[c] not in ColumnCollapser.VALID_VALUES:
                problematic_method_assignments[c] = self._method_dict[c]

        if len(problematic_method_assignments.keys()) > 0:
            raise ValueError("Invalid column collapsing specified: " + str(problematic_method_assignments))


    def _collapse_column(self, column_name, current_value):
        """Here is where values are collapsed based on the input values.  Assumed t be separated by a "|"
        current_value is assumed to be a string.
        """
        if column_name not in self._columns_to_collapse:
            return current_value

        list_vals = current_value.split("|")
        no_blank_vals = [v for v in list_vals if v.strip() != "" and v.strip() != "."]

        # if no "." is found in the value, assume that this should be rendered as an int, if MIN was chosen
        #   we assume it is an int if no "." is found.
        is_assuming_int = all([v.find(".") == -1 for v in no_blank_vals])

        try:
            if is_assuming_int:
                final_vals = [int(v) for v in no_blank_vals]
            else:
                final_vals = [float(v) for v in no_blank_vals]
        except ValueError as te:
            logging.getLogger(__name__).warn("Could not collapse " + column_name + ":" + current_value + " into one number.  Returning the input value.")
            return current_value

        if len(final_vals) == 0:
            return ""

        if self._method_dict[column_name] == ColumnCollapser.MIN:
            return str(min(final_vals))

        elif self._method_dict[column_name] == ColumnCollapser.MEAN:
            return str(mean(final_vals))

    def _collapse_columns(self, mut):
        """ Collapse the annotations in the mutation.

        :param mut:
        :return: dict with key value pairs of annotation_name:collapsed value.  This can then be used to update the mutation
        """
        result = dict()
        relevant_annotations = self._get_relevant_annotations(mut)
        for annot in relevant_annotations:
            result[annot] = self._collapse_column(annot, mut[annot])
        return result

    def update_mutation(self, mut, new_annotation_source=None, copy_old_suffix=None):
        """ Given a mutation, update the relevant annotations with new values.  Please note that this updates in place.

        :param copy_old_suffix: Create another annotation with the old value in it.  Use this suffix  If None, do not create the annotation.
        :param new_annotation_source: string giving the annotation source that should be used for this updating of the mutation
            If set to None, leave the old annotation source in place (use sparingly).
        :param mut:
        :return: None
        """
        update_dict = self._collapse_columns(mut)
        for u in update_dict.keys():
            annotation = mut.getAnnotation(u)
            if copy_old_suffix is not None:
                mut.createCopyAnnotation(annotation, u + copy_old_suffix)
            if new_annotation_source is not None:
                annotation.setDatasource(new_annotation_source)
            annotation.setValue(update_dict[u])

    def _get_relevant_annotations(self, mut):
        """

        :param mut:
        :return: list of annotations that will be collapsed given the mutation
        """

        return list(set(mut.keys()).intersection(self._columns_to_collapse))

    def retrieve_new_annotations_added(self, mut, copy_old_suffix=None):
        """
        This method simply returns a list of annotations that would be created if run through update_mutation.

        :param mut:annotated mutation with numeric fields to collapse (presumably)
        :type mut MutationData
        :param copy_old_suffix: suffix for the new annotations (that will hold the old values).  Please note that
        this method will return an empty list if None (default) is specified
        :return: list of str with annotation names.
        """
        if copy_old_suffix is None:
            return []
        relevant_annotations = self._get_relevant_annotations(mut)
        return [a + copy_old_suffix for a in relevant_annotations]
