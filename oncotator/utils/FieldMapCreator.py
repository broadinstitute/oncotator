from oncotator.utils.ConfigUtils import ConfigUtils


class FieldMapCreator(object):
    """Generic class for creating a mapping of field name to an annotation name.  Under the hood this is a simple dict.
    """

    @staticmethod
    def choose_best_annotation(header, m, alternative_dict, deprioritize_input_annotations):
        """  Choose the best annotation in mut for the given header.
        :param header:
        :param m:
        :param alternative_dict:
        :param deprioritize_input_annotations:
        :return: if no suitable annotation, None
        """

        annotation_names = sorted(list(set(m.keys() + m.getAttributeNames())))
        hdr = header.lower()
        alternatives = alternative_dict.get(hdr, [])
        full_possibilities = [header] + alternatives
        result = None

        if not deprioritize_input_annotations:
            for p in full_possibilities:
                if p in annotation_names:
                    result = p
                    break
        else:
            second_choice = None
            is_choosing_second_choice = True
            for p in full_possibilities:
                if p in annotation_names and m.getAnnotation(p).getDatasource() != "INPUT":
                    result = p
                    is_choosing_second_choice = False
                    break
                if second_choice is None and p in annotation_names:
                    second_choice = p

            if is_choosing_second_choice:
                result = second_choice

        return result


    @staticmethod
    def create_field_map(headers, m, alternative_dict, is_render_internal_fields=True,
                            exposed_fields=None, prepend="i_", deprioritize_input_annotations=False):
        """
        Create a mapping for output header to the best input annotation.

        This can handle prepend fields (attach the prepend to internal fields), exposed fields (ones not in the list of headers, but should not have a prepend),

        :param is_render_internal_fields: Whether annotations not assigned to headers (or superseded by other annotations) should be included in this map.
        :type is_render_internal_fields bool
        :param exposed_fields: list of fields that, if found, should never receive a prepend.
        :param prepend: The prepend to put on internal fields, if any are detected.  If is_render_internal_fields is False, this parameter does nothing.
        :param deprioritize_input_annotations: If an annotation with the exact name of the header is found AND it has a datasource of "INPUT",
            use one the annotations instead.  This is useful in cases where we want to reannotate.  This effectively handles aliases.
        :param headers:  List of headers that need to be populated for rendering.  For example, the columns in a TCGA MAF
        :param m: MutationData to scrape available annotations
        :param alternative_dict: Dictionary of header to list of annotations.  Usually, created from a config file.
        :return: dict of header:annotation name (one annotation name) that should be used for this output rendering.
        """
        result = dict()
        if prepend is None:
            prepend = ""

        if exposed_fields is None:
            exposed_fields = set()

        # Process each header and find the first alias.  If an annotation exists with the exact same name as the header
        #   use that unless deprioritization is in effect.
        annotation_names = sorted(list(set(m.keys() + m.getAttributeNames())))

        for h in headers:
            choice = FieldMapCreator.choose_best_annotation(h, m, alternative_dict, deprioritize_input_annotations)
            if choice is None:
                choice = h
            result[h] = choice

        # Now populate internal fields, if requested.
        if is_render_internal_fields:

            annotation_names_used = result.values()
            internal_field_dict = dict()
            sAnnotations = set(annotation_names)
            internal_fields = sAnnotations.difference(annotation_names_used)

            # Create a dict to do a lookup of annotation to the column to use.
            reverseAlternativeDict = ConfigUtils.buildReverseAlternativeDictionary(alternative_dict)

            for i in internal_fields:
                if i.startswith('_') or i == "transcripts":
                    continue

                no_prepend_name = i
                if prepend != "" and i.startswith(prepend):
                    no_prepend_name = i.replace(prepend, "")

                field_alt_dict = {i: [prepend+i, no_prepend_name]}
                choice = FieldMapCreator.choose_best_annotation(i, m, field_alt_dict, deprioritize_input_annotations)
                key_to_use = reverseAlternativeDict.get(i,i)
                if prepend.strip() == "" or i.startswith(prepend) or i in exposed_fields:
                    internal_field_dict[key_to_use] = choice
                else:
                    internal_field_dict[prepend + key_to_use] = choice

            result.update(internal_field_dict)

        return result
