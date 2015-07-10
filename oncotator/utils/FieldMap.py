from oncotator.utils.ConfigUtils import ConfigUtils


class FieldMap(object):
    """Generic class for storing a mapping of field name to an annotation name.  Under the hood this is a simple dict.
    """


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
            chosen = h
            is_second_choice = h in m.keys() and \
                               (deprioritize_input_annotations and m.getAnnotation(h).getDatasource() == "INPUT")
            if h in annotation_names and not is_second_choice:
                # Map the header to the header, but that has already been done
                pass
            else:
                if h in annotation_names:
                    second_choice = h
                else:
                    second_choice = None

                # Search for the first alias and use that.  If it is an INPUT and we care, keep looking.
                hdr = h.lower()
                if hdr in alternative_dict.keys():
                    alternatives = alternative_dict[hdr]

                    if not deprioritize_input_annotations:
                        for alternative in alternatives:
                            if alternative in annotation_names:
                                chosen = alternative
                                break
                    else:
                        # We want an alias that is present and (if requested by caller) not an input annotation
                        is_take_second_choice = True
                        for alternative in alternatives:
                            if alternative in annotation_names and m.getAnnotation(h).getDatasource() != "INPUT":
                                chosen = alternative
                                is_take_second_choice = False
                                break
                            elif alternative in annotation_names:
                                if second_choice is None:
                                    second_choice = alternative

                        if is_take_second_choice and second_choice is not None:
                            chosen = second_choice
            result[h] = chosen

        # TODO: What about internal fields?  What about internal fields with a prepend that needs to be comapred against something without a prepend?
        # Now populate internal fields, if requested.
        if is_render_internal_fields:
            annotation_names_used = result.values()

            internal_field_dict = dict()

            # Create a dict to do a lookup of annotation to the column to use.
            reverseAlternativeDict = ConfigUtils.buildReverseAlternativeDictionary(alternative_dict)

            sAnnotations = set(annotation_names)
            internalFields = sAnnotations.difference(annotation_names_used)
            for i in internalFields:
                if not i.startswith('_') and i is not "transcripts":
                    key_to_use = reverseAlternativeDict.get(i,i)
                    if prepend.strip() == "" or i.startswith(prepend) or i in exposed_fields:
                        internal_field_dict[key_to_use] = i
                    else:
                        internal_field_dict[prepend + key_to_use] = i
            result.update(internal_field_dict)


        return result
