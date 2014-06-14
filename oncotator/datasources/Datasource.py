# LICENSE_GOES_HERE

import logging
import abc

try:
    import pysam
except ImportError:
    warningString = (
        'ERROR: Could not load pysam.  Some features will be disabled (e.g. COSMIC annotations) and may cause Oncotator to fail.')
    print(warningString)
    logging.getLogger('').error(warningString)


class Datasource(object):
    """
    An individual datasource used for annotation.  This serves as a base class for attributes
    and methods that are common to all types of datasources.  Subclasses of Datasource will
    define behavior for more specific data types.

    src_file
        Absolute path to the annotation datasource file.

    title
        Title string to be appended to each datasource header.

    version
        Version of datasource. Will be displayed in comments of output file.


    Note for developers:  When creating new types of datasources, you may need to be aware of the sorting position.
    For example, datasources that rely on the gene annotation must be placed after the GAF (or other transcript
        datasource).  See DatasourceCreator::sortDatasources()
    """

    def __init__(self, src_file, title='', version=None, default_missing_value=''):
        """
        This method is used to load/parse the annotation datasource.  Should be overridden
        for custom parsing/loading

        The title option should be used to label output headers. Please see the Generic_Gene_DataSource
        and Generic_GenomicPosition_DataSource Datasource subclasses as examples.

        If added annotation fields are to be outputted by AnnotationSet.write_data method,
        then a output_headers list must be instantiated.  Please see the Generic_Gene_DataSource
        and Generic_GenomicPosition_DataSource Datasource subclasses as examples.

        TODO: Update this documentation

        """
        self.src_file = src_file
        self.title = title
        self.version = str(version)
        self.output_headers = []
        self.hashcode = ""

    def attach_to_class(self, cls, name, options):
        """
        Do not touch.  This method is necessary for the AnnotationSetMeta class to work.

        This coded lets AnnotationSet become aware of their constituent Datasource classes.
        """
        self.cls = cls
        self.name = name
        self.options = options
        if self.title is None:
            self.title = name
        options.add_datasource(self)

    @abc.abstractmethod
    def annotate_mutation(self, mutation):
        # The default implementation raise a NotImplementedError
        # to ensure that any subclasses must override this method.
        raise NotImplementedError

    def set_hashcode(self, value):
        self.hashcode = value

    def get_hashcode(self):
        return self.hashcode



