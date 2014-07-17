
class RunSpecification(object):
    """ This class contains a specification for an oncotator run.

    This class is quite simple and simply holds parameters.  Annotator classes will be expected to be able to initialize and annotate
        from an instance of this class.

        inputCreator -- instance
        outputCreator -- instance
        manualAnnotations -- A dict of name-value pairs that will be annotated onto all mutations.
            name is the name of the annotation.
            value is the value of the annotation.
            The source will always be listed as MANUAL
            TODO: The Annotator takes care of some of this.  Move the relevant documentation.
        datasources -- A list of datasources (instance of Datasource).
        isMulticore -- use multicore processing, where available (True/False)
        numCores -- number of cores to use if isMulticore is True.  Otherwise, this is ignored.
        cache_url -- if None, implies that there is no cache.
        """

    ANNOTATE_SEGMENTS = "segments"
    ANNOTATE_MUTATIONS = "mutations"

    def __init__(self):
        self.__inputCreator = None
        self.__outputRenderer = None
        self.__inputFilename = None
        self.__outputFilename = None
        self.__manualAnnotations = None
        self.__defaultAnnotations = None
        self.__datasources = None
        self.__isMulticore = False
        self.__numCores = None
        self.__cache_url = None
        self.__is_read_only_cache=True
        self.__is_skip_no_alts = False
        pass

    def get_cache_url(self):
        return self.__cache_url

    def set_cache_url(self, value):
        self.__cache_url = value

    def del_cache_url(self):
        del self.__cache_url

    def set_is_read_only_cache(self, value):
        self.__is_read_only_cache = value

    def get_is_read_only_cache(self):
        return self.__is_read_only_cache

    def del_is_read_only_cache(self):
        del self.__is_read_only_cache

    def get_is_multicore(self):
        return self.__isMulticore


    def get_num_cores(self):
        return self.__numCores


    def set_is_multicore(self, value):
        self.__isMulticore = value


    def set_num_cores(self, value):
        self.__numCores = value


    def del_is_multicore(self):
        del self.__isMulticore


    def del_num_cores(self):
        del self.__numCores


    def get_datasources(self):
        return self.__datasources


    def set_datasources(self, value):
        self.__datasources = value


    def del_datasources(self):
        del self.__datasources


    def get_input_creator(self):
        return self.__inputCreator


    def get_output_renderer(self):
        return self.__outputRenderer


    def get_manual_annotations(self):
        return self.__manualAnnotations


    def set_input_creator(self, value):
        self.__inputCreator = value


    def set_output_renderer(self, value):
        self.__outputRenderer = value

    def set_manual_annotations(self, value):
        self.__manualAnnotations = value

    def del_input_creator(self):
        del self.__inputCreator

    def del_output_renderer(self):
        del self.__outputRenderer

    def del_manual_annotations(self):
        del self.__manualAnnotations

    def get_default_annotations(self):
        return self.__defaultAnnotations

    def set_default_annotations(self, value):
        self.__defaultAnnotations = value

    def del_default_annotations(self):
        del self.__defaultAnnotations

    def get_is_skip_no_alts(self):
        return self.__is_skip_no_alts

    def set_is_skip_no_alts(self, value):
        self.__is_skip_no_alts = value

    def del_is_skip_no_alts(self):
        del self.__is_skip_no_alts

    def get_annotating_type(self):
        return self.__annotating_type

    def set_annotating_type(self, value):
        self.__annotating_type = value

    def del_annotating_type(self):
        del self.__annotating_type

    def initialize(self, inputCreator, outputRenderer, manualAnnotations=None, datasources=None, isMulticore=False, numCores=4, defaultAnnotations=None, cacheUrl=None, read_only_cache=True, is_skip_no_alts=False, annotating_type=None):
        self.inputCreator = inputCreator
        self.outputRenderer = outputRenderer
        self.manualAnnotations = manualAnnotations if manualAnnotations is not None else dict()
        self.datasources = datasources if datasources is not None else []
        self.isMulticore = isMulticore
        self.numCores = numCores
        self.defaultAnnotations = defaultAnnotations if defaultAnnotations is not None else dict()
        self.cacheUrl = cacheUrl
        self.isReadOnlyCache = read_only_cache
        self.isSkipNoAlts = is_skip_no_alts
        self.annotating_type = annotating_type if annotating_type is not None else RunSpecification.ANNOTATE_MUTATIONS


    inputCreator = property(get_input_creator, set_input_creator, del_input_creator, "inputCreator's docstring")
    outputRenderer = property(get_output_renderer, set_output_renderer, del_output_renderer, "outputRenderer's docstring")
    manualAnnotations = property(get_manual_annotations, set_manual_annotations, del_manual_annotations, "manualAnnotations's docstring")
    defaultAnnotations = property(get_default_annotations, set_default_annotations, del_default_annotations, "Annotations that are populated only when the annotation does not exist or is empty ('' or None)")
    datasources = property(get_datasources, set_datasources, del_datasources, "datasources's docstring")
    isMulticore = property(get_is_multicore, set_is_multicore, del_is_multicore, "isMulticore's docstring")
    numCores = property(get_num_cores, set_num_cores, del_num_cores, "numCores's docstring")
    cacheUrl = property(get_cache_url, set_cache_url, del_cache_url, "cacheUrl's docstring")
    isReadOnlyCache = property(get_is_read_only_cache, set_is_read_only_cache, del_is_read_only_cache, "isReadOnlyCache's docstring")
    isSkipNoAlts = property(get_is_skip_no_alts, set_is_skip_no_alts, del_is_skip_no_alts, "isSkipNoAlts's docstring")
    annotating_type = property(get_annotating_type, set_annotating_type, del_annotating_type, "annotating type is static value indiciating what type of mutation to annotate.")
