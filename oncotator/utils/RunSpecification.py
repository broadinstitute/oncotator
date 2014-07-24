
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

    def _init_(self):
        self._inputCreator = None
        self._outputRenderer = None
        self._inputFilename = None
        self._outputFilename = None
        self._manualAnnotations = None
        self._defaultAnnotations = None
        self._datasources = None
        self._isMulticore = False
        self._numCores = None
        self._cache_url = None
        self._isReadOnlyCache=True
        self._isSkipNoAlts = False
        pass

    def get_cache_url(self):
        return self._cache_url

    def set_cache_url(self, value):
        self._cache_url = value

    def del_cache_url(self):
        del self._cache_url

    def set_is_read_only_cache(self, value):
        self._isReadOnlyCache = value

    def get_is_read_only_cache(self):
        return self._isReadOnlyCache

    def del_is_read_only_cache(self):
        del self._isReadOnlyCache

    def get_is_multicore(self):
        return self._isMulticore


    def get_num_cores(self):
        return self._numCores


    def set_is_multicore(self, value):
        self._isMulticore = value


    def set_num_cores(self, value):
        self._numCores = value


    def del_is_multicore(self):
        del self._isMulticore


    def del_num_cores(self):
        del self._numCores


    def get_datasources(self):
        return self._datasources


    def set_datasources(self, value):
        self._datasources = value


    def del_datasources(self):
        del self._datasources


    def get_input_creator(self):
        return self._inputCreator


    def get_output_renderer(self):
        return self._outputRenderer


    def get_manual_annotations(self):
        return self._manualAnnotations


    def set_input_creator(self, value):
        self._inputCreator = value


    def set_output_renderer(self, value):
        self._outputRenderer = value

    def set_manual_annotations(self, value):
        self._manualAnnotations = value

    def del_input_creator(self):
        del self._inputCreator

    def del_output_renderer(self):
        del self._outputRenderer

    def del_manual_annotations(self):
        del self._manualAnnotations

    def get_default_annotations(self):
        return self._defaultAnnotations

    def set_default_annotations(self, value):
        self._defaultAnnotations = value

    def del_default_annotations(self):
        del self._defaultAnnotations

    def get_is_skip_no_alts(self):
        return self._isSkipNoAlts

    def set_is_skip_no_alts(self, value):
        self._isSkipNoAlts = value

    def del_is_skip_no_alts(self):
        del self._isSkipNoAlts

    def get_annotating_type(self):
        return self._annotating_type

    def set_annotating_type(self, value):
        self._annotating_type = value

    def del_annotating_type(self):
        del self._annotating_type

    def initialize(self, inputCreator, outputRenderer, manualAnnotations=None, datasources=None, isMulticore=False, numCores=4, defaultAnnotations=None, cacheUrl=None, read_only_cache=True, is_skip_no_alts=False, annotating_type=None):
        self._inputCreator = inputCreator
        self._outputRenderer = outputRenderer
        self._manualAnnotations = manualAnnotations if manualAnnotations is not None else dict()
        self._datasources = datasources if datasources is not None else []
        self._isMulticore = isMulticore
        self._numCores = numCores
        self._defaultAnnotations = defaultAnnotations if defaultAnnotations is not None else dict()
        self._cache_url = cacheUrl
        self._isReadOnlyCache = read_only_cache
        self._isSkipNoAlts = is_skip_no_alts
        self._annotating_type = annotating_type if annotating_type is not None else RunSpecification.ANNOTATE_MUTATIONS

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
