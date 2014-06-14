# LICENSE_GOES_HERE
from oncotator.Annotation import Annotation
from oncotator.cache.CacheManager import CacheManager
from oncotator.utils.Hasher import Hasher


"""
Created on Nov 2, 2012

@author: lichtens
"""
import logging
from utils.version import VERSION


class Annotator(object):
    """
    The Annotator is the entry point to actually perform the annotating of mutations.  The Annotator contains one input creator (IC), one output creator (OC), and a list of datasources.   
    
    This class is responsible for the coordination of the annotating process, not the annotations themselves (this is handled by the datasources).
    
    For information on how to initialize the input and output creator, please see the documentation of those classes.

    See the RunSpecification class, which allows for more control
        of an annotator.

    Example usage (with RunSpec and no multicore usage):
    # Create a run configuration to pass to the Annotator class.  See OncotatorCLIUtils.getSupportedOutputFormats()
    #   and OncotatorCLIUtils.getSupportedInputFormats() for allowed inputFormat and outputFormat strings.
    manualOverrides = {'fake_annotation':'picard', 'fake_annotation2':'worf'}
    runConfig = OncotatorCLIUtils.createRunConfig(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations=manualOverrides, datasourceDir="/home/onco/dbs", isMulticore=False)

    annotator = Annotator()
    annotator.initialize(runConfig)
    annotator.annotate()

    Example usage (used in testing, without RunSpec):
        # Assumed myIC and myOC have been initialized as the proper Input and Output Creators, respectively.
        # 1) Initialize the Annotator
        annotator = Annotator()
        annotator.setInputCreator(myIC)
        annotator.setOutputCreator(myOC)
        # 1a) For each datasource (instance of a datasource class), add it to the annotator.
        for datasource in myDataSources:
            annotator.addDatasource(datasource)
        # 2)  Produce the output
        filePointer = annotator.annotate()
    
    NOTE:  While multicore information is passed into the Annotator, currently, nothing is implemented that uses multicore.    
    """

    def __init__(self):
        """
        options should contain the following name-value pairs as a dict: 
        
        Create a new instance of Annotator.
        
        In order to specify the input and output creators and datasources, use the set and addDatasource methods.
        
        """
        self._inputCreator = None
        self._outputRenderer = None
        self._datasources = []
        self.logger = logging.getLogger(__name__)
        self._manualAnnotations = dict()
        self._defaultAnnotations = dict()
        self._isMulticore = None
        self._numCores = None
        self._cacheManager = CacheManager()
        self._cacheManager.initialize(None, "not_used")
        self._cache_stats = {"miss": 0, "hit":0}
        self._is_skip_no_alts = False
        pass

    def getIsMulticore(self):
        return self.__isMulticore

    def getNumCores(self):
        return self.__numCores

    def setIsMulticore(self, value):
        self.__isMulticore = value

    def setNumCores(self, value):
        self.__numCores = value

    def setInputCreator(self, inputCreator):
        self._inputCreator = inputCreator
        
    def setOutputRenderer(self, outputCreator):
        self._outputRenderer = outputCreator
    
    def setManualAnnotations(self, value):
        self._manualAnnotations = value

    def setDefaultAnnotations(self, value):
        self._defaultAnnotations = value

    def create_db_dir_key(self):
        """Create the db_dir_key for this annotation configuration.  Requires the datasources."""
        self.logger.info("Generating db-dir key from datasources...")
        hasher = Hasher()
        for ds in self._datasources:
            self.logger.info(ds.title + " " + ds.version + " md5: " + ds.get_hashcode())
            hasher.update(ds.get_hashcode())
        db_dir_key = Hasher.md5_hash(hasher.hexdigest())
        self.logger.info("Final db-dir md5: " + db_dir_key)
        return db_dir_key

    def create_db_dir_key_simple(self):
        """Create the db_dir_key for this annotation configuration.  Requires the datasources."""
        db_dir_key = Hasher.md5_hash(self.createHeaderString(False))
        return db_dir_key

    def initialize_cache_manager(self, runSpec):
        """Do not bother calculating the db_dir_key if the cache is not being used. """
        cache_url = runSpec.get_cache_url()
        if cache_url is not None and cache_url != "":
            db_dir_key = self.create_db_dir_key()
        else:
            db_dir_key = "never_used"
        self._cacheManager = CacheManager()
        self._cacheManager.initialize(cache_url, db_dir_key, is_read_only=runSpec.get_is_read_only_cache())

    def initialize(self,runSpec):
        """ Given a RunSpecification instance, initialize self properly.  Do not start annotation.
        """
        self.setInputCreator(runSpec.inputCreator)
        self.setOutputRenderer(runSpec.outputRenderer)
        self.setManualAnnotations(runSpec.manualAnnotations)
        self.setDefaultAnnotations(runSpec.defaultAnnotations)
        self._datasources = runSpec.datasources
        self.setIsMulticore(runSpec.get_is_multicore())
        self.setNumCores(runSpec.get_num_cores())
        self._cache_stats = {"miss": 0, "hit":0}
        self._is_skip_no_alts = runSpec.get_is_skip_no_alts()
        self.initialize_cache_manager(runSpec)

    def addDatasource(self, datasource):
        self._datasources.append(datasource)

    def _createMetadata(self):
        metadata = self._inputCreator.getMetadata()
        metadata.update(self._createManualAnnotationsForMetadata(self._manualAnnotations))
        return metadata

    def _createComments(self):
        comments = self._inputCreator.getComments()
        comments.append(self.createHeaderString())
        return comments

    def annotate_mutations(self, mutations):
        mutations = self._annotate_mutations_using_datasources(mutations)
        if mutations is None:
            self.logger.warn("Mutation list points to None after annotation.")

        mutations = self._applyDefaultAnnotations(mutations, self._defaultAnnotations)
        if mutations is None:
            self.logger.warn("Mutation list points to None after default annotations.")

        mutations = self._applyManualAnnotations(mutations, self._manualAnnotations)
        if mutations is None:
            self.logger.warn("Mutation list points to None after manual annotations.")

        return mutations

    def annotate(self):
        """
        Annotate the given mutations specified in the input.

        Call this after the input, output, and datasources have been set.

        :return: outputFilename
        """
        self.logger.info("Annotating with " + str(len(self._datasources)) + " datasources: " + self.createHeaderString())
        
        mutations = self._inputCreator.createMutations()
        if mutations is None: 
            self.logger.warn("Mutation list points to None after creation.")

        mutations = self.annotate_mutations(mutations)

        comments = self._createComments()
        metadata = self._createMetadata()

        filename = self._outputRenderer.renderMutations(mutations, metadata=metadata, comments=comments)

        self.logger.info("Closing cache: (misses: " + str(self._cache_stats['miss']) + "  hits: " + str(self._cache_stats['hit']) + ")")
        self._cacheManager.close_cache()

        return filename
    
    def _applyManualAnnotations(self, mutations, manualAnnotations):
        manualAnnotationKeys = manualAnnotations.keys()
        for m in mutations:
            for k in manualAnnotationKeys:
                # newRequired = False allows this call to overwrite the previous value.
                m.createAnnotation(k, manualAnnotations[k], annotationSource="MANUAL", newRequired=False)
            yield m

    def _applyDefaultAnnotations(self, mutations, defaultAnnotations):
        defaultAnnotationsKeys = defaultAnnotations.keys()
        for m in mutations:
            mKeys = m.keys()
            for k in defaultAnnotationsKeys:
                if k not in mKeys:
                    m.createAnnotation(k, defaultAnnotations[k], annotationSource="DEFAULT")
                if m[k] == "":
                    m.getAnnotation(k).setDatasource("DEFAULT")
                    m.getAnnotation(k).setValue(defaultAnnotations[k])
            yield m

    def _createManualAnnotationsForMetadata(self, manualAnnotations):
        result = {}
        manualAnnotationKeys = manualAnnotations.keys()
        for k in manualAnnotationKeys:
            result[k] = Annotation(manualAnnotations[k], datasourceName="MANUAL")
        return result

    def createHeaderString(self, is_giving_oncotator_version=True):
        """
        Create a default header string that lists version of Oncotator and datasource information.

        :return: string
        """
        onco_string = ""
        if is_giving_oncotator_version:
            onco_string = "Oncotator " +  VERSION + " |"

        datasourceStrings = []
        for ds in self._datasources:
            datasourceStrings.append(" " + ds.title + " " + ds.version + " ")
        
        return onco_string + "|".join(datasourceStrings)
    
    def _annotate_mutations_using_datasources(self, mutations):
        if len(self._datasources) == 0:
            self.logger.warn("THERE ARE NO DATASOURCES REGISTERED")
        for m in mutations:

            # If the alt_allele_seen annotation is present and False, skip this mutation
            if self._is_skip_no_alts and m.get("alt_allele_seen", "True") == "False":
                continue

            annot_dict = self._cacheManager.retrieve_cached_annotations(m)

            if annot_dict is None:
                self._cache_stats['miss'] += 1
                for datasource in self._datasources:
                    m = datasource.annotate_mutation(m)
                self._cacheManager.store_annotations_in_cache(m)
            else:
                self._cache_stats['hit'] += 1
                m.addAnnotations(annot_dict)
            yield m