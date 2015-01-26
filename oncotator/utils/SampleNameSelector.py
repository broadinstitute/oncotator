import logging

from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.MutUtils import MutUtils

from enum import Enum, unique

__author__ = 'louisb'



class SampleNameSelector(object):
    """
    Selects sample names from MutationData object.
    It reads a configuration file to determine which columns to gather the name from.
    It must be initialized with a representative mutation which will set the state of the selector to work on subsequent mutations.
    """
    @unique
    class Source(Enum):
        SAMPLE = 1
        TUMOR_NORMAL = 2
        TUMOR = 3
        NORMAL = 4


    def __init__(self, mut, configFile="sample_name_selection.config", section="SAMPLE_NAME"):
        config = ConfigUtils.createConfigParser(configFile)
        self.logger = logging.getLogger(__name__)
        aliases = ConfigUtils.buildAlternateKeyDictionaryFromConfig(config, section)
        self.configFile=configFile
        sampleAnnotation = self._getAnnotationFromAliases(mut, aliases["sample_name"])
        tumorAnnotation = self._getAnnotationFromAliases(mut, aliases["sample_tumor_name"])
        normalAnnotation = self._getAnnotationFromAliases(mut, aliases["sample_normal_name"])
        source_column = self._getSourceColumn(sampleAnnotation,tumorAnnotation,normalAnnotation)
        self._logSampleNameColumnDescription(source_column, sampleAnnotation, tumorAnnotation, normalAnnotation)
        self.sampleNameGrabber = self._getSampleNameGrabber(source_column, sampleAnnotation, tumorAnnotation, normalAnnotation)
        self.outputAnnotationName = self._deriveOutputAnnotationName(sampleAnnotation)
        self.annotationSource = self._deriveAnnotationSource(source_column)

    def getSampleName(self, mut):
        """
        Return the correct sample name from the MutationData
        :param mut: a MutationData
        :return: the correct sample name
        """
        try:
            return self.sampleNameGrabber(mut)
        except KeyError as e:
            self.logger.fatal("This MutationData is missing the expected sample name annotation.")
            raise e

    @staticmethod
    def _getAnnotationFromAliases(mut, aliases):
        if aliases is not None and mut is not None:
            for alias in aliases:
                if alias in mut:
                    return alias
        return None

    def _deriveAnnotationSource(self,source_column):
        sources = {
            self.Source.SAMPLE: "INPUT",
            self.Source.TUMOR_NORMAL: "OUTPUT",
            self.Source.TUMOR: "INPUT",
            self.Source.NORMAL: "INPUT",
            None: None
        }
        return sources[source_column]

    def _deriveOutputAnnotationName(self, sampleAnnotation):
        if sampleAnnotation is not None:
            return sampleAnnotation
        else:
            return MutUtils.SAMPLE_NAME_ANNOTATION_NAME

    def getOutputAnnotationName(self):
        return self.outputAnnotationName

    def getAnnotationSource(self):
        return self.annotationSource

    def _logSampleNameColumnDescription(self,source_column, sampleAnnotation, tumorAnnotation, normalAnnotation):
        result = {
            self.Source.SAMPLE: "Sample name is in the %s column." % sampleAnnotation,
            self.Source.TUMOR_NORMAL: "Sample name is the concatenation of the %s and %s columns." % (normalAnnotation, tumorAnnotation),
            self.Source.TUMOR: "Sample name is in the %s column." % tumorAnnotation,
            self.Source.NORMAL: "Sample name is in the %s column." % normalAnnotation,
            None: "Unable to figure out a sample name from the input columns."
        }[source_column]
        self.logger.info(result)

    def _getSourceColumn(self, sampleAnnotation, tumorAnnotation, normalAnnotation):
        """choose the right column to use as the sample annotation"""
        if sampleAnnotation is not None:
            source_column = self.Source.SAMPLE
        # Both, tumor and normal sample name annotations are present
        elif tumorAnnotation is not None and normalAnnotation is not None:
            source_column = self.Source.TUMOR_NORMAL
        elif tumorAnnotation is not None:
            source_column = self.Source.TUMOR
        # Only normal sample name is present
        elif normalAnnotation is not None:
            source_column = self.Source.NORMAL
        else:
            source_column = None
        return source_column

    def _getSampleNameGrabber(self, source_column, sampleAnnotation, tumorAnnotation, normalAnnotation):
        """returns the appropriate function to get the sample name"""
        sampleNameGenerators = {
            self.Source.SAMPLE: lambda mut: mut[sampleAnnotation],
            self.Source.TUMOR_NORMAL: lambda mut: "-".join([mut[normalAnnotation],mut[tumorAnnotation]]),
            self.Source.TUMOR: lambda mut: mut[tumorAnnotation],
            self.Source.NORMAL: lambda mut: mut[normalAnnotation],
            None: lambda mut: None
        }
        return sampleNameGenerators[source_column]

