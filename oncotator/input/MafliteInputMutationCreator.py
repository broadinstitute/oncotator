# LICENSE_GOES_HERE
from oncotator.Annotation import Annotation
from oncotator.Metadata import Metadata


"""
Created on Nov 9, 2012

@author: lichtens
"""

from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.MutationData import MutationData
from InputMutationCreator import InputMutationCreator
from MafliteMissingRequiredHeaderException import MafliteMissingRequiredHeaderException
import logging
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.MutUtils import MutUtils


class MafliteInputMutationCreator(InputMutationCreator):
    """
    A maflite file is a simple tsv file

    See the config file maflite_input.config for aliases and required headers.

    Additional columns can be included and will be annotate to the mutation using the header name.

    IMPORTANT NOTE: maflite will look at all aliases for alt_allele (see maflite_input.config) and choose the first that does not match the ref_allele
    """

    def __init__(self, filename, configFile='maflite_input.config', genomeBuild="hg19", other_options=None):
        """
        Constructor

        Currently, this InputCreator does not support any other options.  The parameter is ignored.

        """
        self.logger = logging.getLogger(__name__)

        self.config = ConfigUtils.createConfigParser(configFile)
        self._tsvReader = GenericTsvReader(filename)
        
        # Key is the required columns and the values are a list of valid alternative headers.
        # Key is column name to an alternative.
        self._alternativeDict = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config)
        self._reverseAlternativeDict = ConfigUtils.buildReverseAlternativeDictionary(self._alternativeDict)
        
        missingRequiredHeaders = []
        specifiedFields = self._tsvReader.getFieldNames()
        requiredColumns = sorted(list(MutationData.attributes))
        self._build = genomeBuild

        for col in requiredColumns:
            if col not in specifiedFields:
                isAltFound = False
                for alt in self._alternativeDict[col]:
                    if alt in specifiedFields:
                        isAltFound = True
                        break
                if not isAltFound:

                    # build is optional.
                    if col != "build":
                        missingRequiredHeaders.append(col)
        missingRequiredHeaders.sort()
        
        self.logger.info("Initializing a maflite file with the following header: " + str(self._tsvReader.getFieldNames()))
        if len(missingRequiredHeaders) > 0:
            raise MafliteMissingRequiredHeaderException("Specified maflite file (" + filename + ") missing required headers: " + ",".join(missingRequiredHeaders)  )

    def getComments(self):
        return self._tsvReader.getCommentsAsList()

    def getMetadata(self):
        result = Metadata()
        fieldNames = self._tsvReader.getFieldNames()
        fieldNameAliases = self._reverseAlternativeDict.keys()
        for fieldName in fieldNames:
            if fieldName in fieldNameAliases:
                fieldName = self._reverseAlternativeDict[fieldName]
            result[fieldName] = Annotation("", datasourceName="INPUT")
        return result

    def _find_alt_allele_in_other_field(self, raw_line_dict, ref_allele):
        """Check all the possible alt allele columns and choose the one that does not match the reference allele. """

        list_alternates = self._alternativeDict.get("alt_allele")

        for candidate_field in list_alternates:
            candidate_value = raw_line_dict.get(candidate_field, "")
            if candidate_value != "" and candidate_value != ref_allele:
                return candidate_value
        return ref_allele

    def createMutations(self):
        """ No inputs.
        Returns a generator of mutations built from the specified maflite file. """

        aliasKeys = self._reverseAlternativeDict.keys()
        allColumns = self._tsvReader.getFieldNames()

        for line in self._tsvReader:

            # We only need to assign fields that are mutation attributes and have a different name in the maflite file.
            mut = MutationData(build=self._build)

            for col in allColumns:
                # Three scenarios:
                #   1) col is name of mutation data field -- simple createAnnotation
                #   2) col name is an alias for a mutation data field -- do lookup then createAnnotation
                #   3) col name is not an alias for a mutation data field -- simple createAnnotation
                if col in aliasKeys:
                    realKey = self._reverseAlternativeDict[col]
                    self.logger.debug(realKey + " found from " + col)
                    val = line[col]
                    if realKey == "chr":
                        val = MutUtils.convertChromosomeStringToMutationDataFormat(line[col])
                    mut.createAnnotation(realKey, val, 'INPUT')
                else:
                    # Scenario 1 and 3
                    # Make sure to convert chromosome values.
                    val = line[col]
                    if col == "chr":
                        val = MutUtils.convertChromosomeStringToMutationDataFormat(line[col])
                    mut.createAnnotation(col, val, 'INPUT') 

            # if the alt allele == ref_allele, check that this is not a case where there is an alt_allele2 that is different.
            if mut.alt_allele == mut.ref_allele:
                mut.alt_allele = self._find_alt_allele_in_other_field(line, mut.ref_allele)

            # FIXME: Support more than one alias in the reverse dictionary.  Then this line can be removed.
            if mut.start is not "" and mut.end is "":
                mut.end = mut.start
            if mut.end is not "" and mut.start is "":
                mut.start = mut.end

            yield mut