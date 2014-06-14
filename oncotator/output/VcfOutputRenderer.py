# LICENSE_GOES_HERE

import collections
import logging
import os
import tempfile
import vcf

from oncotator.utils.MutUtils import MutUtils
from OutputRenderer import OutputRenderer
import oncotator.utils.ConfigUtils
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.output.RecordBuilder import RecordBuilder
from oncotator.output.OutputDataManager import OutputDataManager
from oncotator.config_tables.ConfigTableCreatorFactory import ConfigTableCreatorFactory
from oncotator.utils.OptionConstants import OptionConstants
import traceback


class VcfOutputRenderer(OutputRenderer):
    """
    The VcfOutputRenderer renders a vcf file from the given mutations.  All annotations are included with real names
    as column headers.

    Header is determined by the first mutation given.

    No attention is paid to order of the headers.
    """
    _vcfAnnotation = collections.namedtuple(typename="Annotation", field_names=["field", "ID", "num", "type"])

    def __init__(self, filename, configFile="vcf.out.config", otherOptions=None):
        """


        :param filename: output filename
        :param configFile: output config file
        """
        self._filename = filename
        self._otherOpts = dict() if otherOptions is None else otherOptions
        self.configFilename = configFile
        self.logger = logging.getLogger(__name__)
        self.config = oncotator.utils.ConfigUtils.ConfigUtils.createConfigParser(configFile, ignoreCase=False)
        self.chromHashCodeTable = None  # maps every chromosome in the mutations to a sortable integer
        self.configTableBuilder = ConfigTableCreatorFactory.getConfigTableCreatorInstance("output_vcf")
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.sampleNames = []  # all sample names in the mutations
        self.chroms = []  # all chromosomes in the mutations

    def renderMutations(self, mutations, metadata=None, comments=None):
        """
        Generate a simple tsv file based on the incoming mutations.
        Assumes that all mutations have the same annotations, even if some are not populated.

        :param mutations: generator of mutations to render
        :param metadata:
        :param comments:
        :return: returns a file name.
        """

        metadata = [] if metadata is None else metadata
        comments = [] if comments is None else comments

        self.logger.info("Rendering VCF output file: %s" % self._filename)

        # Initialize the config table
        self.configTable = self.configTableBuilder.getConfigTable(configFilename=self.configFilename)

        path = os.getcwd()

        dataManager = OutputDataManager(self.configTable, mutations, comments, metadata, path)
        self.sampleNames = dataManager.getSampleNames()

        # Write the header
        tempTemplateFile = tempfile.NamedTemporaryFile(dir=path, delete=False)
        filePointer = file(tempTemplateFile.name, 'w')
        header = dataManager.getHeader()
        filePointer.write(header)
        filePointer.close()

        # Sort the given TSV file
        sortedTempTsvFileName = dataManager.getSortedTsvFilename(path)

        try:
            inferGenotypes = self._otherOpts[OptionConstants.VCF_OUT_INFER_GENOTYPES]
        except (KeyError, TypeError):
            inferGenotypes = False

        self.logger.info("Render starting...")
        self._renderSortedTsv(tempTemplateFile.name, self._filename, sortedTempTsvFileName, self.sampleNames,
                              dataManager, inferGenotypes)

        # Remove template filename
        os.remove(tempTemplateFile.name)

        # Remove sorted tsv filename
        os.remove(sortedTempTsvFileName)

        self.logger.info("Rendered all mutations.")

        return self._filename

    def _isNewVcfRecordNeeded(self, curChrom, prevChrom, curPos, prevPos, curRefAllele, prevRefAllele):
        """
        Determines whether the current chromosome and position is the same as the previous chromosome and position.

        :param curChrom: current chromosome
        :param prevChrom: previous chromosome
        :param curPos: current position
        :param prevPos: previous position
        :param curRefAllele: current reference allele
        :param prevRefAllele: previous reference allele
        :return: true or false
        """
        if curChrom != prevChrom:
            return True
        if curPos != prevPos:
            return True
        if curRefAllele != prevRefAllele:
            return True
        return False

    def _renderSortedTsv(self, templateFilename, vcfFilename, tsvFilename, sampleNames, dataManager, inferGenotypes):
        """


        :param templateFilename:
        :param vcfFilename:
        :param tsvFilename:
        :param sampleNames:
        :param dataManager:
        """
        tempVcfReader = vcf.Reader(filename=templateFilename, strict_whitespace=True)
        pointer = file(vcfFilename, "w")
        vcfWriter = vcf.Writer(pointer, tempVcfReader, self.lineterminator)
        tsvReader = GenericTsvReader(tsvFilename, delimiter=self.delimiter)
        index = 0
        nrecords = 1000
        chrom = None
        pos = None
        refAllele = None
        recordBuilder = None

        ctr = 0
        m = None
        try:
            for m in tsvReader:
                ctr += 1
                isNewRecord = self._isNewVcfRecordNeeded(chrom, m["chr"], pos, m["start"], refAllele, m["ref_allele"])
                if isNewRecord:
                    if recordBuilder is not None:
                        record = recordBuilder.createRecord()
                        vcfWriter.write_record(record)
                        index += 1
                        if index % nrecords == 0:
                            self.logger.info("Rendered " + str(index) + " vcf records.")
                            vcfWriter.flush()

                    chrom = m["chr"]
                    if chrom.startswith("GL"):
                        chrom = "<" + chrom + ">"
                    pos = m["start"]
                    refAllele = m["ref_allele"]

                    recordBuilder = RecordBuilder(chrom, int(pos), refAllele, sampleNames)

                recordBuilder = self._parseRecordBuilder(m, recordBuilder, dataManager, inferGenotypes)

            if recordBuilder is not None:
                record = recordBuilder.createRecord()
                vcfWriter.write_record(record)
            vcfWriter.close()

        except Exception as e:
            self.logger.error(traceback.format_exc())
            self.logger.error("Error at mutation " + str(ctr) + " " + str([m["chr"], m["start"], m["end"]]) + ": ")

        self.logger.info("Rendered all " + str(index) + " vcf records.")

    def _parseRecordBuilder(self, m, recordBuilder, dataManager, inferGenotype):
        """
        Parse the input mutation object.
        First, this method

        :param m: mutation object
        :param recordBuilder:
        :param dataManager:
        :return:
        """
        idAnnotationNames = dataManager.getAnnotationNames("ID")
        qualAnnotationNames = dataManager.getAnnotationNames("QUAL")
        filterAnnotationNames = dataManager.getAnnotationNames("FILTER")
        infoAnnotationNames = dataManager.getAnnotationNames("INFO")
        formatAnnotationNames = dataManager.getAnnotationNames("FORMAT")
        sampleNameAnnotationNames = dataManager.getAnnotationNames("SAMPLE_NAME")

        if len(sampleNameAnnotationNames) != 0:
            sampleNameAnnotationName = sampleNameAnnotationNames[0]
        else:
            sampleNameAnnotationName = MutUtils.SAMPLE_NAME_ANNOTATION_NAME
        sampleName = m.get(sampleNameAnnotationName, None)

        altAllele = m["alt_allele"]

        recordBuilder.addAlt(altAllele)

        for name in idAnnotationNames:
            val = m.get(name, "")
            recordBuilder.addID(val)

        if len(qualAnnotationNames) == 0:
            qual = "qual"
        else:
            qual = qualAnnotationNames[0]
        recordBuilder.addQual(m.get(qual, "."))

        for name in filterAnnotationNames:
            ID = dataManager.getFieldID(name)
            val = m.get(name, "")
            recordBuilder.addFilter(ID, val)

        for name in infoAnnotationNames:
            annotation = dataManager.getOutputAnnotation(name)
            ID = annotation.getID()
            num = annotation.getNumber()
            dataType = annotation.getDataType()
            isSplit = annotation.isSplit()
            val = m.get(name, "")
            recordBuilder.addInfo(sampleName, ID, num, dataType, val, isSplit)

        for name in formatAnnotationNames:
            annotation = dataManager.getOutputAnnotation(name)
            ID = annotation.getID()
            num = annotation.getNumber()
            dataType = annotation.getDataType()
            isSplit = annotation.isSplit()
            val = m.get(name, "")
            if num == 0 or dataType == "Flag":
                msg = "%s is of data type Flag. Only Integer, Float, Character, and String data types are permissible" \
                      " in the Format field." % name
                logging.warn(msg)
            else:
                recordBuilder.addFormat(sampleName, ID, num, dataType, val, isSplit, inferGenotype)

        return recordBuilder