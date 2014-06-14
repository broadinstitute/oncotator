# LICENSE_GOES_HERE
import re
import logging
import string
import copy

import vcf

from oncotator.Metadata import Metadata
from oncotator.input.InputMutationCreator import InputMutationCreatorOptions
from oncotator.utils.TagConstants import TagConstants
from InputMutationCreator import InputMutationCreator
from oncotator.utils.MutUtils import MutUtils
from oncotator.Annotation import Annotation
from oncotator.config_tables.ConfigTableCreatorFactory import ConfigTableCreatorFactory


class VcfInputMutationCreator(InputMutationCreator):
    """
    TODO: Finish documentation
    
    Convert a VCF 4.1 into a MutationData generator.
            Adds the following annotations (as INPUT annotation source):
                sampleName -- the name of the sample in the VCF
                altAlleleSeen -- whether the alternate allele was seen in the mutation.  This is whether a "1" appears in the GT field.
                
    
    """

    def __init__(self, filename, configFile='vcf.in.config', genomeBuild="hg19", other_options=None):
        """

        :param filename:
        :param configFile:
        """
        self.filename = filename
        self.build = genomeBuild
        self.configFilename = configFile
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)
        self.configTableBuilder = ConfigTableCreatorFactory.getConfigTableCreatorInstance("input_vcf")
        self.isTagSplit = dict()
        self.logger = logging.getLogger(__name__)

        if other_options is None:
            other_options = {}

        self._is_skipping_no_alts = other_options.get(InputMutationCreatorOptions.IS_SKIP_ALTS, False)

    def _addGenotypeData2Mutation(self, mutation, record, index):
        """


        :param mutation: input mutation object
        :param record:
        :param index:
        """
        IDs = self.configTable.getFormatFieldIDs()
        genotypeData = None

        if len(IDs) != 0:
            sampleNameAnnotation = mutation.getAnnotation(MutUtils.SAMPLE_NAME_ANNOTATION_NAME)
            sampleName = sampleNameAnnotation.getValue()
            genotypeData = record.genotype(sampleName)

        if record.FORMAT is not None:
            for ID in IDs:
                val = ""
                dataType = self.vcf_reader.formats[ID].type

                name = self.configTable.getFormatFieldName(ID)
                num = self.vcf_reader.formats[ID].num

                tags = []

                if (genotypeData is not None) and (hasattr(genotypeData.data, ID)):
                    if name not in self.isTagSplit:
                        isTagSplit = self._determineIsSplit(ID, num, "FORMAT")
                        self.isTagSplit[name] = isTagSplit
                    else:
                        isTagSplit = self.isTagSplit[name]

                    if isTagSplit and isinstance(genotypeData[ID], list):
                        val = genotypeData[ID][index]
                    else:
                        val = genotypeData[ID]

                    if isTagSplit:
                        tags = [TagConstants.FORMAT, TagConstants.SPLIT]
                    else:
                        tags = [TagConstants.FORMAT, TagConstants.NOT_SPLIT]

                if (val is None) or (val == ""):
                    if dataType == "Flag":
                        val = "False"
                    else:
                        val = ""
                elif isinstance(val, list):
                    val = string.join(["" if v is None else str(v) for v in val], ",")
                else:
                    val = str(val)
                if name in mutation:
                    name = string.join(words=[name, "__FORMAT__"], sep="")
                mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.formats[ID].desc, tags=tags,
                                          number=num)

        return mutation

    def _addInfoData2Mutation(self, mutation, record, index):
        """
        This method

        :param mutation:
        :param record:
        :param index:
        :return:
        """
        IDs = self.configTable.getInfoFieldIDs()

        for ID in IDs:
            val = ""
            dataType = self.vcf_reader.infos[ID].type

            num = self.vcf_reader.infos[ID].num
            name = self.configTable.getInfoFieldName(ID)
            tags = []

            if ID in record.INFO:
                if name not in self.isTagSplit:
                    isTagSplit = self._determineIsSplit(ID, num, "INFO")
                    self.isTagSplit[name] = isTagSplit
                else:
                    isTagSplit = self.isTagSplit[name]

                if isTagSplit and isinstance(record.INFO[ID], list):
                    val = record.INFO[ID][index]
                else:
                    val = record.INFO[ID]

                if isTagSplit:
                    tags = [TagConstants.INFO, TagConstants.SPLIT]
                else:
                    tags = [TagConstants.INFO, TagConstants.NOT_SPLIT]

            if (val is None) or (val == ""):
                if dataType == "Flag":
                    val = "False"
                else:
                    val = ""
            elif isinstance(val, list):
                val = string.join(["" if v is None else str(v) for v in val], ",")
            else:
                val = str(val)

            mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.infos[ID].desc, tags=tags,
                                      number=num)

        return mutation

    def _determineIsSplit(self, ID, num, fieldType):
        """

        :param ID:
        :param num:
        :param fieldType:
        :return:
        """
        if num == -2:  # by the number of samples
            if fieldType == "FORMAT":
                if self.configTable.isFieldIDInSplitSet(fieldType, ID):
                    return True
            else:
                return False
        elif num == -1:  # by the number of alternates
            if self.configTable.isFieldIDInNotSplitSet(fieldType, ID):  # override the default using the config file
                return False
            else:
                return True
        elif num == 0:
            return False
        elif num is None:
            if self.configTable.isFieldIDInSplitSet(fieldType, ID):  # override the default using the config file
                return True
            else:
                return False

        return False

    def _determineAltSeen(self, sample_gt_str, index):
        """Look at the genotype string to see if the alternate is present.  GT of ./. is considered a 'yes'"""
        # TODO: Confirm that it is necessary to take into account the index.  Take into account the index.
        is_alt_seen = "True"
        if sample_gt_str is not None:
            # Split genotype field into number of entries as ploidy.  Using chars '/' or '|'
            gt_haploid_list = re.split('/|\|', sample_gt_str)
            if all([haploid != str(index) for haploid in gt_haploid_list]):
                is_alt_seen = "False"
        else:
            # GT is ./.
            is_alt_seen = "False"
        return is_alt_seen

    def createMutations(self):
        """ Creates a mutation for each mutation by each sample, regardless of allelic depth, etc.

            Annotations that this will generate (as source = "INPUT"):
                sampleName
                isCalled
                isAltAllele
                allelic_depth -- DP in the vcf

            TODO: Complete documentation
        """
        self.configTable = self.configTableBuilder.getConfigTable(filename=self.filename,
                                                                  configFilename=self.configFilename)
        build = self.build
        is_tn_vcf_warning_delivered = False
        for record in self.vcf_reader:
            for index in range(len(record.ALT)):
                if len(record.samples) <= 0:
                    yield self._createMutation(record, index, build)
                else:
                    sampleRecList = record.samples
                    sample_names = [s.sample for s in sampleRecList]
                    is_tumor_normal_vcf = "NORMAL" in sample_names and len(sample_names) == 2
                    if not is_tn_vcf_warning_delivered and is_tumor_normal_vcf:
                        is_tn_vcf_warning_delivered = True
                        logging.getLogger(__name__).warn("Tumor-Normal VCF detected.  The Normal will assume GT= 0/0, unless GT field specified otherwise.")

                    for sample in sampleRecList:

                        sample_name = sample.sample

                        #TODO: Confirm that alt_allele_seen will be False in all cases of GT = ./.
                        genotype = "GT"
                        is_alt_seen = "True"
                        if genotype in sample.data._fields:
                            is_alt_seen = self._determineAltSeen(sample.data.GT, index + 1)

                        # HACK: If the sample name is NORMAL, there is more than one sample, and
                        # there is no GT field (or GT is ./.) then assume that this is alt_allele_seen of False
                        if is_tumor_normal_vcf and sample_name == "NORMAL" and (genotype not in sample.data._fields):
                            is_alt_seen = "False"

                        # Check to see if we should even render this mutation at all
                        if self._is_skipping_no_alts and not (is_alt_seen == "True"):
                            continue

                        sampleMut = self._createMutation(record, index, build)

                        if is_tumor_normal_vcf and sample_name != "NORMAL":
                            sampleMut.createAnnotation("tumor_barcode", sample_name, "INPUT")
                        sampleMut.createAnnotation(MutUtils.SAMPLE_NAME_ANNOTATION_NAME, sample_name, "INPUT")

                        sampleMut["alt_allele_seen"] = is_alt_seen
                        sampleMut = self._addGenotypeData2Mutation(sampleMut, record, index)

                        yield sampleMut

    def _createMutation(self, record, alt_index, build):
        chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
        startPos = int(record.POS)
        endPos = int(record.POS)
        ref = record.REF
        ref = "" if ref == "." else ref

        alt = ref
        if not record.is_monomorphic:
            alt = str(record.ALT[alt_index])

        mut = MutUtils.initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build)
        ID = "" if record.ID is None else record.ID
        mut.createAnnotation("id", ID, "INPUT", tags=[TagConstants.ID])
        mut.createAnnotation("qual", str(record.QUAL), "INPUT", tags=[TagConstants.QUAL])
        mut.createAnnotation("alt_allele_seen", str(True), "INPUT")
        mut = self._addFilterData2Mutation(mut, record)
        mut = self._addInfoData2Mutation(mut, record, alt_index)
        return mut

    def _addFilterData2Mutation(self, mut, record):
        for flt in self.vcf_reader.filters:  # for each filter in the header
            description = self.vcf_reader.filters[flt].desc  # parse the description
            if record.FILTER is None:  # information is missing
                mut.createAnnotation(flt, "", "INPUT", annotationDescription=description,
                                     tags=[TagConstants.FILTER])
            elif flt in record.FILTER:  # filters the site failed are listed
                mut.createAnnotation(flt, "FAIL", "INPUT", annotationDescription=description,
                                     tags=[TagConstants.FILTER])
            else:  # site passed all filters
                mut.createAnnotation(flt, "PASS", "INPUT", annotationDescription=description,
                                     tags=[TagConstants.FILTER])
        return mut

    def reset(self):
        """ Resets the internal state, so that mutations can be generated. """
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)

    def _parseMiscellaneousMetadata(self):
        """


        :return:
        """
        comments = []
        keys = self.vcf_reader.metadata.keys()
        for key in keys:
            vals = self.vcf_reader.metadata[key]
            if isinstance(vals, list):
                for itm in vals:
                    if isinstance(itm, dict):
                        val = [string.join([str(k), str(v)], "=") for k, v in itm.iteritems()]
                        val = string.join(val, ",")
                        val = string.join(["<", ">"], val)
                    else:
                        val = str(itm)

                    comment = string.join([key, val], "=")
                    comment = string.join(["##", comment], "")
                    comments.append(comment)
            elif isinstance(vals, dict):
                val = [string.join([str(k), str(v)], "=") for k, v in vals.iteritems()]
                val = string.join(val, ",")
                val = string.join(["<", ">"], val)

                comment = string.join([key, val], "=")
                comment = string.join(["##", comment], "")
                comments.append(comment)
            else:
                val = str(vals)
                comment = string.join([key, val], "=")
                comment = string.join(["##", comment], "")
                comments.append(comment)
        return comments

    def _parseContigsMetadata(self):
        comments = []
        keys = self.vcf_reader.contigs.keys()
        for key in keys:
            val = self.vcf_reader.contigs[key]
            ID = val.id
            length = str(val.length)
            val = string.join(["##contig=<ID=", ID, ",length=", length, ">"], "")
            comments.append(val)
        return comments

    def _parseAltsMetadata(self):
        comments = []
        keys = self.vcf_reader.alts.keys()
        for key in keys:
            val = self.vcf_reader.alts[key]
            ID = val.id
            desc = val.desc
            val = string.join(["##ALT=<ID=", ID, ",Description=\"", desc, "\">"], "")
            comments.append(val)
        return comments

    def getComments(self):
        """ Comments often need to be passed into the output.  Get the comments from the input file."""
        comments = self._parseMiscellaneousMetadata()
        comments += self._parseContigsMetadata()
        comments += self._parseAltsMetadata()
        return comments

    def _addFormatFields2Metadata(self, metadata):
        """

        :param metadata:
        :return:
        """
        for ID in self.configTable.getFormatFieldIDs():
            name = self.configTable.getFormatFieldName(ID)
            num = self.vcf_reader.formats[ID].num
            tags = [TagConstants.FORMAT]
            isSplitTag = self._determineIsSplit(ID, num, "FORMAT")
            if isSplitTag:
                tags += [TagConstants.SPLIT]
            if name in metadata:
                name = string.join(words=[name, "__FORMAT__"], sep="")
            metadata[name] = Annotation("", "INPUT", self.vcf_reader.formats[ID].type, self.vcf_reader.formats[ID].desc,
                                        tags=tags, number=num)
        return metadata

    def _addInfoFields2Metadata(self, metadata):
        """
        Add INFO field meta-information to metadata.

        :param metadata:
        :return:
        """
        for ID in self.configTable.getInfoFieldIDs():
            name = self.configTable.getInfoFieldName(ID)
            num = self.vcf_reader.infos[ID].num
            tags = [TagConstants.INFO]
            isSplitTag = self._determineIsSplit(ID, num, "INFO")
            if isSplitTag:
                tags += [TagConstants.SPLIT]
            metadata[name] = Annotation("", "INPUT", self.vcf_reader.infos[ID].type, self.vcf_reader.infos[ID].desc,
                                        tags=tags, number=num)
        return metadata

    def _addFilterFields2Metadata(self, metadata):
        """


        :param metadata:
        :return: modified metadata
        """
        for filt in self.vcf_reader.filters:  # for each filter in the header
            metadata[filt] = Annotation("", "INPUT", "String", self.vcf_reader.filters[filt].desc,
                                        tags=[TagConstants.FILTER])
        return metadata

    def _createMetadata(self):
        """


        :return: metadata
        """
        self.configTable = self.configTableBuilder.getConfigTable(filename=self.filename,
                                                                  configFilename=self.configFilename)
        metadata = Metadata()
        metadata = self._addFilterFields2Metadata(metadata)
        metaData = self._addInfoFields2Metadata(metadata)
        metadata = self._addFormatFields2Metadata(metadata)

        metaData["id"] = Annotation("", "INPUT", "String", "", [TagConstants.ID])
        metaData["qual"] = Annotation("", "INPUT", "String", "", [TagConstants.QUAL])

        return metadata

    def getMetadata(self):
        metadata = self._createMetadata()
        return metadata
