# LICENSE_GOES_HERE

import copy
import logging
import string
import vcf
from oncotator.utils.TagConstants import TagConstants
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.MutUtils import MutUtils
import operator


class IndexedVcfDatasource(Datasource):
    """
    A datasource derived from a VCF file.  Expects a bgzipped vcf using Tabix.

    Instructions on how to index file using Tabix prior Oncotator:

    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.


    The following VCF columns will be added to the output.
        REF, ALT and INFO


    Multiple annotation data for the same mutation will be delimited by "|".

    TODO: Support for FORMAT and sample columns in output

    """
    def __init__(self, src_file, title='', version=None, match_mode="exact"):
        # all of the info is coming from config file
        super(IndexedVcfDatasource, self).__init__(src_file, title=title, version=version)
        self.vcf_reader = vcf.Reader(filename=src_file, strict_whitespace=True)
        self.vcf_info_headers = self.vcf_reader.infos.keys()
        # Annotation name from INFO columns
        self.output_vcf_headers = dict([(ID, '_'.join([self.title, ID])) for ID in self.vcf_info_headers])
        # Type information
        self.output_vcf_types = dict([(ID, self.vcf_reader.infos[ID].type) for ID in self.vcf_info_headers])
        # Number information
        self.output_vcf_nums = dict([(ID, self.vcf_reader.infos[ID].num) for ID in self.vcf_info_headers])
        # Descriptions
        self.output_vcf_descs = dict([(ID, self.vcf_reader.infos[ID].desc) for ID in self.vcf_info_headers])
        self.match_mode = match_mode

    def _determine_tags(self):
        """
        Determine tag constants associated with annotation IDs.

        :return: map of tag IDs and respective tag constants
        """
        tagsMap = {}
        for ID in self.vcf_info_headers:
            num = self.output_vcf_nums[ID]
            if num == -1:  # the num field decides whether to split or not
                tagsMap[ID] = [TagConstants.INFO, TagConstants.SPLIT]
            else:
                tagsMap[ID] = [TagConstants.INFO, TagConstants.NOT_SPLIT]
        return tagsMap

    def _determine_info_annotation_value(self, vcf_record, ID, alt_index):
        """
        Determine the appropriate value that corresponds to a given ID.
        This method checks whether ID is present in the input Vcf record or not. If it is not found, then an appropriate
        missing value is computed.

        :param alt_index:
        :param vcf_record: input Vcf record
        :param ID: ID
        :return: value that corresponds to a given ID
        """
        if alt_index is None:
            val = self._determine_missing_value(ID)
            return val

        if ID in vcf_record.INFO:
            vals = vcf_record.INFO[ID]
            num = self.output_vcf_nums[ID]

            if isinstance(vals, list):
                if num == -1 and alt_index >= 0 and alt_index is not None:  # split by alternative allele
                    val = [vals[alt_index]]
                else:
                    val = vals
            else:  # force it to be a list
                val = [vals]
        else:
            val = self._determine_missing_value(ID)

        return val

    def _determine_missing_value(self, ID):
        """
        Determine the missing value that corresponds to a given ID.
        This method takes into account the number field and computes an appropriate missing value for a given ID.

        :param ID: ID
        :return: value that corresponds to a given ID
        """
        nsamples = len(self.vcf_reader.samples)
        num = self.output_vcf_nums[ID]
        if num == -2:
            return [""]*nsamples
        elif num == -1:
            return [""]
        elif num == 0:
            return [False]
        elif num is None:
            return [""]

        return [""]*num

    def _determine_matching_alt_indices(self, mut, record, build):
        """

        :param mut:
        :param record:
        :return:
        """
        indices = []
        if record.is_monomorphic:
            chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
            startPos = record.POS
            endPos = record.POS
            ref_allele = record.REF

            if self.match_mode == "exact":
                if mut.chr == chrom and mut.ref_allele == ref_allele:
                    indices = [-1]
            else:
                if mut.chr == chrom and int(mut.start) <= startPos and int(mut.end) >= endPos:
                    indices = [-1]
        else:
            # Iterate over all alternates in the record
            for index in xrange(0, len(record.ALT)):
                chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
                startPos = record.POS
                endPos = record.POS
                ref = str(record.REF)
                alt = str(record.ALT[index])
                ds_mut = MutUtils.initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build)

                if self.match_mode == "exact":
                    if mut.chr == ds_mut.chr and mut.ref_allele == ds_mut.ref_allele \
                        and mut.alt_allele == ds_mut.alt_allele and int(mut.start) == int(ds_mut.start) \
                        and int(mut.end) == int(ds_mut.end):
                        indices += [index]
                else:  # cases whether the match mode isn't exact
                    if mut.chr == ds_mut.chr and int(mut.start) == int(ds_mut.start) and int(mut.end) == int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) >= int(ds_mut.start) \
                        and int(mut.end) >= int(ds_mut.end) and int(mut.start) <= int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) <= int(ds_mut.start) and int(mut.end) >= int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) <= int(ds_mut.start) \
                        and int(mut.end) <= int(ds_mut.end) and int(mut.end) >= int(ds_mut.start):
                        indices += [index]

        # if len(indices) == 0:
        #     indices = [None]

        return indices

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from Tabix indexed Vcf file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """

        # Get all the records corresponding to the given mutation
        mut_start = int(mutation.start)
        mut_end = int(mutation.end)
        vals = dict()
        record = None

        try:
            vcf_records = self.vcf_reader.fetch(mutation.chr, mut_start - 1, mut_end)  # query database for records
        except ValueError as ve:
            logging.getLogger(__name__).debug("Exception when looking for vcf records. Empty set of records being "
                                              "returned: " + repr(ve))
        else:
            # Process values.
            build = "hg19"
            for record in vcf_records:
                alt_indices = self._determine_matching_alt_indices(mutation, record, build)
                for ID in self.vcf_info_headers:
                    for alt_index in alt_indices:
                        val = self._determine_info_annotation_value(record, ID, alt_index)
                        if ID not in vals:  # dictionary of list of lists
                            vals[ID] = [val]
                        else:
                            vals[ID] += [val]

        if len(vals) == 0:
            for ID in self.vcf_info_headers:
                vals[ID] = [[""]]

        if record is None:
            msg = "Exception when looking for tsv records for chr%s:%s-%s. " \
                  "Empty set of records being returned." % (mutation.chr, mutation.start, mutation.end)
            logging.getLogger(__name__).debug(msg)

        tags = self._determine_tags()
        for ID in self.vcf_info_headers:
            if ID not in vals:  # this happens in cases where no matching records are found
                val = self._determine_missing_value(ID)
                vals[ID] = [val]

            if self.match_mode == "exact":
                # TODO: monomorphic records?
                val = []
                for v in vals[ID]:  # list of lists
                    v = string.join([str(v) if v is not None else "" for v in v], ",")
                    if v:
                        val += [v]
                val = string.join(val, "|")
            elif self.match_mode == "avg":
                if self.output_vcf_types[ID] in ("Integer", "Float"):  # integers and floats
                    num = self.output_vcf_nums[ID]
                    val = ""
                    if num in (None, -1, 0, 1,):
                        v = reduce(operator.add, vals[ID])  # flattens list of lists to a list
                        v = [v for v in v if v is not None and v != ""]  # discard non-numerics
                        if len(v) > 1:
                            val = str(float(sum(v))/len(v))
                        elif len(v) == 1:
                            val = str(float(v[0]))
                    else:
                        if len(vals[ID][0]) != num:  # addresses cases where no records were found
                            val = [""]
                        else:
                            nvals = len(vals[ID])  # number of lists
                            val = [[]]*num
                            for i in xrange(num):  # iterate over all numbers
                                for j in xrange(nvals):  # iterate over each list
                                    v = vals[ID][j][i]
                                    if v not in (None, "", "."):
                                        val[i] = val[i] + [v]
                                if len(val[i]) >= 1:
                                    val[i] = str(float(sum(val[i]))/len(val[i]))
                                else:
                                    val[i] = ""
                        val = string.join(val, ",")
                    self.output_vcf_types[ID] = "Float"
                elif self.output_vcf_types[ID] == "Character":
                    val = []
                    for v in vals[ID]:
                        val += [map(str, [string.join([v if v is not None else "" for v in v], ",")])]
                    val = string.join(val, "|")
                    self.output_vcf_types[ID] = "String"
                else:
                    val = []
                    for v in vals[ID]:
                        v = string.join(map(str, [v if v is not None else "" for v in v]), ",")
                        if v:
                            val += [v]
                    val = string.join(val, "|")
                    self.output_vcf_types[ID] = "String"
                    if self.output_vcf_nums[ID] == 0:
                        self.output_vcf_nums[ID] = None
            else:
                self.output_vcf_types[ID] = "String"  # everything is a String
                if self.output_vcf_nums[ID] == 0:
                    self.output_vcf_nums[ID] = None
                # values are delimited by pipe, and thus, the type is forced to be a String
                val = []
                for v in vals[ID]:  # list of lists
                    v = string.join([str(v) if v is not None else "" for v in v], ",")
                    if v:
                        val += [v]
                val = string.join(val, "|")

            mutation.createAnnotation(self.output_vcf_headers[ID], val, self.title, self.output_vcf_types[ID],
                                      self.output_vcf_descs[ID], tags=copy.copy(tags[ID]),
                                      number=self.output_vcf_nums[ID])

        return mutation