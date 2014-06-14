# LICENSE_GOES_HERE

import logging
import pysam
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.TagConstants import TagConstants
import string
from oncotator.utils.MutUtils import MutUtils
from oncotator.MutationData import MutationData


class IndexedTsvDatasource(Datasource):
    """
    A datasource derived from a TSV file.  Expects a bgzipped vcf using Tabix.

    Instructions on how to index file using Tabix prior Oncotator:

    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.

    Multiple annotation data for the same mutation will be delimited by "|".
    """
    def __init__(self, src_file, title, version, colNames, indexColNames, annotationColNames, match_mode, colDataTypes):
        super(IndexedTsvDatasource, self).__init__(src_file, title=title, version=version)
        self.tsv_reader = pysam.Tabixfile(filename=src_file)  # initialize the tsv reader
        # Index the column names in tsv header
        self.tsv_headers = dict([(colName, index) for (index, colName) in enumerate(colNames)])
        # Annotation column names are the only columns that are used for annotating a mutation object.
        self.output_tsv_headers = dict([(colName, "_".join([self.title, colName])) for colName in annotationColNames])
        self.tsv_index = dict()
        self.tsv_index["chrom"] = self.tsv_headers[indexColNames[0]]
        self.tsv_index["start"] = self.tsv_headers[indexColNames[1]]
        self.tsv_index["end"] = self.tsv_headers[indexColNames[2]]
        if len(indexColNames) == 5:
            self.tsv_index["ref"] = self.tsv_headers[indexColNames[3]]
            self.tsv_index["alt"] = self.tsv_headers[indexColNames[4]]

        self.dataTypes = colDataTypes
        self.match_mode = match_mode

    def _is_matching(self, mut, tsv_record):

        chrom = tsv_record[self.tsv_index["chrom"]]
        startPos = tsv_record[self.tsv_index["start"]]
        endPos = tsv_record[self.tsv_index["end"]]
        build = "hg19"

        if self.match_mode == "exact":
            if "ref" in self.tsv_index and "alt" in self.tsv_index:  # ref and alt information is present
                ref = tsv_record[self.tsv_index["ref"]]
                alt = tsv_record[self.tsv_index["alt"]]
                if ref == "-" or alt == "-":  # addresses Mutation Annotation Format based tsv records
                    ds_mut = MutationData(chrom, startPos, endPos, ref, alt, build)
                else:  # addresses tsv records where the input isn't a Mutation Annotation Format file
                    ds_mut = MutUtils.initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build)

                if mut.chr == ds_mut.chr and mut.ref_allele == ds_mut.ref_allele \
                    and mut.alt_allele == ds_mut.alt_allele and int(mut.start) == int(ds_mut.start) \
                    and int(mut.end) == int(ds_mut.end):
                    return True
            else:  # do not use ref and alt information
                if mut.chr == chrom and int(mut.start) == int(startPos) and int(mut.end) == int(endPos):
                    return True
        else:
            if mut.chr == chrom and int(mut.start) == int(startPos) and int(mut.end) == int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) >= int(startPos) and int(mut.end) >= int(endPos) \
                and int(mut.start) <= int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) <= int(startPos) and int(mut.end) >= int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) <= int(startPos) and int(mut.end) <= int(endPos) \
                and int(mut.end) >= int(startPos):
                return True

        return False

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from a Tabix indexed TSV file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """
        vals = {}
        mut_start = int(mutation.start)
        mut_end = int(mutation.end)
        try:
            # tabix needs position - 1
            tsv_records = self.tsv_reader.fetch(mutation.chr, mut_start - 1, mut_end, parser=pysam.asTuple())
            for tsv_record in tsv_records:
                if not tsv_record:  # skip in case no records are found
                    continue

                # Determine whether the new tsv record matches mutation or not
                if self._is_matching(mutation, tsv_record):
                    for colName in self.output_tsv_headers:
                        val = tsv_record[self.tsv_headers[colName]]
                        if colName not in vals:
                            vals[colName] = [val]
                        else:
                            vals[colName] += [val]

        except ValueError as ve:
            msg = "Exception when looking for tsv records. Empty set of records being returned: " + repr(ve)
            logging.getLogger(__name__).debug(msg)

        for colName in self.output_tsv_headers:
            if colName in vals:  # this case happens when there are no matching records to be found
                val = string.join(vals[colName], "|")
            else:
                val = ""

            ds_type = self.dataTypes[colName]
            if self.match_mode == "exact":
                if "ref" not in self.tsv_index or "alt" not in self.tsv_index:
                    ds_type = "String"
            elif self.match_mode == "avg":
                if ds_type == "Integer" or ds_type == "Float":
                    ds_type = "Float"
                    if len(vals) != 0:
                        try:
                            v = [float(val) for val in vals[colName] if val not in ("", ".", "-",)]
                            if len(v) != 0:
                                val = str(sum(v)/len(v))
                            else:
                                val = ""
                        except ValueError:  # in cases where it's unsuccessful, default to overlap behavior
                            msg = "Exception when trying to cast %s value to float." % colName
                            logging.getLogger(__name__).warn(msg)
                            val = string.join(map(str, vals[colName]), "|")
                            ds_type = "String"
                else:
                    ds_type = "String"
            else:  # default: overlap
                ds_type = "String"

            if len(vals) == 0:
                msg = "Exception when looking for tsv records for chr%s:%s-%s. " \
                      "Empty set of records being returned." % (mutation.chr, mutation.start, mutation.end)
                logging.getLogger(__name__).debug(msg)

            mutation.createAnnotation(self.output_tsv_headers[colName], val, self.title, annotationDataType=ds_type,
                                      tags=[TagConstants.INFO, TagConstants.NOT_SPLIT])

        return mutation