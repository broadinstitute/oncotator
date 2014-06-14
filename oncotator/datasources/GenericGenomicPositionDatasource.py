# LICENSE_GOES_HERE

from oncotator.utils.db import get_db_data, get_binned_data, get_overlapping_records, get_summary_output_string
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.MutUtils import MutUtils


class GenericGenomicPositionDatasource(Datasource):
    """
    A datasource derived from a generic TSV file in which the first three columns are 'chr',
    'start', and 'end'.  All other columns will be used for annotation.

    use_binary
        if True, existing indexed binary will be used or created for future use.

    """
    def __init__(self, src_file, title='', version=None, use_binary=True, gpColumnNames="chr,start,end"):
        super(GenericGenomicPositionDatasource, self).__init__(src_file, title=title, version=version)

        index_mode = 'genomic_pos'
        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode, gpColumnNames)

    def annotate_mutation(self, mutation):
        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))

        if all([mutation.chr, mutation.start, mutation.end]):
            chrom, start, end = mutation.chr, mutation.start, mutation.end
            records = get_binned_data(self.db_obj, chrom, int(start), int(end))
            records = get_overlapping_records(records, int(start), int(end))
            if records:
                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)

        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)

        return mutation