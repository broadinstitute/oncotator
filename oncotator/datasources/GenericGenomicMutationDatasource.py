# LICENSE_GOES_HERE

from oncotator.datasources.GenericGenomicPositionDatasource import GenericGenomicPositionDatasource
from oncotator.utils.db import get_binned_data, get_overlapping_records, get_summary_output_string


class GenericGenomicMutationDatasource(GenericGenomicPositionDatasource):
    """
    This datasource extends Generic_GenomicPosition_DataSource to also match on ref_allele
    and alt_allele columns.  All other columns will be used for annotation.
    """
    def __init__(self, src_file, title='', version=None, use_complementary_strand_alleles_for_negative_strand_transcripts=False, **kwargs):
        super(GenericGenomicMutationDatasource, self).__init__(src_file, title=title, version=version, **kwargs)

        self.ref_allele_fieldname, self.alt_allele_fieldname = None, None
        for header in self.output_headers:
            if header.endswith('_ref_allele'):
                self.ref_allele_fieldname = header
            if header.endswith('_alt_alele'):
                self.alt_allele_fieldname = header

        if self.ref_allele_fieldname is None or self.alt_allele_fieldname is None:
            raise Exception('Unable to determine ref_allele or alt_allele column name.')

        self.output_headers = [h for h in self.output_headers if not h.endswith('_ref_allele') and not h.endswith('_alt_alele')]
        self.use_complementary_strand_alleles_for_negative_strand_transcripts = use_complementary_strand_alleles_for_negative_strand_transcripts

    def annotate_mutation(self, mutation):
        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))

        if all(field in mutation for field in ['chr','start','end', 'ref_allele', 'alt_allele', 'strand']):
            chr, start, end = mutation.chr, mutation.start, mutation.end
            ref_allele, alt_allele = str(mutation.ref_allele), str(mutation.alt_allele) #changed to str incase vcf.model._Substitution object is being used
            if self.use_complementary_strand_alleles_for_negative_strand_transcripts:
                if mutation.get('strand') == '-':
                    ref_allele, alt_allele = Seq.reverse_complement(ref_allele), Seq.reverse_complement(alt_allele)

            records = get_binned_data(self.db_obj, chr,int(start), int(end))
            records = get_overlapping_records(records, int(start), int(end))
            records = [rec for rec in records if rec[self.ref_allele_fieldname] == ref_allele and rec[self.alt_allele_fieldname] == alt_allele]
            if records:
                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)

        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)

        return mutation
