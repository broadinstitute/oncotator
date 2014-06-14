# LICENSE_GOES_HERE

import logging
import vcf
from oncotator.datasources.Datasource import Datasource


class dbSNP(Datasource):
    """
    dbSNP Datasource expecting a Tabix compressed and indexed dbsnp vcf file.

    Output columns are:

        dbSNP_RS
            RS id of overlapping dbSNPs.  Multiple records are delimited by '|'
        dbSNP_Val_Status
            'byFrequency' if dbSNP is validated by frequency.
            'b10000genomes' if dbSNP is found in 1000genomes dataset.

    """
    def __init__(self, src_file, title='dbSNP', version=None):
        self.title = title

        super(dbSNP, self).__init__(src_file, title=self.title, version=version)
        self.vcf_reader = vcf.Reader(filename=src_file, strict_whitespace=True)
        self.output_headers = ['dbSNP_RS', 'dbSNP_Val_Status']
        self.logger = logging.getLogger(__name__)

    def annotate_mutation(self, mutation):
        chrom, start, end = mutation.chr, int(mutation.start), int(mutation.end)

        final_values = {}
        for header in self.output_headers:
            final_values[header] = set()

        overlapping_vcf_records = []
        try:
            overlapping_vcf_records = self.vcf_reader.fetch(chrom, start, end)
        except ValueError as ve:
            self.logger.debug("Exception when looking for vcf records.  Empty set of records being returned: " + repr(ve))

        for vcf_rec in overlapping_vcf_records:
            for header in self.output_headers:
                final_values['dbSNP_RS'].add(vcf_rec.ID)
                if 'VLD' in vcf_rec.INFO:
                    final_values['dbSNP_Val_Status'].add('byFrequency')
                if 'KGVAL' in vcf_rec.INFO or 'KGPROD' in vcf_rec.INFO or 'KGPilot1' in vcf_rec.INFO or 'KGPilot123' in vcf_rec.INFO:
                    final_values['dbSNP_Val_Status'].add('by1000genomes')

        for header in self.output_headers:
            mutation.createAnnotation(header, '|'.join(final_values[header]), self.title)

        return mutation