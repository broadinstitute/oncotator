# LICENSE_GOES_HERE

import logging
import os
import re
import pysam
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.dbNSFP import dbnsfp_fieldnames


class dbNSFP(Datasource):
    """
    dbNSFP datasource expecting a directory of Tabix compressed and indexed tsv files.

    Example Tabix cmds:
        bgzip dbNSFP2.0_variant.chr1
        tabix -s 1 -b 2 -e 2 dbNSFP2.0_variant.chr1.gz

    This datasource requires the following annotations to be populated:
        variant_classification
        protein_change

    See the 'dbNSFP2.0.readme.txt' readme file in the datasource directory for field descriptions.
    """
    def __init__(self, src_dir, title='dbNSFP', version=None):
        self.title = title
        self.prot_regexp = re.compile('p\.([A-Z])\d+([A-Z])')

        super(dbNSFP, self).__init__(src_dir, title=self.title, version=version)
        self.db_obj = self._load_tabix_file_objs(src_dir)
        self.output_headers = [dbnsfp_fieldnames[14], dbnsfp_fieldnames[17]] + dbnsfp_fieldnames[21:] #only fields of interest
        self.logger = logging.getLogger(__name__)

    def annotate_mutation(self, mutation):
        if mutation['variant_classification'] == 'Missense_Mutation':
            chrn = mutation.chr
            start = mutation.start
            ref_allele, alt_allele = mutation.ref_allele, mutation.alt_allele
            regex_res = self.prot_regexp.match(mutation['protein_change'])

            if regex_res is None:
                logging.getLogger(__name__).error("Could not find a protein change for this missense mutation.")
                dbnsfp_results =[]
            else:
                ref_prot_allele, alt_prot_allele = regex_res.group(1), regex_res.group(2)
                dbnsfp_results = self._query_dbnsfp(chrn, start, start, self.db_obj)

            dbnsfp_data = None
            for res in dbnsfp_results:
                if all((res['ref'] == ref_allele, res['alt'] == alt_allele, res['aaref'] == ref_prot_allele, res['aaalt'] == alt_prot_allele)):
                    dbnsfp_data = res

            if dbnsfp_data is not None:
                for output_field in self.output_headers:
                    mutation.createAnnotation(output_field, dbnsfp_data[output_field], self.title)
            else:
                for output_field in self.output_headers:
                    mutation.createAnnotation(output_field, '', self.title)
        else:
            for output_field in self.output_headers:
                mutation.createAnnotation(output_field, '', self.title)

        return mutation

    def _load_tabix_file_objs(self, dbnsfp_dir):
        """ Returns dictionary of Tabixfile objects keyed by chromosome. """
        tabix_fnames = [f for f in os.listdir(dbnsfp_dir) if f.endswith('.gz')]
        tabix_objs = dict()
        for tabix_fname in tabix_fnames:
            contig = tabix_fname.replace('.gz', '').replace('dbNSFP2.0_variant.chr', '')
            tabix_objs[contig] = pysam.Tabixfile(os.path.join(dbnsfp_dir, tabix_fname))

        return tabix_objs

    def _get_dbnsfp_data_from_rows(self, dbnsfp_rows):
        return [dict(zip(dbnsfp_fieldnames, r.strip().split('\t'))) for r in dbnsfp_rows]

    def _query_dbnsfp(self, chromosome, start, end, tabix_objs):
        try:
            tabix_file = tabix_objs[chromosome]
        except KeyError:
            logging.getLogger(__name__).info("There is no dbNSFP info available for chromosome " + str(chromosome))
            return []
        query_str = '{}:{}-{}'.format(chromosome, start, end)
        res = tabix_file.fetch(region=query_str)
        return self._get_dbnsfp_data_from_rows(res)