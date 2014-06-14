# LICENSE_GOES_HERE

import logging
import re
from oncotator.MissingAnnotationException import MissingAnnotationException
from oncotator.datasources.GenericGenomicPositionDatasource import GenericGenomicPositionDatasource
from oncotator.utils.db import get_db_data, get_binned_data, get_overlapping_records, get_summary_output_string


class GenericGeneProteinPositionDatasource(GenericGenomicPositionDatasource):
    """ For annotating protein positions and changes within a gene.
        In order for this datasource to function properly, input mutations must already
          be annotated with gene and protein_change.
          (For example, GafDatasource does this).

        Additionally, this datasource assumes that a range of positions are going to be in one annotation
            separated by '_'

    TODO: allow required annotations from the config file.
    TODO: Range is still unsupported.  53_54
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, proteinPositionAnnotation="protein_change", gpColumnNames="gene,start_AA,end_AA"):
        # In this case, we want to initialize with the Datasource class
        super(GenericGenomicPositionDatasource, self).__init__(src_file, title=title, version=version)
        self.proteinPositionAnnotation = proteinPositionAnnotation
        self.proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)")
        index_mode = 'gene_protein_pos'
        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode, gpColumnNames)

    def _extractProteinPositions(self, otherTranscriptList):
        transcriptProteinPos = []
        # Parse the other_transcript annotation
        for ot in otherTranscriptList:
            # There was no protein change here (Intron, etc)
            otList = ot.split("p.")
            if len(otList) < 2:
                continue
            if len(otList) > 2:
                logging.getLogger(__name__).warn("More than one protein change detected in an other_transcript annotation.")
            proteinChange = otList[1]
            proteinPosition = self.proteinRegexp.match(proteinChange)
            if proteinPosition is not None:
                transcriptProteinPos.append(proteinPosition.group(1))

        return transcriptProteinPos

    def annotate_mutation(self, mutation):
        requiredAnnotations = ['gene', self.proteinPositionAnnotation]
        if all(field in mutation for field in requiredAnnotations):
            gene = mutation['gene']
            p = None
            transcriptProteinPos = []
            proteinChange = mutation[self.proteinPositionAnnotation].replace("p.", "")
            proteinPosition = self.proteinRegexp.match(proteinChange)
            if proteinPosition is not None:
                p = proteinPosition.group(1)

            if p is not None:
                records = get_binned_data(self.db_obj, gene, int(p), int(p))
                records = get_overlapping_records(records, int(p), int(p))

                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)
        else:
            missingList = []
            for r in requiredAnnotations:
                if r not in mutation:
                    missingList.append(r)
            raise MissingAnnotationException("Missing required annotation: " + str(missingList))

        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)

        return mutation