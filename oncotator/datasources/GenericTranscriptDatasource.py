# LICENSE_GOES_HERE

from oncotator.datasources.GenericGeneDatasource import GenericGeneDatasource


class GenericTranscriptDatasource(GenericGeneDatasource):
    """ Used for generic TSV that is indexed by transcript ID. """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='transcript_id'):
        super(GenericTranscriptDatasource,self).__init__(src_file, title, version, use_binary, geneColumnName)

    def annotate_mutation(self, mutation):
        return super(GenericTranscriptDatasource,self).annotate_mutation(mutation,'transcript_id')