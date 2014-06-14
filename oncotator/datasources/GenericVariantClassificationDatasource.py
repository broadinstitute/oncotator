# LICENSE_GOES_HERE

from oncotator.datasources.GenericGeneDatasource import GenericGeneDatasource


class GenericVariantClassificationDatasource(GenericGeneDatasource):
    """ Used for generic TSV that is indexed by variant classification. """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='variant_classification'):
        super(GenericVariantClassificationDatasource,self).__init__(src_file, title, version, use_binary, geneColumnName)

    def annotate_mutation(self, mutation):
        return super(GenericVariantClassificationDatasource,self).annotate_mutation(mutation,'variant_classification')