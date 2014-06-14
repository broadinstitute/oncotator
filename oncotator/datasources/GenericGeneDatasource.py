# LICENSE_GOES_HERE

from oncotator.datasources.GenericGeneDataSourceException import GenericGeneDataSourceException
from oncotator.utils.db import get_db_data
from oncotator.datasources.Datasource import Datasource


class GenericGeneDatasource(Datasource):
    """
    A datasource derived from a generic TSV file in which the first column is a HUGO gene
    symbol.  First header value must be specified in the constructor (default: 'gene').  All other columns will be used for
    annotation.
    TODO: This is no longer true.  Actually, you can specify any column as the gene column in the config file.  Update documentation.
    use_binary
        if True, existing indexed binary will be used or created for future use.

    """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='gene'):
        super(GenericGeneDatasource, self).__init__(src_file, title=title, version=version)

        index_mode = 'gene'

        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode,indexColumnNames=geneColumnName)

    def annotate_mutation(self, mutation, index_field='gene'):

        if index_field not in mutation:
            raise GenericGeneDataSourceException("Index field (" + index_field + ") not found.  Remember that datasources must be ordered.  Please put a datasource that provides the 'gene' annotation in front of this datasource.")

        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))

        gene = mutation[index_field]
        if gene in self.db_obj:
            annotations = self.db_obj[gene]
            for k in annotations.keys():
                mutation.createAnnotation(k, self.db_obj[gene][k], annotationSource=self.title)
        else:
            for h in self.output_headers:
                mutation.createAnnotation(h, '', annotationSource=self.title)

        return mutation