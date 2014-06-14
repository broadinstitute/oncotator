# LICENSE_GOES_HERE

from GenericTsvDatasourceCreator import GenericTsvDatasourceCreator
from TabixIndexedTsvDatasourceCreator import TabixIndexedTsvDatasourceCreator
from TabixIndexedVcfDatasourceCreator import TabixIndexedVcfDatasourceCreator


class DatasourceBuilderFactory():

    @staticmethod
    def getDatasourceCreatorInstance(dsType):
        if dsType in ("gp_tsv", "gene_tsv", "transcript_tsv", "gpp_tsv"):
            return GenericTsvDatasourceCreator()
        elif dsType in ("indexed_vcf",):
            return TabixIndexedVcfDatasourceCreator()
        elif dsType in ("indexed_tsv",):
            return TabixIndexedTsvDatasourceCreator()