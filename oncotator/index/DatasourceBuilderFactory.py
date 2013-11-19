from oncotator.index.GenericTsvDatasourceBuilder import GenericTsvDatasourceBuilder
from oncotator.index.TabixIndexedTsvDatasourceBuilder import TabixIndexedTsvDatasourceBuilder
from oncotator.index.TabixIndexedVcfDatasourceBuilder import TabixIndexedVcfDatasourceBuilder


class DatasourceBuilderFactory():

    @staticmethod
    def getDatasourceBuilderInstance(dsType):
        if dsType in ("gp_tsv", "gene_tsv", "transcript_tsv", "gpp_tsv"):
            return GenericTsvDatasourceBuilder()
        elif dsType in ("indexed_vcf",):
            return TabixIndexedVcfDatasourceBuilder()
        elif dsType in ("indexed_tsv",):
            return TabixIndexedTsvDatasourceBuilder()