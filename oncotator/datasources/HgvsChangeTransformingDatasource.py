from oncotator.datasources.ChangeTransformingDatasource import ChangeTransformingDatasource
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils

class HgvsChangeTransformingDatasource(ChangeTransformingDatasource):

    def __init__(self, src_file, title='', version=None):
        super(ChangeTransformingDatasource, self).__init__(src_file, title=title, version=version)
        self.output_headers = ['HGVS_genomic_change', 'HGVS_coding_DNA_change', 'HGVS_protein_change']

        gencode_ds = EnsemblTranscriptDatasource(src_file)

        self.vcer = VariantClassifier()
        
    def annotate_mutation(self, mutation):
        # TODO
        # TODO: When annotating, you should probably give the annotation a name like:  self.title + "_" + output_annotation_name
        # You can assume that the annotations from the transript datasource have already been populated (such as protein_change and genome_change)
        if mutation['start'] == '80529551' and mutation['chr'] == '2':
            from IPython import embed
            embed()

        for output_field in self.output_headers:
        	mutation.createAnnotation(output_field, 'TEST', self.title)

        return mutation

### code for annotation intronic mutation on plus strand

#gencode_ds = EnsemblTranscriptDatasource('/Users/aramos/Desktop/ot_hgvs/oncotator_v1_ds_Feb2014/gencode_out2/hg19/gencode.v18.annotation.gtf')
#res = gencode_ds.transcript_db['ENST00000402739.4']
#self.vcer = VariantClassifier()
#
#nearest_exon = vcer._determine_closest_exon(res, int(mutation['start']), int(mutation['end']))
#
#dist_to_exon = vcer._get_splice_site_coordinates(res, int(mutation['start']), int(mutation['end']), nearest_exon)
#
#tx_exons = res.get_exons()
#
#cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_exons[6][0], tx_exons[6][0], res)
#cds_position_of_nearest_exon = cds_position_of_nearest_exon[0] + 1 #why add 1?
#dist_to_exonv= dist_to_exon - 1 #why substract 1?