import logging
from shove.core import Shove
from oncotator.Annotation import Annotation
from oncotator.datasources.Datasource import Datasource
import leveldb


class LevelDbDatasource(Datasource):
    """Datasource backed by LEvelDB (no Shove usage).

    Initialization is done using a file location.

    For example:
    /absolute/path/to/a/folder/
    relative/path/to/a/folder/

    Index columns must match exactly in order to receive a value from the ShoveDatasource

    Note:  This was pulled out of Shove due to CPU constraints.
        Shove was executing a lot of extraneous sync statements which are unnecessary, since this datasource is read-only

    IMPORTANT:  This datasource only supports SNPs.
    IMPORTANT:  This datasource only supports exact matches to chromosome, start, end, ref, and alt.
        Typically, start = end.
    IMPORTANT:  This datasource only supports hg19 chr1-22,X,Y

    """

    def __init__(self, src_file, title, version, index_cols, annotation_columns):
        """

        Example:  LevelDbDatasource("/path/to/some/file", "dbNSFP", "2.5", ["#chr", "pos(1-coor)", "pos(1-coor)", "ref", "alt"], [...])

        Example config file:
        [general]
        version = TEST
        title = dbNSFP_chr1_cut
        type = leveldb
        src_file = dbNSFP2.4_variant.tabix_indexed.chr1_cut.gz
        annotation_column_names = aaref,aaalt,hg18_pos(1-coor),genename,Uniprot_acc,Uniprot_id,Uniprot_aapos,Interpro_domain,cds_strand,refcodon,SLR_test_statistic,codonpos,fold-degenerate,Ancestral_allele,Ensembl_geneid,Ensembl_transcriptid,aapos,aapos_SIFT,aapos_FATHMM,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,RadialSVM_score,RadialSVM_rankscore,RadialSVM_pred,LR_score,LR_rankscore,LR_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF
        index_column_names = #chr,pos(1-coor),pos(1-coor),ref,alt


        :param src_file: LevelDB dir to initialize in this instance.
        :param title:
        :param version:
        :param index_cols:
        :param annotation_columns: list of columns, in order, to annotate with.
        """
        super(LevelDbDatasource, self).__init__(src_file, title=title, version=version)

        # Initialize a level db datasource w/ 100MB of memory cache
        self._db_store = leveldb.LevelDB(src_file, block_cache_size=(100 * (2 << 20)), create_if_missing=False)
        self._annotation_columns = annotation_columns
        self._blank_annotations = {k:Annotation("", datasourceName=self.title, dataType="String", description="") for k in self._annotation_columns}
        self._preload_start = -1
        self._preload_end = -1

    def _alter_preload(self, h):
        if h >=  self._preload_start and h <= self._preload_end:
            pass
        else:
            self._db_store.RangeIter(key_from = h, key_to = str(int(h) + (1000000 << 6)), include_value = True, verify_checksums = False, fill_cache = True)


    def annotate_mutation(self, mutation):
        """ Mutations are annotated only with exact matches of chr, start, end, ref, and alt.
        """
        if mutation.ref_allele == "" or mutation.ref_allele == "-" or mutation.alt_allele == "" or mutation.alt_allele == "-" or len(mutation.ref_allele) > 1 or len(mutation.alt_allele) > 1:
            # TODO: Above statement also needs to check for chromosomes not in [1-22, X, Y]

            # We do not have anything to do here, since it is not a SNP.
            # Populate all annotations with blank values
            mutation.addAnnotations(self._blank_annotations)

        else:
            mutation = self._perform_annotate_mutation(mutation)

        return mutation

    def _perform_annotate_mutation(self, mutation):
        # create hash for this mutation
        h = LevelDbDatasource.generate_hash(mutation.chr, mutation.start, mutation.end, mutation.ref_allele, mutation.alt_allele)

        # extract value for this hash from the db
        annotations_list = []
        try:
            self._alter_preload(h)
            annotations_list = self._db_store.Get(h).split(",")

            # # Annotate
            # for i,col in enumerate(self._annotation_columns):
                # if len(annotations_list) <= i:
                #     # TODO: Throw exception here instead?
                #     logging.getLogger(__name__).error("Disconcordant length of annotation columns between datasource config file and the actual data.")
                #     mutation.createAnnotation(col, "", self.title)
                # else:
                #     mutation.createAnnotation(col, annotations_list[i], self.title)


            [mutation.createAnnotation(col, annotations_list[i], self.title)  for i,col in enumerate(self._annotation_columns)]
        except KeyError:
            # do nothing
            mutation.addAnnotations(self._blank_annotations)



        return mutation

    @staticmethod
    def generate_hash(m):
        return LevelDbDatasource.generate_hash(m.chr, m.start, m.end, m.ref_allele, m.alt_allele)

    @staticmethod
    def _generate_int_from_base(b):
        if b == "A": return 0
        if b == "C": return 1
        if b == "G": return 2
        if b == "T": return 3
        else: return 4

    @staticmethod
    def generate_hash(chrom, start, end, ref, alt):

        #TODO: Non-hg19 builds
        #TODO: support for chomosomes other than chr1-22,X,Y
        #TODO: Convert to genomic space using actual length of contigs, not just 300M

        if chrom == "X":
            raw_chrom = 23
        elif chrom == "Y":
            raw_chrom = 24
        else:
            raw_chrom = int(chrom)
        chrom_pos_offset = ((raw_chrom-1) * 300000000) + int(start)
        chrom_pos_offset = chrom_pos_offset << 6

        b_ref = LevelDbDatasource._generate_int_from_base(ref)
        b_alt = LevelDbDatasource._generate_int_from_base(alt)
        chrom_pos_offset += (b_ref << 3)
        chrom_pos_offset += b_alt

        return str(chrom_pos_offset)

    def get_stats(self):
        """ get stats for the level db database
        """
        return self._db_store.GetStats()