
class Transcript(object):
    """Simple class to hold transcript information in a standard way across Transcript Providing datasources

    transcript_id -- unique id from the transcript (e.g. uc010ajp.1 for GAF) (string)
    exons are a list of (start, end) in genomic coordinates (integers)

    """
    def __init__(self, transcript_id, gene, contig, gene_id = "", exons = [], cds = [], seq=""):
        self._transcript_id = transcript_id
        self._exons = exons
        self._cds = cds
        self._gene = gene
        self._gene_id = gene_id
        self._seq = seq
        self._contig = contig

    def add_exon(self, start, end):
        self._exons.append((start, end))

    def add_cds(self, start, end):
        self._cds.append((start, end))

    def set_transcript_id(self, id):
        self._transcript_id = id

    def set_gene(self, gene):
        #should be self._gene = gene ?
        self._transcript_id = id

    def get_transcript_id(self):
        return self._transcript_id

    def get_gene(self):
        return self._gene

    def get_seq(self):
        return self._seq

    def set_seq(self, seq):
        self._seq = seq

    def get_exons(self):
        return self._exons

    def get_cds(self):
        return self._cds

    def get_start(self):
        exon_starts = [exon[0] for exon in self._exons]
        return min(exon_starts)

    def get_end(self):
        exon_ends = [exon[1] for exon in self._exons]
        return min(exon_ends)

    def get_contig(self):
        return self._contig

    def set_contig(self, value):
        self._contig = value