from Bio import Seq
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils


class Transcript(object):
    """Simple class to hold transcript information in a standard way across Transcript Providing datasources

    transcript_id -- unique id from the transcript (e.g. uc010ajp.1 for GAF) (string)
    exons are a list of (start, end) in genomic coordinates (integers)

    """
    def __init__(self, transcript_id, gene, contig, gene_id="", seq="", strand="+", start_codon=None, stop_codon=None, gene_type=""):
        self._transcript_id = transcript_id
        self._exons = []
        self._cds = []
        self._gene = gene
        self._gene_id = gene_id
        self._seq = seq
        self._contig = contig
        self._strand = strand
        self._start_codon = start_codon
        self._stop_codon = stop_codon
        self._other_attributes = {}
        self._gene_type = gene_type
        self._protein_seq = None

    def add_exon(self, start, end, exon_number):
        self._exons.append((start, end, exon_number))

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
        return max(exon_ends)

    def get_contig(self):
        return self._contig

    def set_contig(self, value):
        self._contig = value

    def get_strand(self):
        return self._strand

    def set_strand(self, value):
        self._strand = value

    def set_start_codon(self, start, end):
        self._start_codon = (start, end)

    def set_stop_codon(self, start, end):
        self._stop_codon = (start, end)

    def get_stop_codon(self):
        return self._stop_codon

    def get_start_codon(self):
        return self._start_codon

    def add_other_attribute(self, key, value):
        self._other_attributes[key] = value

    def get_other_attributes(self):
        return self._other_attributes

    def get_gene_type(self):
        return self._gene_type

    def set_gene_type(self, value):
        self._gene_type = value

    # TODO: Reduce code duplication
    def determine_transcript_start(self):
        """Includes UTR but not padding
        Returns the start location in genomic space.  This will be the highest position if strand is negative."""
        all_locations_start = [s[0] for s in self._exons]
        all_locations_end = [s[1] for s in self._exons]
        all_locations = []
        all_locations.extend(all_locations_start)
        all_locations.extend(all_locations_end)
        if self._strand == "-":
            return max(all_locations)
        else:
            return min(all_locations)

    def determine_transcript_stop(self):
        """Includes UTR but not padding
        Returns the stop location in genomic space.  This will be the highest position if strand is positive."""
        all_locations_start = [s[0] for s in self._exons]
        all_locations_end = [s[1] for s in self._exons]
        all_locations = []
        all_locations.extend(all_locations_start)
        all_locations.extend(all_locations_end)
        if self._strand == "-":
            return min(all_locations)
        else:
            return max(all_locations)

    def determine_cds_start(self):
        all_locations_start = [s[0] for s in self._cds]
        all_locations_end = [s[1] for s in self._cds]
        all_locations = []
        all_locations.extend(all_locations_start)
        all_locations.extend(all_locations_end)
        if len(all_locations) == 0:
            return -1
        if self._strand == "-":
            return max(all_locations)
        else:
            return min(all_locations)

    def determine_cds_stop(self):
        all_locations_start = [s[0] for s in self._cds]
        all_locations_end = [s[1] for s in self._cds]
        all_locations = []
        all_locations.extend(all_locations_start)
        all_locations.extend(all_locations_end)
        if len(all_locations) == 0:
            return -1
        if self._strand == "-":
            return min(all_locations)
        else:
            return max(all_locations)

    def get_protein_seq(self):
        return self._protein_seq

    def set_protein_seq(self, value):
        self._protein_seq = value

    def determine_cds_footprint(self):
        """ Returns the cds in genomic space.  Note that strand is ignored, so the first return value is always lower.
        :param tx:
        :return: cds_start, cds_stop in genomic coordinates.
        """
        s = cds_start = self.determine_cds_start()
        e = cds_stop = self.determine_cds_stop()
        if cds_stop < cds_start:
            s = cds_stop
            e = cds_start
        return s, e

