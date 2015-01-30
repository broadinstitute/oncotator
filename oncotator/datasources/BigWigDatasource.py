import numpy as np

from oncotator.datasources.Datasource import Datasource

class BigWigDatasource(Datasource):
    """
    A datasource derived from a BigWig file.  For variants spanning a genomic range (i.e. non SNVs),
    the median of values from the BigWig are returned.
    """
    def __init__(self, src_file, title='', version=None):
        # only necessary to import ngslib if instance of BigWigDatasource is created
        # This should not run on OS X machines
        from ngslib import BigWigFile

        super(BigWigDatasource, self).__init__(src_file, title=title, version=version)

        self.output_headers = [title + '_score']
        self.bigwig_fh = BigWigFile(src_file)
        self.has_chr = True if self.bigwig_fh.chroms[0].startswith('chr') else False

    def annotate_mutation(self, mutation):
        if self.has_chr and not mutation.chr.startswith('chr'):
            chrn = 'chr' + mutation.chr
        else:
            chrn = mutation.chr

        variant_start, variant_end = int(mutation.start) - 1, int(mutation.end) #start - 1 because bigwig format is zero-based coords

        scores = [r[2] for r in self.bigwig_fh.fetch(chrom=chrn, start=variant_start, stop=variant_end)]

        if not scores:
            final_score = None
        elif len(scores) == 1:
            final_score = scores[0]
        else:
            final_score = np.median(scores)

        mutation.createAnnotation(self.output_headers[0], final_score, annotationSource=self.title)
        return mutation

    def close(self):
        self.bigwig_fh.close()