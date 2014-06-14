# LICENSE_GOES_HERE


from oncotator.output.OutputRenderer import OutputRenderer
import logging
import csv


class SimpleBedOutputRenderer(OutputRenderer):
    """
    Responsible for conversion of MutationData to a minimal BED File.
    
    This class has very limited functionality compared to other output renderers.

    Any config file specified during initialization is ignored.
    
    IMPORTANT NOTES: 
    1)  This class does not handle any optional BED fields.  Functionality in this
        class is mostly meant for generating simple liftover files.
    2)  Comments are ignored and not rendered.
    3)  Annotations are discarded.
    4)  If chr, start, and end are not defined for a mutation, do not expect this
        class to behave well.
    """

    def __init__(self, filename, configFile="", otherOptions=None):
        """
        Constructor

        :param filename: output filename
        :param configFile: configuration filename
        """
        self._outputFilename = filename
        self.logger = logging.getLogger(__name__)
        self.logger.warn("BED Renderer initialized, please note that only the most basic BED files are produced.")

    def renderMutations(self, mutations, metadata=None, comments=None):
        """
        Generate a simple ssv file that only has three columns for chrom, start, and end.  Zero-indexing conversion is
        done here.

        :param mutations: generator of MutationData objects
        :param metadata:
        :param comments:
        :return:
        """

        self.logger.warn("BED Rendering started, please note that much BED-related functionality is not implemented.")
        f = file(self._outputFilename, 'w')
        headers = ['chr', 'start', 'end']
        writer = csv.DictWriter(f, headers, delimiter=' ', lineterminator='\n', extrasaction='ignore')
        row = dict()
        for m in mutations:
            row['chr'] = 'chr' + m.chr
            row['start'] = str((int(m.start) - 1))
            row['end'] = m.end
            writer.writerow(row)
        f.close()