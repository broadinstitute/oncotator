# LICENSE_GOES_HERE

import logging
import os
from oncotator.datasources.Datasource import Datasource


class ReferenceDatasource(Datasource):
    """ Reference annotations.  A custom datasource initialized by genome flat files.

    Genome flat files are a simple format.  TODO: Finish this blurb.

    All annotations are with respect to the reference.

    Provides the following annotations:
        ref_context -- Small window into the reference at the variant.  The center position should be the same as Reference_Allele.  The total string should be of odd length and have a minimum length of 3.
            For example (SNV): Reference Allele is G, Chromosome is 1, Start_position and End_position are 120906037:  ref_context is CTTTTTTCGCGCAAAAATGCC  (string size is 21, in this case)


        gc_content --

        TODO: Finish the documentation.
    """

    def __init__(self, src_dir, title="Flat File Reference", version='hg19', windowSizeRef=10, windowSizeGCContent=100):
        """ Constructor
        src_dir is the parent directory of the chrXXX.txt files.
        """

        self.directoryName = src_dir
        self.windowSizeRef = int(windowSizeRef)
        self.windowSizeGC = windowSizeGCContent
        self.logger = logging.getLogger(__name__)
        super(ReferenceDatasource, self).__init__(src_dir, title=title, version=version)
        self._filePointers = dict()

    def annotate_mutation(self, m):

        # TODO: Error checking
        iStart = int(m['start'])
        iEnd = int(m['end'])
        m.createAnnotation('ref_context', self.getRange(m['chr'], iStart - (self.windowSizeRef+1), iEnd + (self.windowSizeRef-1)), annotationSource=self.title)

        # Populate gc content
        gcWindow = self.getRange(m['chr'], iStart - (self.windowSizeGC+1), iEnd + (self.windowSizeGC-1))
        gcCountInWindow = gcWindow.count("C") + gcWindow.count("G") + gcWindow.count("c") + gcWindow.count("g")

        if len(gcWindow) != 0:
            gc_content = "%0.3f" % (float(gcCountInWindow)/float(len(gcWindow)))
        else:
            gc_content = "0"
        m.createAnnotation('gc_content', gc_content, self.title)
        return m

    def _getFilePointer(self, fullChrFilename):
        """ Implements lazy loading of file pointers. """
        if fullChrFilename not in self._filePointers.keys():
            self._filePointers[fullChrFilename] = file(fullChrFilename, 'r')
        return self._filePointers[fullChrFilename]

    def getRange(self, chr, start, end):
        """ Returns string of the reference genome in the range requested. """

        chrFilename = self.convertMutationChrToFilename(chr)

        # TODO: We need a master lookup table.  Chip Stewart may have one for all reference builds.

        iEnd = int(end)
        iStart = max(int(start),0)
        if iEnd < iStart:
            return ""

        fullChrFilename = self.directoryName + "/" + chrFilename

        if not (os.path.exists(fullChrFilename)):
            self.logger.warn(fullChrFilename + " not found.  Please add it.")
            return ""
        else:
            fp = self._getFilePointer(fullChrFilename)

        fp.seek(iStart, 0) # Second parameter indicates "from the beginning of the file"
        result = fp.read(iEnd-iStart+1)
        return result

    # TODO: Low priority: Need cleanup of file pointers

    def convertMutationChrToFilename(self, chr):
        """ Convert the standard mutation chromosome convention to the convention used for the filenames.

        Examples (for hg19):
        GL000209.1 --> chr19_gl000209_random.txt
        X --> chrX.txt
        20 --> chr20.txt
        """
        result = "chr" + str(chr) + ".txt"

        # If we have a chr that starts with GL:
        #    1) Make all lowercase
        #    2) Strip away anything after '.'
        #    3) Find the file that contains the result of step 1 and 2.
        if chr.startswith("GL"):
            tmp = chr.lower()
            if tmp.find(".") != -1:
                tmp = tmp[0:tmp.find(".")]

            allFiles = os.listdir(self.directoryName)
            for f in allFiles:
                if f.find(tmp)<> -1:
                    result = f
                    break

        return result