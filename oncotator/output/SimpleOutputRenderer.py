# LICENSE_GOES_HERE


from OutputRenderer import OutputRenderer
import csv
import logging
from oncotator.utils.MutUtils import MutUtils
import itertools
import copy


class SimpleOutputRenderer(OutputRenderer):
    """
    The SimpleOutputRenderer renders a basic tsv file from the given mutations.  All annotations are included with real names as column headers.
    
    Header is determined by the first mutation given.
    
    No attention is paid to order of the headers.    
    
    Any headers starting with "_" are removed and not rendered.

    There is no config file needed for initialization of this class.  If specified, it is ignored.
    """

    def __init__(self, filename, configFile="", otherOptions=None):
        """
        Constructor

        Config file parameter is ignored.
        :param filename:
        :param configFile:
        """
        self._filename = filename
        self._logger = logging.getLogger(__name__)
        self._lineterminator = "\n"
        self._delimiter = "\t"

    def _determineHeaders(self, mut, metadata):
        headers = MutUtils.getAllAttributeNames(mut) if not None else []
        if len(headers) == 0:
            headers = metadata.keys()

        # Remove headers that start with "_"
        for header in headers:
            if header.startswith("_"):
                headers.remove(header)

        return headers

    def renderMutations(self, mutations, metadata=None, comments=None):
        """ Generate a simple tsv file based on the incoming mutations.
        
        Assumes that all mutations have the same annotations, even if some are not populated.
        
        Any annotation name starting with "_" will be ignored.
        
        Returns a file name. """
        
        self._logger.info("Simple rendering output file: " + self._filename)
        self._logger.info("Render starting...")

        comments = [] if comments is None else comments

        ctr = 0
        fptr = file(self._filename, 'w')
        if len(comments) != 0:
            fptr.write('#' + "\n# ".join(comments) + "\n")
            
        mut = None
        for mutation in mutations:
            mut = copy.deepcopy(mutation)
            lst = [mutations, (mut for mut in [mut])]
            mutations = itertools.chain(*lst)

        headers = self._determineHeaders(mut, metadata)

        writer = csv.DictWriter(fptr, headers, delimiter=self._delimiter, lineterminator=self._lineterminator,
                                extrasaction="ignore")
        writer.writeheader()

        for m in mutations:
            writer.writerow(m)
            # Update mutation count and log every 1000 mutations
            ctr += 1
            if (ctr % 1000) == 0:
                self._logger.info("Rendered " + str(ctr) + " mutations.")

        fptr.close()
        
        self._logger.info("Rendered " + str(ctr) + " mutations into " + self._filename + ".")
        return self._filename