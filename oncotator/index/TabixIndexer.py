# LICENSE_GOES_HERE


'''
Created on Jan 10, 2013

@author: lichtens
'''
import shutil
from oncotator.utils.GenericTsvReader import GenericTsvReader
import csv
from oncotator.utils.MutUtils import MutUtils
import pysam
from oncotator.utils.TsvFileSorter import TsvFileSorter
import os
import tempfile
import string
from oncotator.index.TabixIndexerFileMissingError import TabixIndexerFileMissingError


class TabixIndexer(object):
    """
    Static class (do not instantiate) that handles the creation of datasource indexes.

    Typically, the index is created once and used in the annotation process.

    In other words, index is called during datasource install, not datasource access.

    IMPORTANT: Currently, only tsv files are supported.

    """

    def __init__(self, params):
        """
        Constructor
        """
        raise NotImplementedError("This is a static class -- do not instantiate")
    
    @staticmethod
    def index(destDir, inputFilename, fileColumnNumList=None, preset=None):
        """
        Create a tabix index file for genomic position datasource tsv files.
        Prerequisites (for genomic position indexed):
            Input file has three columns that can be mapped to chromosome, start position, and end position without any modification.
                For example, ['hg19.oreganno.chrom', 'hg19.oreganno.chromStart', 'hg19.oreganno.chromEnd'] in oreganno.hg19.txt

        This will overwrite an existing index (since the force parameter is set to True in pysam.tabix_index() call).
        Also, in cases where the inputFilename doesn't end with a ".gz", the a compressed file will be created and indexed.

        :param destDir: destination directory
        :param ds_foldername: destination folder name
        :param fileColumnNumList: ordered list.  This list contains the corresponding entries (column numbers)
            in the tsv file. Typically, this would be [chr,start,end]  or [gene, startAA, endAA]
        :param inputFilename: tsv file input
        :param preset: if preset is provided, the column coordinates are taken from a preset. Valid values for preset
        are "gff", "bed", "sam", "vcf", "psltbl", and "pileup".
        """
        fileColumnNumList = [] if fileColumnNumList is None else fileColumnNumList
        inputFilename = os.path.abspath(inputFilename)
        fileDir = os.path.dirname(inputFilename)
        fileName, fileExtension = os.path.splitext(os.path.basename(inputFilename))

        if fileExtension in (".gz",):
            # Ensure .gz.tbi file is there as well
            inputIndexFilename = os.path.join(fileDir, string.join([inputFilename, "tbi"], "."))
            if not os.path.exists(inputIndexFilename):
                msg = "Missing tabix index file %s." % inputIndexFilename
                raise TabixIndexerFileMissingError(msg)

            outputFilename = os.path.join(destDir, string.join([fileName, "gz"], "."))
            shutil.copyfile(inputFilename, outputFilename)

            outputIndexFilename = os.path.join(destDir, string.join([fileName, "gz", "tbi"], "."))
            shutil.copyfile(inputIndexFilename, outputIndexFilename)

            return outputFilename

        outputFilename = os.path.join(destDir, string.join([fileName, ".tabix_indexed", fileExtension], ""))
        # Copy the input file to output file.
        shutil.copyfile(inputFilename, outputFilename)

        # Load the file into a tsvReader.
        if preset in ("gff", "bed", "sam", "vcf", "psltbl", "pileup"):
            tabix_index = pysam.tabix_index(filename=outputFilename, force=True, preset=preset)
        else:
            tabix_index = pysam.tabix_index(filename=outputFilename, force=True, seq_col=fileColumnNumList[0],
                                            start_col=fileColumnNumList[1], end_col=fileColumnNumList[2])

        return tabix_index


    @staticmethod
    def indexGeneProteinPosition(geneColumn, proteinInfoColumn, inputFilename, outputFilename):
        """
        Creates an intermediate temporary file that includes two additional columns, startAA and endAA,
        sorts the file, writes thee sorted file to outputFilename, and then indexes the sorted file.

        :param geneColumn: name of the gene column in the inputFilename
        :param proteinInfoColumn: name of the protein change or position column. Can be of formats: p.K128_R130del
        (position 128 through 130) For more examples, see MutUtilsTest.testProteinChange()
        :param inputFilename: input tsv filename
        :param outputFilename: output filename
        """
        startAACol = "startAA"
        endAACol = "endAA"

        # Create intermediate file.  Do not use '#' for comments, since header can start with '#'
        tsvReader = GenericTsvReader(inputFilename, commentPrepend=";")

        # These are the outputHeaders for the intermediate file.
        headers = tsvReader.getFieldNames()

        if startAACol not in headers:
            headers += [startAACol]
        if endAACol not in headers:
            headers += [endAACol]

        # Write to the intermediate temporary file.
        # This file is created in the current working directory."
        temp = tempfile.NamedTemporaryFile()
        csvfile = file(temp.name, 'w')

        # Initialize the intermediate file's header.
        tsvWriter = csv.DictWriter(csvfile, headers, delimiter='\t', lineterminator='\n')
        # If the headers have a leading '#', get rid of it.
        for i in range(0, len(headers)):
            header = headers[i]
            if header.startswith("#"):
                headers[i] = header.replace("#", "")
        tsvWriter.writeheader()

        # Get indices of relevant columns.
        gene_i = headers.index(geneColumn)
        startAA_i = headers.index(startAACol)
        endAA_i = headers.index(endAACol)

        # Write each line of the intermediate file.
        for row in tsvReader:
            protein = row[proteinInfoColumn]
            if protein is None or not protein.strip():
                continue
            [startAA, endAA] = MutUtils.extractProteinPosition(protein)
            if not startAA.strip() or not endAA.strip():
                continue
            row[startAACol] = startAA
            row[endAACol] = endAA
            tsvWriter.writerow(row)
        csvfile.flush()
        csvfile.close()

        # Sort the intermediate tsv file.
        tsvSorter = TsvFileSorter(temp.name)
        func = lambda val: ((val["Gene name"]).lower(), int(val["startAA"]), int(val["endAA"]))

        # Use the whole file path name.
        outputFilename = os.path.abspath(outputFilename)
        tsvSorter.sortFile(outputFilename, func)

        return TabixIndexer.index(destDir=os.path.dirname(os.path.abspath(outputFilename)),
                                  inputFilename=outputFilename, fileColumnNumList=[gene_i, startAA_i, endAA_i])
