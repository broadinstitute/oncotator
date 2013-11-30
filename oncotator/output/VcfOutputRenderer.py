"""
# By downloading the PROGRAM you agree to the following terms of use:
#
# BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
# FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
#
# This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
# WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
# WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
# NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
#
# 1. DEFINITIONS
# 1.1	"PROGRAM" shall mean copyright in the object code and source code known as Oncotator and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.
#
# 2. LICENSE
# 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
# LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
# The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
# 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
# 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
#
# 3. OWNERSHIP OF INTELLECTUAL PROPERTY
# LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
#
# Copyright 2012 Broad Institute, Inc.
# Notice of attribution:  The Oncotator program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc.
#
# LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
#
# 4. INDEMNIFICATION
# LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
#
# 5. NO REPRESENTATIONS OR WARRANTIES
# THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
# IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
#
# 6. ASSIGNMENT
# This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
#
# 7. MISCELLANEOUS
# 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
# 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
# 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
# 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
# 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
# 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
# 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
#"""

import collections
import logging
import os
import tempfile
import vcf

from oncotator.utils.MutUtils import MutUtils
from OutputRenderer import OutputRenderer
import oncotator.utils.ConfigUtils
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.output.RecordBuilder import RecordBuilder
from oncotator.output.OutputDataManager import OutputDataManager
from oncotator.config_tables.ConfigTableCreatorFactory import ConfigTableCreatorFactory


class VcfOutputRenderer(OutputRenderer):
    """
    The VcfOutputRenderer renders a vcf file from the given mutations.  All annotations are included with real names
    as column headers.

    Header is determined by the first mutation given.

    No attention is paid to order of the headers.
    """
    _vcfAnnotation = collections.namedtuple(typename="Annotation", field_names=["field", "ID", "num", "type"])

    def __init__(self, filename, datasources=None, configFile='vcf.out.config'):
        """


        :param filename:
        :param datasources:
        :param configFile: output config file
        """
        self._filename = filename
        self.configFilename = configFile
        self.logger = logging.getLogger(__name__)
        self._datasources = [] if datasources is None else datasources
        self.config = oncotator.utils.ConfigUtils.ConfigUtils.createConfigParser(configFile, ignoreCase=False)
        self.chromHashCodeTable = None  # maps every chromosome in the mutations to a sortable integer
        self.configTableBuilder = ConfigTableCreatorFactory.getConfigTableCreatorInstance("output_vcf")
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.sampleNames = []  # all sample names in the mutations
        self.chroms = []  # all chromosomes in the mutations

    def renderMutations(self, mutations, metadata=None, comments=None):
        """ Generate a simple tsv file based on the incoming mutations.
        Assumes that all mutations have the same annotations, even if some are not populated.
        :param mutations:
        :param metadata:
        :param comments:
        Returns a file name. """

        metadata = [] if metadata is None else metadata
        comments = [] if comments is None else comments

        self.logger.info("Rendering VCF output file: " + self._filename)

        # Initialize config table
        self.configTable = self.configTableBuilder.getConfigTable(configFilename=self.configFilename)

        path = os.getcwd()

        dataManager = OutputDataManager(self.configTable, mutations, comments, metadata, path)
        self.sampleNames = dataManager.getSampleNames()

        # Write the header
        tempTemplateFile = tempfile.NamedTemporaryFile(dir=path, delete=False)
        filePointer = file(tempTemplateFile.name, 'w')
        header = dataManager.getHeader()
        filePointer.write(header)
        filePointer.close()

        sortedTempTsvFileName = dataManager.getSortedTsvFilename(path)

        self.logger.info("Render starting...")
        self._renderSortedTsv(tempTemplateFile.name, self._filename, sortedTempTsvFileName, self.sampleNames,
                              dataManager)

        # Remove template filename
        os.remove(tempTemplateFile.name)

        # Remove sorted tsv filename
        os.remove(sortedTempTsvFileName)

        self.logger.info("Rendered all mutations.")

        return self._filename

    def _isNewVcfRecordNeeded(self, curChrom, prevChrom, curPos, prevPos):
        """

        :param curChrom:
        :param prevChrom:
        :param curPos:
        :param prevPos:
        :return:
        """
        isNew = False
        if curChrom != prevChrom:
            isNew = True
        if curPos != prevPos:
            isNew = True
        return isNew

    def _renderSortedTsv(self, templateFilename, vcfFilename, tsvFilename, sampleNames, dataManager):
        """

        :param vcfFilename:
        :param tsvFilename:
        :param sampleNames:
        :param dataManager:
        """
        tempVcfReader = vcf.Reader(filename=templateFilename, strict_whitespace=True)
        pointer = file(vcfFilename, "w")
        vcfWriter = vcf.Writer(pointer, tempVcfReader, self.lineterminator)
        tsvReader = GenericTsvReader(tsvFilename, delimiter=self.delimiter)
        index = 0
        nrecords = 1000
        chrom = None
        pos = None
        recordBuilder = None
        for m in tsvReader:
            isNewRecord = self._isNewVcfRecordNeeded(chrom, m["chr"], pos, m["start"])
            if isNewRecord:
                if recordBuilder is not None:
                    record = recordBuilder.createRecord()
                    vcfWriter.write_record(record)
                    index += 1
                    if index % nrecords == 0:
                        self.logger.info("Rendered " + str(index) + " vcf records.")
                        vcfWriter.flush()

                chrom = m["chr"]
                if chrom.startswith("GL"):
                    chrom = "<" + chrom + ">"
                pos = m["start"]
                refAllele = m["ref_allele"]

                recordBuilder = RecordBuilder(chrom, int(pos), refAllele, sampleNames)

            recordBuilder = self._parseRecordBuilder(m, recordBuilder, dataManager)

        if recordBuilder is not None:
            record = recordBuilder.createRecord()
            vcfWriter.write_record(record)

        vcfWriter.close()
        self.logger.info("Rendered all " + str(index) + " vcf records.")

    def _parseRecordBuilder(self, m, recordBuilder, dataManager):
        """
        Parse the input mutation object.
        First, this method

        :param m: mutation object
        :param recordBuilder:
        :param dataManager:
        :return:
        """

        idAnnotationNames = dataManager.getAnnotationNames("ID")
        qualAnnotationNames = dataManager.getAnnotationNames("QUAL")
        filterAnnotationNames = dataManager.getAnnotationNames("FILTER")
        infoAnnotationNames = dataManager.getAnnotationNames("INFO")
        formatAnnotationNames = dataManager.getAnnotationNames("FORMAT")
        sampleNameAnnotationNames = dataManager.getAnnotationNames("SAMPLE_NAME")

        if len(sampleNameAnnotationNames) != 0:
            sampleNameAnnotationName = sampleNameAnnotationNames[0]
        else:
            sampleNameAnnotationName = MutUtils.SAMPLE_NAME_ANNOTATION_NAME
        sampleName = m.get(sampleNameAnnotationName, None)

        altAllele = m["alt_allele"]

        recordBuilder.addAlt(altAllele)

        for name in idAnnotationNames:
            val = m.get(name, "")
            recordBuilder.addID(val)

        if len(qualAnnotationNames) == 0:
            qual = "qual"
        else:
            qual = qualAnnotationNames[0]
        recordBuilder.addQual(m.get(qual, "."))

        for name in filterAnnotationNames:
            ID = dataManager.getFieldID(name)
            val = m.get(name, "")
            recordBuilder.addFilter(ID, val)

        for name in infoAnnotationNames:
            annotation = dataManager.getOutputAnnotation(name)
            ID = annotation.getID()
            num = annotation.getNumber()
            dataType = annotation.getDataType()
            isSplit = annotation.isSplit()
            val = m.get(name, "")
            recordBuilder.addInfo(sampleName, ID, num, dataType, val, isSplit)

        for name in formatAnnotationNames:
            annotation = dataManager.getOutputAnnotation(name)
            ID = annotation.getID()
            num = annotation.getNumber()
            dataType = annotation.getDataType()
            isSplit = annotation.isSplit()
            val = m.get(name, "")
            if num == 0 or dataType == "Flag":
                msg = "%s is of data type Flag. Only Integer, Float, Character, and String data types are permissible" \
                      " in the Format field." % name
                logging.warn(msg)
            else:
                recordBuilder.addFormat(sampleName, ID, num, dataType, val, isSplit)

        return recordBuilder