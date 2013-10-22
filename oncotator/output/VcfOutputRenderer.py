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

from OutputRenderer import OutputRenderer
from oncotator.utils.ConfigUtils import ConfigUtils
import csv
import collections
import logging
import os
import tempfile
import itertools
import string
import vcf
from oncotator.input.ConfigInputIncompleteException import ConfigInputIncompleteException
from oncotator.utils.TsvFileSorter import TsvFileSorter
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.output.RecordFactory import RecordFactory
from oncotator.output.OutputDataManager import OutputDataManager
from oncotator.utils.ConfigTable import ConfigTable
from oncotator.utils.MutUtils import MutUtils


class VcfOutputRenderer(OutputRenderer):
    """
    The VcfOutputRenderer renders a vcf file from the given mutations.  All annotations are included with real names
    as column headers.

    Header is determined by the first mutation given.

    No attention is paid to order of the headers.
    """
    _vcfAnnotation = collections.namedtuple(typename="Annotation", field_names=["field", "ID", "num", "type"])

    def __init__(self, filename, datasources=[], configFile='vcf.out.config'):
        """
        Constructor
        :param filename:
        :param datasources:
        :param configFile:
        """
        self._filename = filename
        self.logger = logging.getLogger(__name__)
        self._datasources = datasources
        self.config = ConfigUtils.createConfigParser(configFile, ignoreCase=False)
        self.chromHashCodeTable = None  # maps every chromosome in the mutations to a sortable integer
        self.configTable = ConfigTable()
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.sampleNames = []  # all sample names in the mutations
        self.chroms = []  # all chromosomes in the mutations
        self.reservedAnnotationNames = ['chr', 'start', 'ref_allele', 'alt_allele', 'end']

    def _writeMuts2Tsv(self, filename, fieldnames, muts):
        """

        :param filename:
        :param fieldnames:
        :param muts:
        """
        sampleNames = set()
        chroms = set()

        with open(filename, 'w') as fptr:
            writer = csv.DictWriter(fptr, fieldnames, extrasaction='ignore', delimiter=self.delimiter,
                                    lineterminator=self.lineterminator)
            writer.writeheader()
            ctr = 0
            for mut in muts:
                if "sampleName" in mut:
                    sampleName = mut['sampleName']
                    sampleNames.add(sampleName)

                # Parse chromosome
                chroms.add(mut.chr)

                if mut.ref_allele == "_" or mut.alt_allele == "-":
                    ref_allele, alt_allele, updated_start = MutUtils.retrievePrecedingBase(mut)
                    mut.start = updated_start
                    mut.ref_allele = ref_allele
                    mut.alt_allele = alt_allele

                writer.writerow(mut)

                ctr += 1
                if (ctr % 1000) == 0:
                    self.logger.info("Wrote " + str(ctr) + " mutations to tsv.")

        if len(sampleNames) > 0:
            self.sampleNames = list(sampleNames)
            self.sampleNames.sort()

        if len(chroms) > 0:
            self.chroms = list(chroms)

    def _getFieldnames(self, mut, md):
        """

        :param mut:
        :param md:
        :return:
        """
        fieldnames = self.reservedAnnotationNames
        if mut is not None:
            fieldnames = set(fieldnames).union(md.keys())
            fieldnames = fieldnames.union(mut.keys())
            fieldnames = fieldnames.difference(["end", "sampleName", "build"])
            if "sampleName" in mut:
                fieldnames = fieldnames.union(["sampleName"])
        return list(fieldnames)

    def _validateOutputConfigFile(self):
        """


        :raise:
        """
        sections = ["INFO", "FORMAT", "OTHER", "SPLIT_TAGS", "NOT_SPLIT_TAGS", "INFO_DESCRIPTION", "FILTER_DESCRIPTION",
                    "FORMAT_DESCRIPTION"]
        for section in sections:
            if not ConfigUtils.hasSectionKey(self.config, section):
                raise ConfigInputIncompleteException("Missing %s section in the output config file." % section)

            if section in ("OTHER",):
                self._doFieldsExist(section, ["ID", "QUAL", "FILTER"])
            elif section in ("NOT_SPLIT_TAGS", "SPLIT_TAGS",):
                self._doFieldsExist(section, ["INFO", "FORMAT"])

    def _doFieldsExist(self, section, fields):
        """

        :param section:
        :param fields:
        :raise:
        """
        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, section)
        for field in fields:
            if field not in table:
                raise ConfigInputIncompleteException("Missing %s field in the %s section of the output config file."
                                                     % (field, section))

    def _parseConfig(self):
        """


        :return:
        """
        configTable = ConfigTable()

        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, "INFO")
        for ID, name in table.items():
            configTable.addInfoFieldID(ID, name)

        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, "FORMAT")
        for ID, name in table.items():
            configTable.addFormatFieldID(ID, name)

        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, "OTHER")
        for ID, name in table.items():
            configTable.addOtherFieldID(ID, name)

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "INFO_DESCRIPTION")
        for ID, desc in table.items():
            configTable.addInfoFieldIDDesc(ID, string.join(desc, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "FORMAT_DESCRIPTION")
        for ID, desc in table.items():
            configTable.addFormatFieldIDDesc(ID, string.join(desc, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "FILTER_DESCRIPTION")
        for ID, desc in table.items():
            configTable.addFilterFieldIDDesc(ID, string.join(desc, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "SPLIT_TAGS")
        for fieldType, IDs in table.items():
            for ID in IDs:
                configTable.addFieldIDToSplit(fieldType, ID)

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "NOT_SPLIT_TAGS")
        for fieldType, IDs in table.items():
            for ID in IDs:
                configTable.addFieldIDToNotSplit(fieldType, ID)
        return configTable

    def renderMutations(self, mutations, metadata=[], comments=[]):
        """ Generate a simple tsv file based on the incoming mutations.
        Assumes that all mutations have the same annotations, even if some are not populated.
        :param mutations:
        :param metadata:
        :param comments:
        Returns a file name. """

        self.logger.info("Rendering VCF output file: " + self._filename)

        # Initialize config table
        self._validateOutputConfigFile()
        self.configTable = self._parseConfig()

        # Initialize the data manager
        mut = None
        for mutation in mutations:
            mut = mutation
            lst = [(mut for mut in [mut]), mutations]
            mutations = itertools.chain(*lst)
            break

        fieldnames = self._getFieldnames(mut, metadata)

        path = os.getcwd()
        tempTsvFile = tempfile.NamedTemporaryFile(dir=path)  # create a temporary file to write tab-separated file
        self.logger.info("Creating intermediate tsv file...")
        self._writeMuts2Tsv(tempTsvFile.name, fieldnames, mutations)
        dm = OutputDataManager(self.configTable, comments, metadata, mut, self.sampleNames)

        self.logger.info("Intermediate tsv file created.")

        # Sort the tsv file
        chrom2HashCode = self._createChrom2HashCodeTable(self.chroms)
        tsvFileSorter = TsvFileSorter(tempTsvFile.name)
        sortedTempTsvFile = tempfile.NamedTemporaryFile(dir=path)
        func = lambda val: (chrom2HashCode[val["chr"]], int(val["start"]), val["alt_allele"])
        tsvFileSorter.sortFile(sortedTempTsvFile.name, func)
        self.logger.info("Intermediate tsv file sorted.")

        # Write the header
        filePointer = file(self._filename, 'w')
        header = dm.getHeader()
        filePointer.write(header)
        filePointer.close()

        self.logger.info("Render starting...")
        self._renderSortedTsv(self._filename, sortedTempTsvFile.name, self.sampleNames, dm)
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

    def _renderSortedTsv(self, vcfFilename, tsvFilename, sampleNames, dataManager):
        """

        :param vcfFilename:
        :param tsvFilename:
        :param sampleNames:
        :param dataManager:
        """
        tempVcfReader = vcf.Reader(filename=vcfFilename, strict_whitespace=True)
        pointer = file(vcfFilename, "w")
        vcfWriter = vcf.Writer(pointer, tempVcfReader, self.lineterminator)
        tsvReader = GenericTsvReader(tsvFilename, delimiter=self.delimiter)
        index = 0
        nrecords = 1000
        chrom = None
        pos = None
        recordFactory = None
        for m in tsvReader:
            isNewRecord = self._isNewVcfRecordNeeded(chrom, m["chr"], pos, m["start"])
            if isNewRecord:
                if recordFactory is not None:
                    record = recordFactory.createRecord()
                    vcfWriter.write_record(record)
                    index += 1
                    if index % nrecords == 0:
                        self.logger.info("Rendered " + str(index) + " vcf records.")
                        vcfWriter.flush()

                chrom = m["chr"]
                pos = m["start"]
                refAllele = m["ref_allele"]

                recordFactory = RecordFactory(chrom, int(pos), refAllele, sampleNames)

            recordFactory = self._parseRecordFactory(m, recordFactory, dataManager)

        if recordFactory is not None:
            record = recordFactory.createRecord()
            vcfWriter.write_record(record)

        vcfWriter.close()
        self.logger.info("Rendered all " + str(index) + " vcf records.")

    def _parseRecordFactory(self, m, recordFactory, dataManager):
        """

        :param m:
        :param recordFactory:
        :param dataManager:
        :return:
        """
        IDs = dataManager.getAnnotationNames("ID")
        quals = dataManager.getAnnotationNames("QUAL")
        filts = dataManager.getAnnotationNames("FILTER")
        infos = dataManager.getAnnotationNames("INFO")
        formats = dataManager.getAnnotationNames("FORMAT")

        altAllele = m["alt_allele"]
        sampleName = m.get("sampleName", None)
        recordFactory.addAlt(altAllele)

        for name in IDs:
            val = m.get(name, "")
            recordFactory.addID(val)

        qual = quals[0]
        recordFactory.addQual(m[qual])

        for name in filts:
            ID = dataManager.getFieldID(name)
            val = m.get(name, "")
            recordFactory.addFilter(ID, val)

        for name in infos:
            if name == "ESP_AvgSampleReadDepth":
                stop = True
            annotation = dataManager.getOutputAnnotation(name)
            ID = annotation.getID()
            num = annotation.getNumber()
            dataType = annotation.getDataType()
            isSplit = annotation.isSplit()
            val = m.get(name, "")
            recordFactory.addInfo(sampleName, ID, num, dataType, val, isSplit)

        for name in formats:
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
                recordFactory.addFormat(sampleName, ID, num, dataType, val, isSplit)

        return recordFactory

    def _createChrom2HashCodeTable(self, chroms):
        """

        :param chroms:
        :return:
        """
        table = dict()
        highestHashCode = 0
        sorted(chroms)
        for chrom in chroms:
            table[chrom] = None
            if chrom.isdigit():
                table[chrom] = int(chrom)
                if highestHashCode < table[chrom]:
                    highestHashCode = table[chrom]
        index = 0
        for chrom in chroms:
            if table[chrom] is None:
                if chrom.upper() == 'X':  # X chromosome
                    table[chrom] = highestHashCode + 1
                elif chrom.upper() == 'Y':  # Y chromosome
                    table[chrom] = highestHashCode + 2
                elif (chrom.upper() == 'M') or (chrom.upper() == 'MT'):  # mitochondrial chromosome
                    table[chrom] = highestHashCode + 3
                else:
                    index += 1
                    table[chrom] = highestHashCode + index + 3
        return table