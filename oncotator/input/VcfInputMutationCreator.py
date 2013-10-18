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
from oncotator.Metadata import Metadata


"""
Created on Oct 23, 2012

@author: gavlee
"""

from ConfigInputIncompleteException import ConfigInputIncompleteException
from oncotator.utils.ConfigUtils import ConfigUtils
from InputMutationCreator import InputMutationCreator
import logging
from oncotator.MutationData import MutationData
from oncotator.utils.MutUtils import MutUtils
import string
from SyntaxException import SyntaxException
import vcf
import copy
from oncotator.Annotation import Annotation
from oncotator.Metadata import Metadata


class VcfInputMutationCreator(InputMutationCreator):
    '''
    TODO: Finish documentation
    
    Convert a VCF 4.1 into a MutationData generator.
            Adds the following annotations (as INPUT annotation source):
                sampleName -- the name of the sample in the VCF
                altAlleleSeen -- whether the alternate allele was seen in the mutation.  This is whether a "1" appears in the GT field.
                
    
    '''

    def __init__(self, filename, configFile='vcf.in.config'):
        '''
        TODO: Pass in options: called-only -- only create mutations for called entries in the vcf
        '''
        self.filename = filename
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)
        self.config = ConfigUtils.createConfigParser(configFile, ignoreCase=False)
        self.configTable = None
        self.tags = {"INFO": ["aggregate"], "FORMAT": ["variant"], "FILTER": ["filter"], "ID": ["identifier"],
                     "QUAL": ["quality"]}
        self.logger = logging.getLogger(__name__)

    def _correctTable(self, table, inputIDs, metaInfoIDs, fieldType=""):
        """

        :param table:
        :param inputIDs:
        :param metaInfoIDs:
        :param fieldType:
        :return: :raise:
        """
        for key, value in table.items():
            if key not in metaInfoIDs:
                del table[key]
        diff = set(inputIDs).difference(set(metaInfoIDs))
        if len(diff) > 0:
            raise SyntaxException("Missing %s ID(s) (%s) in the meta-information specification of the VCF file, %s."
                                  % (fieldType, string.join(diff, ","), self.filename))
        for ID in inputIDs:
            if ID not in table:
                table[ID] = ID
        return table

    def _addGenotypeDataToMutation(self, mutation, record, index):
        """


        :param mutation:
        :param record:
        :param index:
        """
        IDs = []
        genotypeData = None

        if "FORMAT" in self.configTable:
            IDs = self.configTable["FORMAT"].keys()
            sampleName = mutation.getAnnotation("sampleName").getValue()
            genotypeData = record.genotype(sampleName)

        if record.FORMAT is not None:
            for ID in IDs:
                val = ""
                dataType = self.vcf_reader.formats[ID].type

                num = self.vcf_reader.formats[ID].num
                tags = copy.copy(self.tags["FORMAT"])

                if (genotypeData is not None) and (hasattr(genotypeData.data, ID)):
                    isSplitTag = self._determineIsSplit(ID, num, "FORMAT")
                    if isSplitTag:
                        val = genotypeData[ID][index]
                    else:
                        val = genotypeData[ID]

                if (val is None) or (val == ""):
                    if dataType == "Flag":
                        val = "False"
                    else:
                        val = ""
                elif isinstance(val, list):
                    val = string.join(["" if v is None else str(v) for v in val], ",")
                else:
                    val = str(val)

                name = self._getAnnotationName("FORMAT", ID)
                mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.formats[ID].desc, tags=tags,
                                          number=num)

        return mutation

    def _addInfoDataToMutation(self, mutation, record, index):
        """

        :param mutation:
        :param record:
        :param index:
        :return:
        """
        IDs = []
        if "INFO" in self.configTable:
            IDs = self.configTable["INFO"].keys()

        for ID in IDs:
            val = ""
            dataType = self.vcf_reader.infos[ID].type

            num = self.vcf_reader.infos[ID].num
            tags = copy.copy(self.tags["INFO"])

            if ID in record.INFO:
                isSplitTag = self._determineIsSplit(ID, num, "INFO")
                if isSplitTag:
                    val = record.INFO[ID][index]
                else:
                    val = record.INFO[ID]

            if (val is None) or (val == ""):
                if dataType == "Flag":
                    val = "False"
                else:
                    val = ""
            elif isinstance(val, list):
                val = string.join(["" if v is None else str(v) for v in val], ",")
            else:
                val = str(val)

            name = self._getAnnotationName("INFO", ID)
            mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.infos[ID].desc, tags=tags,
                                      number=num)

        return mutation

    def _determineIsSplit(self, ID, num, fieldType):
        if num == -2:  # by the number of samples
            split = False
            if ID in self.configTable["SPLIT_TAGS"]["FORMAT"]:  # override the default using the config file
                split = True
        elif num == -1:  # by the number of alternates
            split = True
            if ID in self.configTable["NOT_SPLIT_TAGS"][fieldType]:  # override the default using the config file
                split = False
        elif num == 0:
            split = False
        elif num is None:
            split = False
            if ID in self.configTable["SPLIT_TAGS"][fieldType]:  # override the default using the config file
                split = True
        else:
            split = False

        return split

    def _createConfigTableKeys(self):
        # Parse fields from INFO section of the config file
        if not ConfigUtils.hasSectionKey(self.config, "INFO"):
            raise ConfigInputIncompleteException("Missing %s section in input config file for the input vcf file %s."
                                                 % ("INFO", self.filename))
        self.configTable["INFO"] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, "INFO")

        # Parse fields from FORMAT section of the config file
        if not ConfigUtils.hasSectionKey(self.config, "FORMAT"):
            raise ConfigInputIncompleteException("Missing %s section in input config file for the input vcf file %s."
                                                 % ("FORMAT", self.filename))
        self.configTable["FORMAT"] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, "FORMAT")

        # Parse fields from NOT_SPLIT_TAGS section of the config file
        if not ConfigUtils.hasSectionKey(self.config, "NOT_SPLIT_TAGS"):
            raise ConfigInputIncompleteException("Missing %s section in input config file for the input vcf file %s."
                                                 % ("NOT_SPLIT_TAGS", self.filename))
        self.configTable["NOT_SPLIT_TAGS"] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config,
                                                                                               "NOT_SPLIT_TAGS")

        # Parse fields from SPLIT_TAGS section of the config file
        if not ConfigUtils.hasSectionKey(self.config, "SPLIT_TAGS"):
            raise ConfigInputIncompleteException("Missing %s section in input config file for the input vcf file %s."
                                                 % ("SPLIT_TAGS", self.filename))
        self.configTable["SPLIT_TAGS"] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, "SPLIT_TAGS")

    def _createConfigTable(self):
        # Iterate over the input VCF file and parse all possible fields in INFO and FORMAT section (excluding the
        # fields listed in the meta-information section)
        infos = set()
        formats = set()
        for variant in self.vcf_reader:  # parse INFO and FORMAT fields per record
            infos = infos.union(set(variant.INFO.keys()))
            for sample in variant.samples:
                formats = formats.union(set(sample.data._fields))

        self.configTable["INFO"] = self._correctTable(self.configTable["INFO"], infos, self.vcf_reader.infos.keys(),
                                                      "INFO")
        self.configTable["FORMAT"] = self._correctTable(self.configTable["FORMAT"], formats,
                                                        self.vcf_reader.formats.keys(), "FORMAT")

    def createMutations(self):
        """ Creates a mutation for each mutation by each sample, regardless of allelic depth, etc.
            
            Annotations that this will generate (as source = "INPUT"):
                sampleName
                isCalled
                isAltAllele
                allelic_depth -- DP in the vcf

            TODO: Complete documentation
        """
        self.reset()
        self.configTable = dict()
        self._createConfigTableKeys()
        self._createConfigTable()

        self.reset()
        for record in self.vcf_reader:
            for flt in record.FILTER:
                if flt not in self.vcf_reader.filters:
                    raise SyntaxException(
                        "Missing FILTER ID, %s, in the meta-information specification of the VCF file, %s."
                        % (flt, self.filename))

            for index in range(len(record.ALT)):
                mut = self._createMutation(record, index)

                if len(record.samples) <= 0:
                    yield mut
                else:
                    for sample in record.samples:
                        # TODO: move this to deep copy
                        sampleMut = self._createMutationCopy(mut)
                        sampleMut.createAnnotation("sampleName", sample.sample, "INPUT")
                        genotype = "GT"
                        if genotype in sample.data._fields:
                            if (sample.data.GT is None) or (sample.data.GT.find("1") == -1):
                                sampleMut["altAlleleSeen"] = False
                        sampleMut = self._addGenotypeDataToMutation(sampleMut, record, index)

                        yield sampleMut

    def _createMutationCopy(self, mutation):
        """
        """
        chrom = mutation["chr"]
        alt = mutation["alt_allele"]
        ref = mutation["ref_allele"]
        startPos = mutation["start"]
        endPos = mutation["end"]
        mut = MutationData(chrom, alt, ref, startPos, endPos, "hg19")
        for annotationName in mutation.annotations:
            annotation = mutation.getAnnotation(annotationName)
            mut.createAnnotation(annotationName, annotation.getValue(), annotation.getDatasource(),
                                 annotation.getDataType(), annotation.getDescription(), True, annotation.getTags(),
                                 annotation.getNumber())
        return mut

    def _createMutation(self, record, index):
        """
        """
        chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)

        alt = record.ALT[index]
        if alt is None:
            alt = ""
        else:
            alt = str(alt)

        endPos = int(record.POS) + len(alt) - 1
        ref = record.REF
        if ref == ".":
            ref = ""
        mut = MutationData(chrom, record.POS, endPos, ref, alt, "hg19")

        ID = record.ID
        if ID is None:
            ID = ""
        mut.createAnnotation("id", ID, "INPUT", tags=copy.copy(self.tags["ID"]))

        mut.createAnnotation("qual", str(record.QUAL), "INPUT", tags=copy.copy(self.tags["QUAL"]))
        for filt in self.vcf_reader.filters:  # for each filter in the header
            description = self.vcf_reader.filters[filt].desc  # parse the description
            if (len(record.FILTER) != 0) and \
                    (filt in record.FILTER):  # if the filter is mentioned for this variant, then it failed
                mut.createAnnotation(filt, "FAIL", "INPUT", annotationDescription=description,
                                     tags=copy.copy(self.tags["FILTER"]))
            else:
                mut.createAnnotation(filt, "PASS", "INPUT", annotationDescription=description,
                                     tags=copy.copy(self.tags["FILTER"]))
        mut.createAnnotation("altAlleleSeen", str(True), "INPUT")
        mut = self._addInfoDataToMutation(mut, record, index)
        return mut

    def _getAnnotationName(self, fieldType, ID):
        if ID not in self.configTable[fieldType]:
            raise SyntaxException(
                "Missing %s ID, %s, in meta-information of VCF file, %s." % (fieldType, ID, self.filename))
        return self.configTable[fieldType][ID]

    def reset(self):
        """ Resets the internal state, so that mutations can be generated. """
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)

    def getComments(self):
        """ Comments often need to be passed into the output.  Get the comments from the input file."""
        comments = list()
        keys = self.vcf_reader.metadata.keys()
        for key in keys:
            val = self.vcf_reader.metadata[key]
            if isinstance(val, list):
                val = string.join(map(str, val), ";")
            if isinstance(val, dict):
                val = string.join(map(str, val), ";")
            comment = string.join([key, val], "=")
            comments.append(comment)
        return comments

    def _addFormatFields2Metadata(self, metadata):
        for ID, annotationName in self.configTable["FORMAT"].iteritems():
            num = self.vcf_reader.formats[ID].num
            tags = copy.copy(self.tags["FORMAT"])
            isSplitTag = self._determineIsSplit(ID, num, "FORMAT")
            if isSplitTag:
                tags += ["SPLIT"]
            metadata[annotationName] = Annotation("", "INPUT", self.vcf_reader.formats[ID].type,
                                                  self.vcf_reader.formats[ID].desc, tags=tags, number=num)
        return metadata

    def _addInfoFields2Metadata(self, metadata):
        for ID, annotationName in self.configTable["INFO"].iteritems():
            num = self.vcf_reader.infos[ID].num
            tags = copy.copy(self.tags["INFO"])
            isSplitTag = self._determineIsSplit(ID, num, "INFO")
            if isSplitTag:
                tags += ["SPLIT"]
            metadata[annotationName] = Annotation("", "INPUT", self.vcf_reader.infos[ID].type,
                                                  self.vcf_reader.infos[ID].desc, tags=tags, number=num)
        return metadata

    def _addFilterFields2Metadata(self, metadata):
        for filt in self.vcf_reader.filters:  # for each filter in the header
            metadata[filt] = Annotation("", "INPUT", "String", self.vcf_reader.filters[filt].desc,
                                        tags=copy.copy(self.tags["FILTER"]))
        return metadata

    def _createMetadata(self):
        self.configTable = dict()
        self._createConfigTableKeys()
        self._createConfigTable()

        metadata = Metadata()
        metadata = self._addFilterFields2Metadata(metadata)
        metadata = self._addFormatFields2Metadata(metadata)
        metaData = self._addInfoFields2Metadata(metadata)

        metaData["id"] = Annotation("", "INPUT", "String", "", copy.copy(self.tags["ID"]))
        metaData["qual"] = Annotation("", "INPUT", "String", "", copy.copy(self.tags["QUAL"]))

        return metadata

    def getMetadata(self):
        metadata = self._createMetadata()
        return metadata
