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


'''
Created on Oct 23, 2012

@author: gavlee
'''

from ConfigInputIncompleteException import ConfigInputIncompleteException
from oncotator.utils.ConfigUtils import ConfigUtils
from InputMutationCreator import InputMutationCreator
import logging
from oncotator.MutationData import MutationData
from oncotator.utils.MutUtils import MutUtils
import string
from SyntaxException import SyntaxException
import vcf
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
        self.tags = {"INFO": ["aggregate"], "FORMAT": ["variant"], "FILTER": ["filter"]}
        self.logger = logging.getLogger(__name__)

    def __getUnspecifiedIDs(self, s, t):
        return set(s).difference(set(t))

    def __createInfoTable(self, infos, infoIDs):
        """
        """
        for k, v in infos.items(): # compare parsed config inputs to meta-information specifications
            if not k in self.vcf_reader.infos: # delete INFO field when not found in the meta-information section of the input VCF file
                del infos[k]
            else:
                infos[k] = v

        unspecifiedInfoIDs = self.__getUnspecifiedIDs(s=infoIDs, t=self.vcf_reader.infos.keys()) # set of IDs that are not present in VCF meta-information
        if len(unspecifiedInfoIDs) > 0:
            raise SyntaxException(
                "Missing INFO ID(s) (%s) in the meta-information specification of the VCF file, " % ",".join(
                    unspecifiedInfoIDs) + self.filename + ".")

        for ID in infoIDs:
            if ID not in infos:
                infos[ID] = ID

        return infos

    def __createFormatTable(self, formats, formatIDs):
        """
        """
        for k, v in formats.items():  # compare parsed config inputs to meta-information specification
            if not k in self.vcf_reader.formats:  # delete FORMAT field when not found in the meta-information section of the input VCF file
                del formats[k]
            else:
                formats[k] = v

        unspecifiedFormatsIDs = self.__getUnspecifiedIDs(s=formatIDs, t=self.vcf_reader.formats.keys())
        if len(unspecifiedFormatsIDs) > 0:  # set of IDs that are not present in VCF meta-information
            raise SyntaxException(
                "Missing FORMAT ID(s) (%s) in the meta-information specification of the VCF file, " % ",".join(
                    unspecifiedFormatsIDs) + self.filename + ".")

        for ID in formatIDs:
            if ID not in formats:
                formats[ID] = ID

        return formats

    def __addGenotypeDataToMutation(self, mutation, record, index):
        """

        """
        formatIDs = []
        sampleName = None
        sampleGenotype = None

        if "FORMAT" in self.configTable:
            formatIDs = self.configTable["FORMAT"].keys()
            sampleName = mutation.getAnnotation("sampleName").getValue()
            sampleGenotype = record.genotype(sampleName)

        if record.FORMAT is not None:
            for ID in formatIDs:
                dataType = self.vcf_reader.formats[ID].type
                value = ""

                number = "."
                if not self.vcf_reader.formats[ID].num is None:
                    number = self.vcf_reader.formats[ID].num

                if ID in sampleGenotype.data._fields:
                    if ("format" in self.configTable["SPLIT_TAGS"]) and \
                            (ID in self.configTable["SPLIT_TAGS"]["format"]):
                        if (isinstance(sampleGenotype[ID], list) and len(record.ALT) != len(sampleGenotype[ID])) or \
                                (not isinstance(sampleGenotype[ID], list) and len(record.ALT) > 1 and
                                    not sampleGenotype[ID] is None):
                            raise Exception("Number of alternative alleles at chromosome %s and position %s for "
                                            "sample %s does not match the number of values specified by the %s "
                                            "FORMAT ID." % (record.CHROM, record.POS, sampleName, ID))
                        elif isinstance(sampleGenotype[ID], list) and (not sampleGenotype[ID][index] is None):
                            value = str(sampleGenotype[ID][index])
                        elif not sampleGenotype[ID] is None:
                            value = str(sampleGenotype[ID])
                    else:
                        if isinstance(sampleGenotype[ID], list):
                            value = string.join(["" if v is None else str(v) for v in sampleGenotype[ID]], ",")
                        elif not sampleGenotype[ID] is None:
                            value = sampleGenotype[ID]
                if (dataType == "Flag") and (not value):
                    value = str(False)

                mutation.createAnnotation(annotationName=self.__getAnnotationName("FORMAT", ID), annotationValue=value,
                                          annotationSource="INPUT", annotationDataType=dataType,
                                          annotationDescription=self.vcf_reader.formats[ID].desc, number=number,
                                          tags=self.tags["FORMAT"])

        return mutation

    def __addInfoDataToMutation(self, mutation, record, index):

        infoIDs = []
        if "INFO" in self.configTable:
            infoIDs = self.configTable["INFO"].keys()

        for ID in infoIDs:
            value = ""
            dataType = self.vcf_reader.infos[ID].type

            number = "."
            if not self.vcf_reader.infos[ID].num is None:
                number = self.vcf_reader.infos[ID].num

            if ID in record.INFO:
                if ("info" in self.configTable["SPLIT_TAGS"]) and (ID in self.configTable["SPLIT_TAGS"]["info"]):
                    if (isinstance(record.INFO[ID], list) and len(record.ALT) != len(record.INFO[ID])) or \
                            (not isinstance(record.INFO[ID], list) and len(record.ALT) > 1):
                        raise Exception("Number of alternative alleles at chromosome %s and position %s does not "
                                        "match the number of values specified by the %s INFO ID."
                                        % (record.CHROM, record.POS, ID))
                    elif isinstance(record.INFO[ID], list) and (not record.INFO[ID][index] is None):
                        value = str(record.INFO[ID][index])
                    elif not record.INFO[ID] is None:
                        value = str(record.INFO[ID])
                else:
                    if isinstance(record.INFO[ID], list):
                        value = string.join(["" if v is None else str(v) for v in record.INFO[ID]], ",")
                    elif not record.INFO[ID] is None:
                        value = str(record.INFO[ID])

            if (dataType == "Flag") and (not value):
                value = str(False)

            mutation.createAnnotation(annotationName=self.__getAnnotationName("INFO", ID), annotationValue=value,
                                      annotationSource="INPUT", annotationDataType=dataType,
                                      annotationDescription=self.vcf_reader.infos[ID].desc, number=number,
                                      tags=self.tags["INFO"])

        return mutation

    def __createConfigTableKeys(self):
        sections = ["INFO", "FORMAT", "SPLIT_TAGS"]
        for section in sections:
            if not ConfigUtils.hasSectionKey(self.config, sectionKey=section):
                raise ConfigInputIncompleteException(
                    "Missing %s section in input config file for input file %s." % (section, self.filename))

        # Parse fields from INFO section of the config file
        self.configTable["INFO"] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser=self.config,
                                                                                           sectionKey="INFO")
        # Parse fields from FORMAT section of the config file
        self.configTable["FORMAT"] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser=self.config,
                                                                                             sectionKey="FORMAT")
        # Parse fields from SPLIT_TAGS section of the config file
        self.configTable["SPLIT_TAGS"] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser=self.config,
                                                                                           sectionKey="SPLIT_TAGS")

    def __createConfigTable(self):
        # Iterate over the input VCF file and parse all possible fields in INFO and FORMAT section (excluding the
        # fields listed in the meta-information section)
        infoIDs = set()
        formatIDs = set()
        for variant in self.vcf_reader:  # parse INFO and FORMAT fields per record
            for fieldID in variant.INFO:
                infoIDs.add(fieldID)
            if variant.FORMAT is not None:
                formatIDs = formatIDs.union(set(variant.FORMAT.split(':')))

        self.configTable["INFO"] = self.__createInfoTable(self.configTable["INFO"], infoIDs)
        self.configTable["FORMAT"] = self.__createFormatTable(self.configTable["FORMAT"], formatIDs)

    def createMutations(self):
        """ Creates a mutation for each mutation by each sample, regardless of allelic depth, etc.
            
            Annotations that this will generate (as source = "INPUT"):
                sampleName
                isCalled
                isAltAllele
                allelic_depth -- DP in the vcf

            TODO: Complete documentation
        """
        if self.configTable is None:
            self.configTable = dict()
            self.__createConfigTableKeys()
            self.__createConfigTable()

        self.reset()
        for record in self.vcf_reader:
            if record.ID is None:
                record.ID = ""
            elif isinstance(record.ID, list):
                record.ID = string.join(map(str, record.ID), ';')

            if record.REF is None:  # replace "None" with empty string for REF
                record.REF = ""

            if record.FILTER is None:  # FILTER is either "PASS" or semicolon separated list of alphanumerics
                record.FILTER = ""
            elif isinstance(record.FILTER, list):
                for ID in record.FILTER:
                    if ID not in self.vcf_reader.filters:
                        raise SyntaxException(
                            "Missing FILTER field ID, %s, in the meta-information specification of the VCF file, " % ID
                            + self.filename + ".")

            for index in range(len(record.ALT)):
                mut = self.__createMutation(record=record, index=index)

                if len(record.samples) <= 0:
                    yield mut
                else:
                    for sample in record.samples:
                        # TODO: move this to deep copy
                        sampleMut = self.__createMutationCopy(mutation=mut)
                        sampleMut.createAnnotation(annotationName="sampleName", annotationValue=sample.sample,
                                                   annotationSource="INPUT")
                        # Check if the alt allele is actually seen
                        genotype = "GT"
                        if hasattr(sample.data, genotype):
                            if (sample.data.GT is None) or (sample.data.GT.find("1") == -1):
                                sampleMut["altAlleleSeen"] = False
                        sampleMut = self.__addGenotypeDataToMutation(mutation=sampleMut, record=record, index=index)

                        yield sampleMut

    def __createMutationCopy(self, mutation):
        """
        """
        chrom = mutation["chr"]
        alt = mutation["alt_allele"]
        ref = mutation["ref_allele"]
        startPos = mutation["start"]
        endPos = mutation["end"]
        mut = MutationData(chr=chrom, alt_allele=alt, ref_allele=ref, start=startPos, end=endPos, build="hg19")
        for annotationName in mutation.annotations:
            annotation = mutation.getAnnotation(annotationName=annotationName)
            mut.createAnnotation(annotationName=annotationName, annotationValue=annotation.getValue(),
                                 annotationSource=annotation.getDatasource(),
                                 annotationDataType=annotation.getDataType(),
                                 annotationDescription=annotation.getDescription(), number=annotation.getNumber(),
                                 tags=annotation.getTags())
        return mut

    def __createMutation(self, record, index):
        """
        """
        chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
        alt = record.ALT[index]
        if alt is None:
            alt = ""
        endPos = int(record.POS) + len(alt) - 1
        mut = MutationData(chr=chrom, alt_allele=alt, ref_allele=record.REF, start=record.POS, end=endPos, build="hg19")
        mut.createAnnotation(annotationName="id", annotationValue=record.ID, annotationSource="INPUT")  # dbSNP variants
        mut.createAnnotation(annotationName="qual", annotationValue=str(record.QUAL), annotationSource="INPUT")
        for fltr in self.vcf_reader.filters:  # for each filter in the header
            description = self.vcf_reader.filters[fltr].desc  # parse the description
            if (len(record.FILTER) != 0) and (
                    fltr in record.FILTER):  # if the filter is mentioned for this variant, then it failed
                mut.createAnnotation(annotationName=fltr, annotationValue="FAIL", annotationSource="INPUT",
                                     annotationDescription=description, tags=self.tags["FILTER"])
            else:
                mut.createAnnotation(annotationName=fltr, annotationValue="PASS", annotationSource="INPUT",
                                     annotationDescription=description, tags=self.tags["FILTER"])
        mut.createAnnotation(annotationName="altAlleleSeen", annotationValue=str(True), annotationSource="INPUT")
        mut = self.__addInfoDataToMutation(mutation=mut, record=record, index=index)
        return mut

    def __getAnnotationName(self, section, ID):
        if ID not in self.configTable[section]:
            raise SyntaxException(
                "Missing %s ID, %s, in meta-information of VCF file, %s." % (section, ID, self.filename))
        return self.configTable[section][ID]

    def reset(self):
        """ Resets the internal state, so that mutations can be generated. """
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)

    def getComments(self):
        """ Comments often need to be passed into the output.  Get the comments from the input file."""
        comments = list()
        for k, v in self.vcf_reader.metadata.items():
            if isinstance(v, list):
                v = ";".join(str(k) for k in v)
            comment = k + "=" + v
            comments.append(comment)
        return comments

    def __addFormatMetaData(self, metaData):
        for ID, annotationName in self.configTable["FORMAT"].iteritems():
            number = "."
            if not self.vcf_reader.formats[ID].num is None:
                number = self.vcf_reader.formats[ID].num
            metaData[annotationName] = Annotation(value="", datasourceName="INPUT",
                                                  dataType=self.vcf_reader.formats[ID].type,
                                                  description=self.vcf_reader.formats[ID].desc, number=number,
                                                  tags=self.tags["FORMAT"])
        return metaData

    def __addInfoMetaData(self, metaData):
        for ID, annotationName in self.configTable["INFO"].iteritems():
            number = "."
            if not self.vcf_reader.infos[ID].num is None:
                number = self.vcf_reader.infos[ID].num
            metaData[annotationName] = Annotation(value="", datasourceName="INPUT",
                                                  dataType=self.vcf_reader.infos[ID].type,
                                                  description=self.vcf_reader.infos[ID].desc, number=number,
                                                  tags=self.tags["INFO"])
        return metaData

    def __addFilterMetaData(self, metaData):
        for fltr in self.vcf_reader.filters:  # for each filter in the header
            metaData[fltr] = Annotation(value="", datasourceName="INPUT",
                                        description=self.vcf_reader.filters[fltr].desc, tags=self.tags["FILTER"])
        return metaData

    def getMetadata(self):
        if self.configTable is None:
            self.configTable = dict()
            self.__createConfigTableKeys()
            self.__createConfigTable()

        metaData = Metadata()
        metaData = self.__addFilterMetaData(metaData=metaData)
        metaData = self.__addFormatMetaData(metaData=metaData)
        metaData = self.__addInfoMetaData(metaData=metaData)

        return metaData
