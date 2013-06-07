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


'''
Created on Nov 7, 2012

@author: lichtens
'''
from OutputRenderer import OutputRenderer
from oncotator.utils.ConfigUtils import ConfigUtils
import copy
import csv
import collections
import logging
import os
import sys
import tempfile
import itertools
import operator
import string
import vcf
import heapq
from oncotator.utils.MutUtils import MutUtils
from oncotator.input.ConfigInputIncompleteException import ConfigInputIncompleteException
from oncotator.utils.TsvFileSorter import TsvFileSorter


class VCFOutputDataManager:
    comments = None

    def __init__(self, configTable, comments=[], metadata=[], mut=None):
        comments = comments
        self.annotationTable = dict()
        self.table = configTable
        pass

    def _resolveFieldType(self, name, tags):
        fieldType = None
        if fieldType is None:
            value = None
            if isinstance(tags, list) and (len(tags) > 0):
                value = tags[0]  # we assume that the first tag corresponds to field type
            if value == "aggregate":
                fieldType = "INFO"
            elif value == "variant":
                fieldType = "FORMAT"
            elif value == "filter":
                fieldType = "FILTER"
            elif value == "identifier":
                fieldType = "ID"
            elif value == "quality":
                fieldType = "QUAL"

        if fieldType is None:
            if name in self.table["INFO"]:
                fieldType = "INFO"
            elif name in self.table["FORMAT"]:
                fieldType = "FORMAT"
            elif name in self.table["OTHER"]:
                fieldType = self.table["OTHER"][name]

        if fieldType is None:
            if name.upper() in vcf.parser.RESERVED_INFO:
                fieldType = "INFO"
            elif name.upper() in vcf.parser.RESERVED_FORMAT:
                fieldType = "FORMAT"

        if fieldType is None:
            fieldType = "FORMAT"

        return fieldType

    def _resolveFieldID(self, fieldType, name):
        ID = None
        if (fieldType in self.table) and (name in self.table[fieldType]):
            ID = self.table[fieldType][ID]
        if ID is None:
            ID = name
        return ID

    def _resolveFieldDataType(self, fieldType, ID, dataType):
        if (dataType == "String") or (dataType == ".") or (dataType == ""):
            if fieldType == "FILTER":
                if ID in vcf.parser.RESERVED_INFO:
                    dataType = vcf.parser.RESERVED_INFO[ID]
            elif fieldType == "FORMAT":
                if ID in vcf.parser.RESERVED_FORMAT:
                    dataType = vcf.parser.RESERVED_FORMAT[ID]
            else:
                dataType = "String"
        return dataType

    def _resolveFieldDescription(self, fieldType, ID, desc):
        # description in the config files overwrite the description in the mutation
        if (fieldType == "FILTER") and (ID in self.table["FILTER_DESCRIPTION"]):
            desc = self.table["FILTER_DESCRIPTION"][ID]
        elif (fieldType == "FORMAT") and (ID in self.table["FORMAT_DESCRIPTION"]):
            desc = self.table["FORMAT_DESCRIPTION"][ID]
        elif (fieldType == "INFO") and (ID in self.table["INFO_DESCRIPTION"]):
            desc = self.table["INFO_DESCRIPTION"][ID]
        if desc == "":
            desc = "Unknown"
        return desc

    def _annotation2str(self, fieldType, ID, desc="Unknown", dataType="String", num="."):
        if (fieldType == "FORMAT") or (fieldType == "INFO"):
            return "%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (fieldType, ID, num, dataType, desc)
        elif fieldType == "FILTER":
            return "%s=<ID=%s,Description=\"%s\">" % (fieldType, ID, desc)
        return ""

    def createMetaInfoHeader(self, names):
        headers = []
        flg = False
        if (not self.comments is None) and (len(self.comments) > 0):
            for idx in xrange(len(self.comments)):
                comment = self.comments[idx]
                if comment.startswith("fileformat=VCFv4."):
                    headers = [string.join(["##", comment], "")]
                    del self.comments[idx]
                    flg = True
                    break
        if not flg:
            headers = ["##fileformat=VCFv4.1"]
        if (not self.comments is None) and (len(self.comments) > 0):
            headers += self.comments[0:len(self.comments)-1]

        if len() == 0:
                annotationNames = metadata.keys()

        for annotationName in annotationNames:
            annotation = metadata[annotationName]
            field = self.__determineField(annotationName=annotationName, annotationTags=annotation.getTags())
            if (field != "ID") and (field != "QUAL"):
                ID = self.__determineFieldID(annotationName=annotationName, field=field)
                num = annotation.getNumber()
                dataType = self.__determineDataType(field=field, ID=ID, annotation=annotation)
                desc = self.__determineDescription(field=field, ID=ID, annotation=annotation)
                headers += [self.__annotation2str(field=field, ID=ID, desc=desc, dataType=dataType, num=num)]

        headers += [string.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'],
                                self.delimiter)]
        return headers



class VcfOutputRenderer(OutputRenderer):
    """
    The SimpleOutputRenderer renders a basic tsv file from the given mutations.  All annotations are included with real names as column headers.
    
    Header is determined by the first mutation given.
    
    No attention is paid to order of the headers.    
    """
    _vcfAnnotation = collections.namedtuple(typename="Annotation", field_names=["field", "ID", "num", "type"])

    def __init__(self,filename, datasources=[], configFile='vcf.out.config'):
        """
        Constructor
        """
        self._filename = filename
        self.logger = logging.getLogger(__name__)
        self._datasources = datasources
        self.config = ConfigUtils.createConfigParser(configFile,ignoreCase=False)
        self.chromHashCodeTable = None
        self.configTable = dict()
        self.delimiter = '\t'

        self.reservedAnnotationNames = ['chr', 'start', 'sampleName', 'ref_allele', 'alt_allele', 'end', 'build',
                                        'altAlleleSeen']

    # def __determineField(self, annotationName, annotationTags):
    #     field = None
    #     if field is None: # first, determine field ID using annotation's tags field
    #         if isinstance(annotationTags, list) and (len(annotationTags) > 0):
    #             field = annotationTags[0]
    #         if (not field is None) and (field == "aggregate"):
    #             field = "INFO"
    #         elif (not field is None) and (field == "variant"):
    #             field = "FORMAT"
    #         elif (not field is None) and (field == "filter"):
    #             field = "FILTER"
    #         elif (not field is None) and (field == "identifier"):
    #             field = "ID"
    #         elif (not field is None) and (field == "quality"):
    #             field = "QUAL"
    #
    #     if field is None: # second, determine field ID using config file
    #         if annotationName in self.configTable["INFO"]:
    #             field = "INFO"
    #         elif annotationName in self.configTable["FORMAT"]:
    #             field = "FORMAT"
    #         elif annotationName in self.configTable["OTHER"]:
    #             field = self.configTable["OTHER"][annotationName]
    #
    #     if field is None: # third, determine field ID using PyVCF defaults
    #         if annotationName.upper() in vcf.parser.RESERVED_INFO:
    #             field = "INFO"
    #         elif annotationName.upper() in vcf.parser.RESERVED_FORMAT:
    #             field = "FORMAT"
    #
    #     if field is None: # default: FORMAT
    #         field = "FORMAT"
    #
    #     return field

    # def __determineFieldID(self, annotationName, field):
    #     ID = annotationName
    #     if (field in self.configTable) and (annotationName in self.configTable[field]):
    #         ID = self.configTable[field][annotationName]
    #     return ID

    # def __determineDataType(self, field, ID, annotation):
    #     dataType = annotation.getDataType()
    #
    #     if (dataType == "String") or (dataType == ".") or (not dataType):
    #         if field == "FILTER":
    #             if ID in vcf.parser.RESERVED_INFO:
    #                 dataType = vcf.parser.RESERVED_INFO[ID]
    #         elif field == "FORMAT":
    #             if ID in vcf.parser.RESERVED_FORMAT:
    #                 dataType = vcf.parser.RESERVED_FORMAT[ID]
    #         else:
    #             dataType = "String"
    #
    #     return dataType
    #
    # def __determineDescription(self, field, ID, annotation):
    #     desc = annotation.getDescription()
    #
    #     # description in the config files overwrite the description in the mutation
    #     if (field == "FILTER") and \
    #             ("FILTER_DESCRIPTION" in self.configTable and ID in self.configTable["FILTER_DESCRIPTION"]):
    #         desc = self.configTable["FILTER_DESCRIPTION"][ID]
    #     elif (field == "FORMAT") and \
    #             ("FORMAT_DESCRIPTION" in self.configTable and ID in self.configTable["FORMAT_DESCRIPTION"]):
    #         desc = self.configTable["FORMAT_DESCRIPTION"][ID]
    #     elif (field == "INFO") and \
    #             ("INFO_DESCRIPTION" in self.configTable and ID in self.configTable["INFO_DESCRIPTION"]):
    #         desc = self.configTable["INFO_DESCRIPTION"][ID]
    #
    #     if desc == "":
    #         desc = "Unknown"
    #     return desc

    def __createAnnotationTableFromMetaData(self, metadata, mutation):
        annotationNames = set(mutation.keys()).intersection(metadata.keys())
        annotationTable = dict()
        for annotationName in annotationNames:
            annotation = metadata[annotationName]
            field = self.__determineField(annotationName=annotationName, annotationTags=annotation.getTags())
            ID = self.__determineFieldID(annotationName=annotationName, field=field)
            dataType = self.__determineDataType(field=field, ID=ID, annotation=annotation)
            num = annotation.getNumber()
            annotationTable[annotationName] = self._vcfAnnotation(field=field, ID=ID, num=num, type=dataType)
        return annotationTable

    def __createMetaInformationHeader(self, metadata, comments, annotationNames):
        headers = []
        flg = False
        if (not comments is None) and (len(comments) > 0):
            for idx in xrange(len(comments)):
                comment = comments[idx]
                if comment.startswith("fileformat=VCFv4."):
                    headers = [string.join(["##", comment], "")]
                    del comments[idx]
                    flg = True
                    break
        if not flg:
            headers = ["##fileformat=VCFv4.1"]
        if (not comments is None) and (len(comments) > 0):
            headers += comments[0:len(comments)-1]

        if len(annotationNames) == 0:
                annotationNames = metadata.keys()

        for annotationName in annotationNames:
            annotation = metadata[annotationName]
            field = self.__determineField(annotationName=annotationName, annotationTags=annotation.getTags())
            if (field != "ID") and (field != "QUAL"):
                ID = self.__determineFieldID(annotationName=annotationName, field=field)
                num = annotation.getNumber()
                dataType = self.__determineDataType(field=field, ID=ID, annotation=annotation)
                desc = self.__determineDescription(field=field, ID=ID, annotation=annotation)
                headers += [self.__annotation2str(field=field, ID=ID, desc=desc, dataType=dataType, num=num)]

        headers += [string.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'],
                                self.delimiter)]
        return headers

    # def __annotation2str(self, field, ID, desc="Unknown", dataType="String", num="."):
    #     if (field == "FORMAT") or (field == "INFO"):
    #         return "%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (field, ID, num, dataType, desc)
    #     elif field == "FILTER":
    #         return "%s=<ID=%s,Description=\"%s\">" % (field, ID, desc)
    #     return ""

    def __writeMutationsToFile(self, filename, mutations, metadata):
        sampleNames = set()
        chroms = set()
        annotationTable = dict()
        annotationNames = []

        fp = open(filename, 'w')
        dw = None
        for mut in mutations:
            # Sample names (set)
            if "sampleName" in mut:
                sampleName = mut.getAnnotation('sampleName').getValue()
                if not sampleName in sampleNames:
                    sampleNames.add(sampleName)

            # Parse chromosome
            chrom = mut.getAnnotation('chr').getValue()
            if not chrom in chroms:
                chroms.add(chrom)

            if (len(annotationTable) == 0) and (len(annotationNames) == 0):
                annotationTable = self.__createAnnotationTableFromMetaData(metadata=metadata, mutation=mut)
                annotationNames = self.reservedAnnotationNames + annotationTable.keys()

            if dw is None: # initialize the CSV dictionary writer
                dw = csv.DictWriter(fp, annotationNames, delimiter=self.delimiter, lineterminator="\n")
                dw.writeheader()

            dw.writerow(mut)
        fp.close()

        return sampleNames, chroms, annotationTable

    def __getAnnotationNamesByField(self, tbl, fld):
        names = []
        for name in tbl:
            annotation = tbl[name]
            if annotation.field == fld:
                names.append(name)
        return names

    def __map(self, func, iterable, bad="."):
        return [func(x) if x != bad else None
                for x in iterable]

    def __parseValue(self, val, dt, num):
        vals = val.split(",")
        vals = ["." if not x else x
                for x in vals]
        if dt == "Integer":
            try:
                val = self.__map(int, vals)
            except ValueError:
                val = self.__map(float, vals)
        elif dt == "Float":
            val = self.__map(float, vals)
        elif (dt == "Flag") or (dt == "Numeric"):
            val = self.__map(MutUtils.str2bool, vals)
            val = val[0]
        elif dt == 'String':
            try:
                val = self.__map(str, vals)
            except IndexError:
                val = True
        try:
            if (num == 1) and (not dt in ("Flag",)):
                val = val[0]
        except KeyError:
            pass

        return val

    def renderMutations(self, mutations, metadata=[], comments=[]):
        ''' Generate a simple tsv file based on the incoming mutations.
        Assumes that all mutations have the same annotations, even if some are not populated.
        Returns a file name. '''
        
        # Parse the config file
        for sectionKey in ["INFO","FORMAT","OTHER"]: # required sections for annotations
            if ConfigUtils.hasSectionKey(configParser=self.config, sectionKey=sectionKey):
                self.configTable[sectionKey] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser=self.config, sectionKey=sectionKey)
            else:
                raise ConfigInputIncompleteException("Missing %s section in the output config file." % sectionKey)

        for sectionKey in ["INFO_DESCRIPTION","FILTER_DESCRIPTION","FORMAT_DESCRIPTION"]: # optional sections for specifying description
            if ConfigUtils.hasSectionKey(configParser=self.config, sectionKey=sectionKey):
                self.configTable[sectionKey] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser=self.config, sectionKey=sectionKey)
                for ID in self.configTable[sectionKey]:
                    self.configTable[sectionKey][ID] = string.join(words=self.configTable[sectionKey][ID], sep=",")

        self.configTable["NOT_SPLIT_TAGS"] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(
            configParser=self.config, sectionKey="NOT_SPLIT_TAGS")

        path = os.getcwd()

        self.logger.info("Rendering VCF output file: " + self._filename)
        self.logger.info("Data sources included: " + str(self._datasources))
        self.logger.info("Render starting...")

        temp = tempfile.NamedTemporaryFile(dir=path) # create a temporary file to write tab-separated file
        sampleNames, chroms, annotationTable = self.__writeMutationsToFile(filename=temp.name, mutations=mutations,
                                                                           metadata=metadata)
        headers = self.__createMetaInformationHeader(metadata=metadata, comments=comments,
                                                     annotationNames=annotationTable.keys())
        header = headers.pop(len(headers)-1)
        # add sampleNames to the header
        if len(sampleNames) > 0:
            sampleNames = list(sampleNames) # convert to list to preserve ordering
            sampleNames.sort() # lexicographic ordering of sampleNames
            header = string.join([header] + sampleNames, self.delimiter)

        infos = self.__getAnnotationNamesByField(tbl=annotationTable, fld="INFO")
        fmts = self.__getAnnotationNamesByField(tbl=annotationTable, fld="FORMAT")
        filts = self.__getAnnotationNamesByField(tbl=annotationTable, fld="FILTER")
        ids = self.__getAnnotationNamesByField(tbl=annotationTable, fld="ID")
        quals = self.__getAnnotationNamesByField(tbl=annotationTable, fld="QUAL")

        # detect genotype (or GT information)
        for i in xrange(len(fmts)):
            annotation = annotationTable[fmts[i]]
            if annotation.ID == "GT":
                tmp = fmts[0]
                fmts[0] = fmts[i]
                fmts[i] = tmp
                break

        fp = file(self._filename,'w')
        fp.write(string.join([string.join(headers,'\n##'), header], "\n#"))
        fp.close()

        temp_vcf_reader = vcf.Reader(filename=self._filename, strict_whitespace=True)

        fp = file(self._filename,'w')
        vcf_writer = vcf.Writer(fp, temp_vcf_reader)

        with open(temp.name, "r") as fp:
            chrom = None
            pos = None
            ref = None
            alts = []
            ID = None # default: None
            qual = "."
            filt = []
            info = collections.OrderedDict()
            fmt = dict()
            sampleIndexes = dict([(x,i) for (i,x) in enumerate(sampleNames)])
            index = 0

            nfmts = len(fmts)
            dr = csv.DictReader(fp, delimiter=self.delimiter)
            for row in dr:
                # new (chromosome, position) pair is encountered; write pair and re-initialize data structures
                if (chrom != row["chr"]) or (pos != int(row["start"])):
                    # write existing data
                    if (not chrom is None) and (not pos is None):
                        record = self.__getRecord(chrom=chrom, pos=pos, ID=ID, ref=ref, alts=alts, qual=qual, filt=filt,
                                                  info=info, infos=infos, fmt=fmt, fmts=fmts,
                                                  sampleIndexes=sampleIndexes, sampleNames=sampleNames,
                                                  annotationTable=annotationTable)
                        vcf_writer.write_record(record)

                        index += 1
                        if index % 1000 == 0:
                            vcf_writer.flush()

                    # chromosome
                    chrom = row["chr"]
                    # position
                    pos = int(row["start"])
                    # reference allele
                    ref = row["ref_allele"]

                    # alternate allele
                    alts = [row["alt_allele"]]

                    # semi-colon separated list of unique identifiers
                    ID = None  # default: None
                    if len(ids) > 0:
                        ID = []
                        for i in xrange(len(ids)):
                            ID.append(row[ids[i]])

                    # phred-scaled quality score for the assertion made in ALT
                    qual = None  # default: None
                    if len(quals) == 1:
                        val = row[quals[0]]
                        if val.isdigit():
                            qual = int(val)

                    filt = []  # parse PASS or list of the names of filter that have failed
                    for i in xrange(len(filts)):
                        if (row[filts[i]] != "PASS") and (row[filts[i]] != "."):
                            annotation = annotationTable[filts[i]]
                            filt.append(annotation.ID)

                    info = collections.OrderedDict()
                    for i in xrange(len(infos)):  # only parsed when a new alternate arrives
                        annotation = annotationTable[infos[i]]
                        info[annotation.ID] = row[infos[i]]

                    fmt = collections.OrderedDict()
                    sampleName = row["sampleName"]
                    fmt[sampleName] = [None]*nfmts
                    for i in xrange(nfmts):
                        fmt[sampleName][i] = row[fmts[i]]
                else:
                    for i in xrange(len(ids)):
                        if not row[ids[i]] in ID:
                            ID.append(row[ids[i]])

                    if len(quals) == 1:
                        val = row[quals[0]]
                        if val.isdigit():
                            if qual != int(val):
                                raise Exception()

                    for i in xrange(len(filts)):
                        if (row[filts[i]] != "PASS") and (row[filts[i]] != "."):
                            annotation = annotationTable[filts[i]]
                            if not annotation.ID in filt:
                                filt.append(annotation.ID)

                    if row["alt_allele"] not in alts:
                        alts.append(row["alt_allele"])
                        for i in xrange(len(infos)):  # only parsed when a new alternate allele is detected
                            annotation = annotationTable[infos[i]]
                            if annotation.ID in self.configTable["SPLIT_TAGS"]["INFO"]:
                                info[annotation.ID] += "," + row[infos[i]]

                    sampleName = row["sampleName"]
                    if sampleName not in fmt:
                        fmt[sampleName] = [None]*nfmts
                        for i in xrange(nfmts):
                            fmt[sampleName][i] = row[fmts[i]]
                    else:
                        for i in xrange(nfmts):
                            annotation = annotationTable[fmts[i]]
                            if annotation.ID in self.configTable["SPLIT_TAGS"]["FORMAT"]:
                                fmt[sampleName][i] += "," + row[fmts[i]]

        record = self.__getRecord(chrom=chrom, pos=pos, ID=ID, ref=ref, alts=alts, qual=qual, filt=filt, info=info,
                                  infos=infos, fmt=fmt, fmts=fmts, sampleIndexes=sampleIndexes, sampleNames=sampleNames,
                                  annotationTable=annotationTable)
        vcf_writer.write_record(record)
        vcf_writer.close()
        return self._filename

    def __getRecord(self, chrom, pos, ID, ref, alts, qual, filt, info, infos, fmt, fmts, sampleIndexes, sampleNames,
                    annotationTable):
        if not ID is None:
            ID = string.join(ID, ";")

        if qual is None:
            qual = "."
            self.logger.warn("")

        nalts = len(alts)
        nfmts = len(fmts)

        samples = [None]*len(sampleNames)
        sampflds = [None]*nfmts # field names for the sample
        samptypes = [None]*nfmts
        sampnums = [None]*nfmts
        for sampleName in sampleNames:
            if not sampleName in fmt:
                fmt[sampleName] = [None]*nfmts
                for i in xrange(nfmts):
                    annotation = annotationTable[fmts[i]]
                    if (annotation.num == ".") and (annotation.ID in self.configTable["SPLIT_TAGS"]["FORMAT"]):
                        fmt[sampleName][i] = string.join(["."]*nalts, ",")
                    elif annotation.num == ".":
                        fmt[sampleName][i] = "."
                    else:
                        if annotation.num == 1:
                            fmt[sampleName][i] = "."
                        else:
                            fmt[sampleName][i] = string.join(["."]*annotation.num, ",")

            for i in xrange(nfmts):
                annotation = annotationTable[fmts[i]]
                sampflds[i] = annotation.ID
                sampnums[i] = annotation.num
                samptypes[i] = annotation.type
                fmt[sampleName][i] = self.__parseValue(val=fmt[sampleName][i], dt=annotation.type, num=annotation.num)

            calldata = vcf.model.make_calldata_tuple(sampflds)
            calldata._types = samptypes
            calldata._nums = sampnums
            samples[sampleIndexes[sampleName]] = calldata(*fmt[sampleName])

        for i in xrange(len(infos)):
            annotation = annotationTable[infos[i]]
            info[annotation.ID] = self.__parseValue(val=info[annotation.ID], dt=annotation.type, num=annotation.num)

        record = vcf.model._Record(chrom, pos, ID, ref, alts, qual, filt, info, string.join(sampflds, ":"),
                                   sampleIndexes)

        for sampleName in sampleNames:
            calldata = samples[sampleIndexes[sampleName]]
            samples[sampleIndexes[sampleName]] = vcf.model._Call(record, sampleName, calldata)

        record.samples = samples

        return record

    def __createChromHashCodeTable(self, chroms):
        chromHashCodeTable = dict()
        highestHashCode = 0
        for chrom in chroms:
            chromHashCodeTable[chrom] = None
            if chrom.isdigit():
                chromHashCodeTable[chrom] = int(chrom)
                if highestHashCode < chromHashCodeTable[chrom]:
                    highestHashCode = chromHashCodeTable[chrom]
        index = 0
        for chrom in sorted(chroms):
            if chromHashCodeTable[chrom] is None:
                if chrom.upper() == 'X':  # X chromosome
                    chromHashCodeTable[chrom] = highestHashCode + 1
                elif chrom.upper() == 'Y':  # Y chromosome
                    chromHashCodeTable[chrom] = highestHashCode + 2
                elif (chrom.upper() == 'M') or (chrom.upper() == 'MT'):  # mitochondrial chromosome
                    chromHashCodeTable[chrom] = highestHashCode + 3
                else:
                    index += 1
                    chromHashCodeTable[chrom] = highestHashCode + index + 3
        return chromHashCodeTable
