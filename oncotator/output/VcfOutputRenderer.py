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
from oncotator.utils.MutUtils import MutUtils
from oncotator.input.ConfigInputIncompleteException import ConfigInputIncompleteException
from oncotator.utils.TsvFileSorter import TsvFileSorter
from oncotator.utils.GenericTsvReader import GenericTsvReader


class OutputDataManager:
    def __init__(self, configTable, comments=[], md=[], mut=None, sampleNames=[]):
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.comments = comments
        self.table = configTable
        self.metadata = md
        self.mutation = mut
        self.sampleNames = sampleNames
        self.sampleNames.sort()
        # TODO: convert to annotation class, ft can be done via lookup to reverseAnnotation table
        self._Annotation = collections.namedtuple(typename="Annotation", field_names=["ft", "ID", "num", "dt", "desc",
                                                                                      "src", "split"])
        # TODO: convert to a annotation table class
        self.annotationTable, self.reverseAnnotationTable = self._createTables(self.metadata, self.mutation)

    def getHeader(self):
        return self._createHeader(self.comments, self.delimiter, self.lineterminator)

    def getFieldType(self, name):
        ft = "FORMAT"
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            ft = annotation.ft
        return ft

    def getFieldID(self, name):
        ID = name
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            ID = annotation.ID
        return ID

    def getFieldDataType(self, name):
        dataType = "String"
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            dataType = annotation.dt
        return dataType

    def getFieldNum(self, name):
        num = "."
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            num = annotation.num
        return num

    def getFieldDesc(self, name):
        desc = "Unknown"
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            desc = annotation.desc
        return desc

    def getFieldDatasource(self, name):
        src = "INPUT"
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            src = annotation.src
        return src

    def isFieldSplit(self, name):
        isSplit = False
        if name in self.annotationTable:
            annotation = self.annotationTable[name]
            isSplit = annotation.split
        return isSplit

    def getAnnotationNames(self, fieldType):
        name = []
        if fieldType in self.reverseAnnotationTable:
            name = self.reverseAnnotationTable[fieldType]
        return name

    def _createHeader(self, comments=[], delimiter="\t", lineterminator="\n"):
        headers = ["##fileformat=VCFv4.1"]
        if (comments is not None) and (len(comments) > 0):
            for i in xrange(len(comments)):
                comment = comments[i]
                if comment.startswith("fileformat=VCFv4."):
                    comments.pop(i)
                    break
        # Last line of the comments ("Oncotator v1.0.0.0rc20|") is NOT included in the header
        headers += [string.join(["##", comment], "") for comment in comments[0:len(comments)-1]]

        annotations = self.annotationTable.values()
        for annotation in annotations:
            headers += [self._annotation2str(annotation.ft, annotation.ID, annotation.desc, annotation.dt,
                                             annotation.num)]
        headers += [string.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] +
                                self.sampleNames, delimiter)]
        header = string.join(filter(None, headers), lineterminator)
        return header

    def _createTables(self, md, mut):
        """
        TODO: comments
        """
        if mut is not None:
            names = set(md.keys())
            names = names.union(mut.keys())
            names = names.difference(['chr', 'start', 'end', 'ref_allele', 'alt_allele', 'altAlleleSeen', 'sampleName',
                                      'build'])
        else:
            names = md.keys()

        table = dict()
        revTable = dict()
        for name in names:
            annotation = None
            num = "."
            tags = []
            dt = "String"
            desc = "Unknown"
            src = "INPUT"

            if name in md:
                annotation = md[name]
            elif name in mut:
                annotation = mut.getAnnotation(name)

            if annotation is not None:
                num = annotation.getNumber()
                tags = annotation.getTags()
                dt = annotation.getDataType()
                desc = annotation.getDescription()
                src = annotation.getDatasource()

            ft = self._resolveFieldType(name, tags)
            ID = self._resolveFieldID(ft, name)
            dt = self._resolveFieldDataType(ft, ID, dt)
            desc = self._resolveFieldDescription(ft, ID, desc)

            isSplitTag = False
            if ft in ("INFO", "FORMAT",):
                isSplitTag = self._determineIsSplit(ID, num, ft, "SPLIT" in annotation.getTags())
            table[name] = self._Annotation(ft, ID, num, dt, desc, src, isSplitTag)
            if ft not in revTable:
                revTable[ft] = [name]
            else:
                revTable[ft] += [name]

        return table, revTable

    def _determineIsSplit(self, ID, num, fieldType, split):
        if num in (-2,):  # by the number of samples
            split = False
            if fieldType in ("FORMAT",):
                if ID in self.table["SPLIT_TAGS"][fieldType]:  # override the default using the config file
                    split = True
        elif num in (-1,):  # by the number of alternates
            split = True
            if ID in self.table["NOT_SPLIT_TAGS"][fieldType]:  # override the default using the config file
                split = False
        elif num in (".",):
            split = False
            if ID in self.table["SPLIT_TAGS"][fieldType]:  # override the default using the config file
                split = True
        else:
            split = False

        return split

    def _resolveFieldType(self, name, tags):
        ft = None

        if ft is None:
            m = {"aggregate": "INFO", "variant": "FORMAT", "filter": "FILTER", "identifier": "ID", "quality": "QUAL"}
            m_keys = m.keys()

            for tag in tags:
                if tags in m_keys:
                    ft = m[tag]

        if ft is None:
            if name in self.table["INFO"]:
                ft = "INFO"
            elif name in self.table["FORMAT"]:
                ft = "FORMAT"
            elif name in self.table["OTHER"]:
                ft = self.table["OTHER"][name]

        if ft is None:
            if name.upper() in vcf.parser.RESERVED_INFO:
                ft = "INFO"
            elif name.upper() in vcf.parser.RESERVED_FORMAT:
                ft = "FORMAT"

        if ft is None:
            ft = "INFO"

        return ft

    def _resolveFieldID(self, ft, name):
        ID = None
        if (ft in self.table) and (name in self.table[ft]):
            ID = self.table[ft][name]
        if ID is None:
            ID = name
        return ID

    def _resolveFieldDataType(self, ft, ID, dt):
        if dt in ("String", ".", "",):
            if ft == "FILTER":
                if ID in vcf.parser.RESERVED_INFO:
                    dt = vcf.parser.RESERVED_INFO[ID]
            elif ft == "FORMAT":
                if ID in vcf.parser.RESERVED_FORMAT:
                    dt = vcf.parser.RESERVED_FORMAT[ID]
            else:
                dt = "String"
        return dt

    def _resolveFieldDescription(self, ft, ID, desc):
        # description in the config file overwrite the description in the mutation
        if (ft == "FILTER") and (ID in self.table["FILTER_DESCRIPTION"]):
            desc = self.table["FILTER_DESCRIPTION"][ID]
        elif (ft == "FORMAT") and (ID in self.table["FORMAT_DESCRIPTION"]):
            desc = self.table["FORMAT_DESCRIPTION"][ID]
        elif (ft == "INFO") and (ID in self.table["INFO_DESCRIPTION"]):
            desc = self.table["INFO_DESCRIPTION"][ID]
        if (desc is None) or (desc == ""):
            desc = "Unknown"
        return desc

    def _annotation2str(self, ft, ID, desc="Unknown", dt="String", num="."):
        if ft in ("FORMAT", "INFO",):
            return "##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (ft, ID, num, dt, desc)
        elif ft in ("FILTER",):
            return "##%s=<ID=%s,Description=\"%s\">" % (ft, ID, desc)
        return ""


class RecordFactory:
    def __init__(self, chrom, pos, ref, sampleNames):
        self._logger = logging.getLogger(__name__)
        self._chrom = chrom
        self._pos = pos
        self._ID = None  # semi-colon separated list of unique identifiers where available
        self._refAllele = ref  # reference base(s)
        self._alts = []  # comma separated list of alternate non-reference alleles called on at least one of the samples
        self._qual = None  # phred-scaled quality score for the assertion made in ALT
        self._filt = []  # PASS if this position has passed all filters
        self._info = collections.OrderedDict()  # additional information
        self._infoFieldProperty = collections.OrderedDict()
        self._fmt = [None]*len(sampleNames)
        self._fmtIDs = []
        self._fmtFieldProperty = collections.OrderedDict()
        self._sampleNames = sampleNames
        self._sampleNameIndexes = dict([(x, i) for (i, x) in enumerate(sampleNames)])
        self._fieldProperty = collections.namedtuple(typename="Property", field_names=["num", "dt", "split"])

    def getAlts(self):
        return self._alts

    def _map(self, func, iterable, bad="."):
        return [func(v) if v != bad else None for v in iterable]

    def _replace_chrs(self, text, skipList=[]):
        if not isinstance(text, str):
            return text

        dic = {",": "|", "=": "~", ";": "|", "\n": "#", "\t": "_", " ": "_", ":": ">"}
        for i, j in dic.iteritems():
            if i not in skipList:
                text = text.replace(i, j)

        return text

    def _determineVal(self, val, dt, num):
        vals = ["." if (not v or v == "None") else v for v in val]

        if dt == "Integer":
            try:
                val = self._map(int, vals)
            except ValueError:
                val = self._map(float, vals)
        elif (dt == "Float") or (dt == "Numeric"):
            val = self._map(float, vals)
        elif dt == "Flag":
            val = self._map(MutUtils.str2bool, vals)
            val = val[0]
        elif dt == "String":
            try:
                val = self._map(str, vals)
            except IndexError:
                val = True
        try:
            if (num == 1) and (not dt in ("Flag",)):
                val = val[0]
        except KeyError:
            pass
        return val

    def _resolveInfo(self):
        nalts = len(self._alts)
        nsamples = len(self._sampleNames)

        IDs = self._info.keys()
        info = collections.OrderedDict()
        for ID in IDs:
            prop = self._infoFieldProperty[ID]
            val = self._info[ID]
            num = prop.num
            dataType = prop.dt
            isSplit = prop.split

            if num == -2:  # num is the number of samples
                if nsamples < len(val):
                    tmp = [None]*nsamples
                    for sampleName in self._sampleNames:
                        sampleNameIndex = self._sampleNameIndexes[sampleName]
                        tmp[sampleNameIndex] = val[sampleNameIndex]
                    val = tmp
                elif nsamples > len(val):
                    tmp = [None]*nsamples
                    for sampleName in self._sampleNames:
                        sampleNameIndex = self._sampleNameIndexes[sampleName]
                        if sampleNameIndex < len(val):
                            tmp[sampleNameIndex] = val[sampleNameIndex]
                    val = tmp
            elif num == -1:  # num is the number of alternative alleles
                if nalts < len(val):
                    val = val[0:nalts]
                elif nalts > len(val):
                    val += [None]*(nalts-len(val))
            elif num == 0:
                val = val[0]
            elif num == ".":  # num is unknown
                if isSplit:  # now, num is the number of alternative alleles
                    if nalts < len(val):
                        val = val[0:nalts]
                    elif nalts > len(val):
                        val += [None]*(nalts-len(val))
            else:
                if num < len(val):
                    val = val[0:num]
                elif num > len(val):
                    val += [None]*(num-len(val))

            val = map(str, val)
            val = self._determineVal(val, dataType, num)

            # TODO: Remove isinstance calls here and for the FORMAT field
            if isinstance(val, bool):
                if val is True:
                    info[ID] = val
            elif isinstance(val, list):  # all values are None
                if len(filter(None, val)) != 0:
                    info[ID] = val
            elif val is None:  # value is None
                pass
            else:
                info[ID] = val

        return info

    def _resolveSamples(self, record):
        nalts = len(self._alts)
        samples = [None]*len(self._sampleNames)

        IDs = [None]*len(self._fmtIDs)
        dataTypes = [None]*len(self._fmtFieldProperty)
        nums = [None]*len(self._fmtFieldProperty)

        for sampleName in self._sampleNames:
            sampleData = [None]*len(self._fmtIDs)
            sampleNameIndex = self._sampleNameIndexes[sampleName]
            data = self._fmt[sampleNameIndex]
            for i in xrange(len(self._fmtIDs)):
                ID = self._fmtIDs[i]
                IDs[i] = self._fmtIDs[i]
                prop = self._fmtFieldProperty[ID]
                dataTypes[i] = prop.dt
                nums[i] = prop.num
                split = prop.split
                val = ["None"]

                if (data is not None) and (ID in data):
                    val = data[ID]
                    if nums[i] in (-2,):
                        if split:  # now, num is the number of alternative alleles
                            if nalts < len(val):
                                val = val[0:nalts]
                            elif nalts > len(val):
                                val += [None]*(nalts-len(val))
                    elif nums[i] in (-1,):
                        if nalts < len(val):
                            val = val[0:nalts]
                        elif nalts > len(val):
                            val += [None]*(nalts-len(val))
                    elif nums[i] in (0,):
                        val = val[0]
                    elif nums[i] in (".",):
                        if split:  # now, num is the number of alternative alleles
                            if nalts < len(val):
                                val = val[0:nalts]
                            elif nalts > len(val):
                                val += [None]*(nalts-len(val))
                    else:
                        if nums[i] < len(val):
                            val = val[0:nums[i]]
                        elif nums[i] > len(val):
                            val += [None]*(nums[i]-len(val))

                val = string.join(map(str, val), ",")
                val = self._determineVal(val, dataTypes[i], nums[i])
                if ID == "GT":
                    if isinstance(val, list):
                        val = val[0]
                sampleData[i] = val

            calldata = vcf.model.make_calldata_tuple(IDs)
            calldata._types = dataTypes
            calldata._nums = nums
            samples[self._sampleNameIndexes[sampleName]] = calldata(*sampleData)

        for sampleName in self._sampleNames:
            samples[self._sampleNameIndexes[sampleName]] = \
                vcf.model._Call(record, sampleName, samples[self._sampleNameIndexes[sampleName]])
        return samples

    def createRecord(self):
        chrom = self._chrom
        pos = self._pos
        refAllele = self._refAllele
        alts = self._alts
        if len(alts) == 0:
            alts = ["."]

        ID = self._ID
        if ID is not None:
            ID = string.join(ID, ";")

        qual = self._qual
        if qual is None:
            self._logger.warn("Variant at chromosome %s and position %s is missing phred-scaled quality score.")

        filt = self._filt
        info = self._resolveInfo()
        fmt = string.join(self._fmtIDs, ":")
        record = vcf.model._Record(chrom, pos, ID, refAllele, alts, qual, filt, info, fmt, self._sampleNameIndexes)
        record.samples = self._resolveSamples(record)
        return record

    def addFormat(self, sampleName, ID, num=".", dt="String", val=None, src="INPUT", split=True):
        if sampleName in self._sampleNames:
            sampleNameIndex = self._sampleNameIndexes[sampleName]
            if self._fmt[sampleNameIndex] is None:
                self._fmt[sampleNameIndex] = collections.OrderedDict()
                self._fmt[sampleNameIndex]["GT"] = None
                self._fmtIDs = ["GT"]
                self._fmtFieldProperty["GT"] = self._fieldProperty(1, "String", False)

            if src not in ("INPUT",):
                val = self._replace_chrs(val)
            else:
                val = self._replace_chrs(val, [","])

            if num in (-2,):  # num is the number of samples
                if split:
                    if ID not in self._fmt[sampleNameIndex]:
                        self._fmt[sampleNameIndex][ID] = [val]
                    else:
                        self._fmt[sampleNameIndex][ID] += [val]
                else:
                    self._fmt[sampleNameIndex][ID] = [val]
            elif num in (-1,):  # num is the number of alternate alleles
                if split:
                    nalts = len(self._alts)
                    if nalts == 1:
                        self._fmt[sampleNameIndex][ID] = [val]
                    elif nalts > 1:
                        vals = self._fmt[sampleNameIndex][ID]
                        if len(vals) < nalts:
                            self._fmt[sampleNameIndex][ID] += [val]
                else:
                    self._fmt[sampleNameIndex][ID] = val.split(",")
            elif num in (".",):  # num is unknown
                if split:  # now, num is the number of alternate alleles
                    nalts = len(self._alts)
                    if nalts == 1:
                        self._fmt[sampleNameIndex][ID] = [val]
                    elif nalts > 1:
                        vals = self._fmt[sampleNameIndex][ID]
                        if len(vals) < nalts:
                            self._fmt[sampleNameIndex][ID] += [val]
                else:
                    self._fmt[sampleNameIndex][ID] = val.split(",")
            else:
                if split:
                    if ID not in self._fmt[sampleNameIndex]:
                        self._fmt[sampleNameIndex][ID] = [val]
                    else:
                        self._fmt[sampleNameIndex][ID] += [val]
                else:
                    self._fmt[sampleNameIndex][ID] = val.split(",")

            if ID not in self._fmtIDs:
                self._fmtIDs += [ID]

            if ID not in self._fmtFieldProperty:
                self._fmtFieldProperty[ID] = self._fieldProperty(num, dt, split)

    def addInfo(self, sampleName, ID, num=".", dt="String", val=None, src="INPUT", split=True):
        if src not in ("INPUT",):
            val = self._replace_chrs(val, [":"])
        else:
            val = self._replace_chrs(val, [":", ","])

        if num in (-2,):  # num is the number of samples
            nsamples = len(self._sampleNames)
            if sampleName in self._sampleNames:
                if ID not in self._info:
                    self._info[ID] = [None]*nsamples
                sampleNameIndex = self._sampleNameIndexes[sampleName]
                self._info[ID][sampleNameIndex] = val
        elif num in (-1,):  # num is the number of alternate alleles
            if split:
                nalts = len(self._alts)
                if nalts == 1:
                    self._info[ID] = [val]
                elif nalts > 1:
                    vals = self._info[ID]
                    if len(vals) < nalts:
                        self._info[ID] += [val]
            else:
                self._info[ID] = val.split(",")
        elif num in (".",):  # num is unknown
            if split:  # now, num is the number of alternate alleles
                nalts = len(self._alts)
                if nalts == 1:
                    self._info[ID] = [val]
                elif nalts > 1:
                    vals = self._info[ID]
                    if len(vals) < nalts:
                        self._info[ID] += [val]
            else:
                self._info[ID] = val.split(",")
        else:
            if split:
                if ID not in self._info:
                    self._info[ID] = [val]
                else:
                    self._info[ID] += [val]
            else:
                self._info[ID] = val.split(",")

        if (ID in self._info) and (ID not in self._infoFieldProperty):
            self._infoFieldProperty[ID] = self._fieldProperty(num, dt, split)

    def addQual(self, qual):
        try:
            self._qual = int(qual)
        except ValueError:
            try:
                self._qual = float(qual)
            except ValueError:
                self._qual = None

    def addID(self, ID):
        if ID not in (".", "",):
            if self._ID is None:
                self._ID = [ID]
            else:
                if ID not in self._ID:
                    self._ID += [ID]

    def addAlt(self, alt):
        if alt in ("",):
            alt = "."
        if alt not in self._alts:
            self._alts += [alt]

    def addFilter(self, filt, val):
        if val not in ("PASS", ".",):
            if filt not in self._filt:
                self._filt += [filt]

    def setChrom(self, chrom):
        self._chrom = chrom

    def setPos(self, pos):
        self._pos = pos

    def setReferenceAllele(self, ref):
        self._refAllele = ref

    def setSampleNames(self, sampleNames):
        self._sampleNames = sampleNames
        self._sampleNameIndex = dict([(x, i) for (i, x) in enumerate(sampleNames)])


class VcfOutputRenderer(OutputRenderer):
    """
    The SimpleOutputRenderer renders a basic tsv file from the given mutations.  All annotations are included with real names as column headers.

    Header is determined by the first mutation given.

    No attention is paid to order of the headers.
    """
    _vcfAnnotation = collections.namedtuple(typename="Annotation", field_names=["field", "ID", "num", "type"])

    def __init__(self, filename, datasources=[], configFile='vcf.out.config'):
        """
        Constructor
        """
        self._filename = filename
        self.logger = logging.getLogger(__name__)
        self._datasources = datasources
        self.config = ConfigUtils.createConfigParser(configFile, ignoreCase=False)
        self.chromHashCodeTable = None  # maps every chromosome in the mutations to a sortable integer
        self.configTable = dict()
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.sampleNames = []  # all sample names in the mutations
        self.chroms = []  # all chromosomes in the mutations
        self.reservedAnnotationNames = ['chr', 'start', 'ref_allele', 'alt_allele', 'end']

    def _writeMuts2Tsv(self, filename, fieldnames, muts):
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
        fieldnames = self.reservedAnnotationNames
        if mut is not None:
            fieldnames = set(fieldnames).union(md.keys())
            fieldnames = fieldnames.union(mut.keys())
            fieldnames = fieldnames.difference(["end", "sampleName", "build"])
            if "sampleName" in mut:
                fieldnames = fieldnames.union(["sampleName"])
        return list(fieldnames)

    def _doFieldsExist(self, sect, fields):
        tbl = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, sect)
        ks = set(tbl.keys())
        for fld in fields:
            if fld not in ks:
                raise ConfigInputIncompleteException("Missing %s field in %s section in the output config file."
                                                     % (fld, sect))

    def _parseConfig(self):
        sects = ["INFO", "FORMAT", "OTHER", "NOT_SPLIT_TAGS", "INFO_DESCRIPTION", "FILTER_DESCRIPTION",
                 "FORMAT_DESCRIPTION", "SPLIT_TAGS"]
        for sect in sects:
            if not ConfigUtils.hasSectionKey(self.config, sect):
                raise ConfigInputIncompleteException("Missing %s section in the output config file." % sect)
            if sect in ("OTHER",):
                reqKs = ["ID", "QUAL", "FILTER"]
                self._doFieldsExist(sect, reqKs)
            elif sect in ("NOT_SPLIT_TAGS", "SPLIT_TAGS",):
                reqKs = ["INFO", "FORMAT"]
                self._doFieldsExist(sect, reqKs)

        table = dict()
        for sect in sects:
            if sect in ("INFO", "FORMAT", "OTHER",):
                table[sect] = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(self.config, sect)
            elif sect in ("INFO_DESCRIPTION", "FILTER_DESCRIPTION", "FORMAT_DESCRIPTION",):
                table[sect] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, sect)
                for k, v in table[sect].items():
                    table[sect][k] = string.join(v, ",")
            elif sect in ("SPLIT_TAGS", "NOT_SPLIT_TAGS",):
                table[sect] = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config, sect)

        return table

    def renderMutations(self, mutations, metadata=[], comments=[]):
        """ Generate a simple tsv file based on the incoming mutations.
        Assumes that all mutations have the same annotations, even if some are not populated.
        Returns a file name. """
        self.logger.info("Rendering VCF output file: " + self._filename)

        # Initialize config table
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
        isNew = False
        if curChrom != prevChrom:
            isNew = True
        if curPos != prevPos:
            isNew = True
        return isNew

    def _renderSortedTsv(self, vcfFilename, tsvFilename, sampleNames, dataManager):
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
        IDs = dataManager.getAnnotationNames("ID")
        quals = dataManager.getAnnotationNames("QUAL")
        filts = dataManager.getAnnotationNames("FILTER")
        infos = dataManager.getAnnotationNames("INFO")
        formats = dataManager.getAnnotationNames("FORMAT")

        altAllele = m["alt_allele"]
        sampleName = m["sampleName"] if "sampleName" in m else None
        recordFactory.addAlt(altAllele)

        for name in IDs:
            val = m[name] if name in m else ""
            recordFactory.addID(val)

        qual = quals[0]
        recordFactory.addQual(m[qual])

        for name in filts:
            ID = dataManager.getFieldID(name)
            val = m[name] if name in m else ""
            recordFactory.addFilter(ID, val)

        for name in infos:
            ID = dataManager.getFieldID(name)
            num = dataManager.getFieldNum(name)
            dt = dataManager.getFieldDataType(name)
            src = dataManager.getFieldDatasource(name)
            isSplit = dataManager.isFieldSplit(name)
            val = m[name] if name in m else ""
            recordFactory.addInfo(sampleName, ID, num, dt, val, src, isSplit)

        for name in formats:
            ID = dataManager.getFieldID(name)
            num = dataManager.getFieldNum(name)
            dt = dataManager.getFieldDataType(name)
            src = dataManager.getFieldDatasource(name)
            isSplit = dataManager.isFieldSplit(name)
            val = m[name] if name in m else ""
            recordFactory.addFormat(sampleName, ID, num, dt, val, src, isSplit)

        return recordFactory

    def _createChrom2HashCodeTable(self, chroms):
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