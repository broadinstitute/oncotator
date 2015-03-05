"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""

import collections
from collections import OrderedDict
import logging
import string
import vcf
from oncotator.utils.MutUtils import MutUtils


class RecordBuilder:

    fieldProperty = collections.namedtuple(typename="Property", field_names=["num", "dataType", "isSplit"])
    logger = logging.getLogger(__name__)

    def __init__(self, chrom, pos, ref, sampleNames):
        """

        :param chrom:
        :param pos:
        :param ref:
        :param sampleNames:
        """
        self._chrom = chrom
        self._pos = pos
        self._ID = None  # semi-colon separated list of unique identifiers where available
        self._refAllele = ref  # reference base(s)
        self._alts = []  # comma separated list of alternate non-reference alleles called on at least one of the samples
        self._qual = None  # phred-scaled quality score for the assertion made in ALT
        self._filt = []  # PASS if this position has passed all filters
        self._info = collections.OrderedDict()  # additional information
        self._infoFieldProperty = collections.OrderedDict()
        self._fmt = [OrderedDict() for x in xrange(0,len(sampleNames))]
        self._fmtIDs = []
        self._fmtFieldProperty = collections.OrderedDict()
        self._sampleNames = sampleNames
        self._sampleNameIndexes = dict([(x, i) for (i, x) in enumerate(sampleNames)])

        if sampleNames is not None and len(sampleNames) != 0:
            self._fmtIDs = ["GT"]
            self._fmtFieldProperty["GT"] = self.fieldProperty(1, "String", False)

    def _map(self, func, iterable, bad=(".", "",)):
        """

        :param func:
        :param iterable:
        :param bad:
        :return:
        """
        return [func(v) if v not in bad else None for v in iterable]

    def _resolveInfo(self):
        """


        :return:
        """
        IDs = self._info.keys()
        info = collections.OrderedDict()

        nalts = len(self._alts)
        nsamples = len(self._sampleNames)
        for ID in IDs:
            val = self._info[ID]
            prop = self._infoFieldProperty[ID]

            if prop.num == -2:
                if len(val) == nsamples and len(filter(None, val)) != 0:
                    info[ID] = val
            elif prop.num == -1:
                if len(val) == nalts and len(filter(None, val)) != 0:
                    info[ID] = val
            elif prop.num == 0:
                if val:
                    info[ID] = val
            elif prop.num is None:
                if prop.isSplit:
                    if len(val) == nalts and len(filter(None, val)) != 0:
                        info[ID] = val
                elif len(filter(None, val)) != 0:
                    info[ID] = val
            else:
                if len(val) == prop.num and len(filter(None, val)) != 0:
                    info[ID] = val

        return info

    def _resolveSamples(self, record):
        """

        :param record:
        :return:
        """
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
                dataTypes[i] = prop.dataType
                nums[i] = prop.num if prop.num is not None else "."
                val = [None]

                if data is not None and ID in data:
                    val = data[ID]

                if prop.num == -2:
                    pass
                elif prop.num == -1:
                    if len(val) != nalts:
                        val = nalts*[None]
                elif prop.num == 0:
                    pass
                elif prop.num is None:
                    if prop.isSplit:
                        if len(val) != nalts:
                            val = nalts*[None]
                else:
                    if len(val) != prop.num:
                        val = abs(prop.num)*[None]

                if ID == "GT":
                    sampleData[i] = val[0]
                else:
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
        """


        :return:
        """
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
            self.logger.warn("Variant at chromosome %s and position %s is missing phred-scaled quality score "
                             "(typically: annotation 'qual')." % (chrom, pos))

        filt = self._filt
        info = self._resolveInfo()
        fmt = string.join(self._fmtIDs, ":")
        record = vcf.model._Record(chrom, pos, ID, refAllele, alts, qual, filt, info, fmt, self._sampleNameIndexes)
        record.samples = self._resolveSamples(record)

        return record

    def _correct(self, iterable, bad=(".", "",)):
        """

        :param iterable:
        :param bad:
        :return:
        """
        return [v if v not in bad else None for v in iterable]

    def _fixVal(self, val, isSplit):
        """

        :param val:
        :param isSplit:
        :return:
        """
        if isSplit:
            val = MutUtils.replaceChrs(val, ",=;\n\t ", "|~|#__")  # exclude ":"
        else:
            val = MutUtils.replaceChrs(val, "=;\n\t :", "~|#__>")  # exclude ":" and ","

        if not isSplit:
            val = self._correct(val.split(","))
        else:
            val = self._correct([val])

        return val

    def _determineVal2FixedNumField(self, data, field_name, num, is_split, val):
        """

        :param data:
        :param field_name:
        :param num:
        :param is_split:
        :param val:
        """
        if field_name not in data:
            data[field_name] = self._fixVal(val, is_split)
        elif is_split and num > 1:
            vals = data[field_name]
            if len(vals) < num:
                vals += self._fixVal(val, is_split)
                data[field_name] = vals

    def addInfo(self, sampleName, ID, num=None, dataType="String", val=None, isSplit=True):
        """

        :param sampleName:
        :param ID:
        :param num:
        :param dataType:
        :param val:
        :param isSplit:
        """
        if num == -2:  # num is the number of samples
            nsamples = len(self._sampleNames)
            if sampleName in self._sampleNames:
                if ID not in self._info:
                    self._info[ID] = [None]*nsamples
                sampleNameIndex = self._sampleNameIndexes[sampleName]
                val = self._fixVal(val, isSplit)
                self._info[ID][sampleNameIndex] = val[0]
        elif num == -1:  # num is the number of alternate alleles
            nalts = len(self._alts)
            self._determineVal2FixedNumField(self._info, ID, nalts, isSplit, val)
        elif num == 0:  # num is either true or false
            if ID not in self._info:

                if dataType == 'Flag' and val == '':
                    # This prevents None from being returned by self._map call
                    val = 'false'

                val = self._map(MutUtils.str2bool, self._fixVal(val, isSplit))  # convert the value to boolean
                self._info[ID] = val[0]
        elif num is None:  # num is unknown
            nalts = len(self._alts)
            self._determineVal2FixedNumField(self._info, ID, nalts, isSplit, val)
        else:
            self._determineVal2FixedNumField(self._info, ID, num, isSplit, val)

        if ID not in self._infoFieldProperty:
            self._infoFieldProperty[ID] = self.fieldProperty(num, dataType, isSplit)

    def _determineGenotype(self):
        nalts = len(self._alts)
        genotype = string.join(map(str, [0, nalts]), "/")  # unphased genotype
        return genotype

    def addGTField(self, sampleName, inferGenotype):
        if sampleName in self._sampleNames:  # FORMAT fields can never have a value of type flag
            sampleNameIndex = self._sampleNameIndexes[sampleName]
            if "GT" not in self._fmt[sampleNameIndex].keys():
                self._fmtFieldProperty["GT"] = self.fieldProperty(1, "String", False)
                if inferGenotype:
                    self._fmt[sampleNameIndex]["GT"] = [self._determineGenotype()]


    def addFormat(self, sampleName, field_name, num=None, dataType="String", val=None, isSplit=True):
        """

        :param sampleName:
        :param field_name:
        :param num:
        :param dataType:
        :param val:
        :param isSplit:
        :param bool inferGenotype:
        """
        if sampleName in self._sampleNames and num != 0:  # FORMAT fields can never have a value of type flag
            sampleNameIndex = self._sampleNameIndexes[sampleName]

            # If GT was specified by the input, then delete the previous copy created in the above few lines.
            if field_name == "GT":
                try:
                    del self._fmt[sampleNameIndex][field_name]
                except KeyError:
                    pass

            nalts = len(self._alts)

            if num == -2 or num == -1 or num is None:
                if not isSplit:
                    self._determineVal2FixedNumField(self._fmt[sampleNameIndex], field_name, nalts, isSplit, val)
                else:
                    # Make sure that all samples get a new entry for the alternate in this field, since it is split.
                    for sample_name in self._sampleNames:
                        if field_name not in self._fmt[self._sampleNameIndexes[sample_name]].keys():
                            self._fmt[self._sampleNameIndexes[sample_name]][field_name] = []
                        sample_field_dict = self._fmt[self._sampleNameIndexes[sample_name]]
                        while len(sample_field_dict[field_name]) < nalts:
                            sample_field_dict[field_name].append(None)
                    this_sample_field_dict = self._fmt[sampleNameIndex]
                    this_sample_field_dict[field_name][-1] = self._fixVal(val, isSplit)[0]
            else:  # num is fixed
                self._determineVal2FixedNumField(self._fmt[sampleNameIndex], field_name, num, isSplit, val)


            if field_name not in self._fmtIDs:
                self._fmtIDs += [field_name]

            if field_name not in self._fmtFieldProperty:
                self._fmtFieldProperty[field_name] = self.fieldProperty(num, dataType, isSplit)

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
        if val == "FAIL":
            if filt not in self._filt:
                self._filt += [filt]
        elif val == "":
            self._filt = None

    def setChrom(self, chrom):
        self._chrom = chrom

    def setPos(self, pos):
        self._pos = pos

    def setReferenceAllele(self, ref):
        self._refAllele = ref

    def setSampleNames(self, sampleNames):
        self._sampleNames = sampleNames
        self._sampleNameIndex = dict([(x, i) for (i, x) in enumerate(sampleNames)])
