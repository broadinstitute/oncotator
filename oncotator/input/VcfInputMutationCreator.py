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
import re
import logging
import string
import copy

import vcf

from oncotator.Metadata import Metadata
from oncotator.utils.TagConstants import TagConstants
from InputMutationCreator import InputMutationCreator
from oncotator.MutationData import MutationData
from oncotator.utils.MutUtils import MutUtils
from oncotator.Annotation import Annotation
from oncotator.config_tables.ConfigTableCreatorFactory import ConfigTableCreatorFactory


class VcfInputMutationCreator(InputMutationCreator):

    def __init__(self, filename, configFile='vcf.in.config'):
        """

        :param filename:
        :param configFile:
        """
        self.filename = filename
        self.configFilename = configFile
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)
        self.configTableBuilder = ConfigTableCreatorFactory.getConfigTableCreatorInstance("input_vcf")
        self.logger = logging.getLogger(__name__)

    def _addGenotypeDataToMutation(self, mutation, record, index):
        """


        :param mutation: input mutation object
        :param record:
        :param index:
        """
        IDs = self.configTable.getFormatFieldIDs()
        genotypeData = None

        if len(IDs) != 0:
            sampleName = mutation.getAnnotation("sample_name").getValue()
            genotypeData = record.genotype(sampleName)

        if record.FORMAT is not None:
            for ID in IDs:
                val = ""
                dataType = self.vcf_reader.formats[ID].type

                num = self.vcf_reader.formats[ID].num
                tags = [TagConstants.FORMAT]

                if (genotypeData is not None) and (hasattr(genotypeData.data, ID)):
                    isSplitTag = self._determineIsSplit(ID, num, "FORMAT")
                    if isSplitTag:
                        val = genotypeData[ID][index]
                    else:
                        val = genotypeData[ID]

                    if isSplitTag:
                        tags += [TagConstants.SPLIT]
                    else:
                        tags += [TagConstants.NOT_SPLIT]

                if (val is None) or (val == ""):
                    if dataType == "Flag":
                        val = "False"
                    else:
                        val = ""
                elif isinstance(val, list):
                    val = string.join(["" if v is None else str(v) for v in val], ",")
                else:
                    val = str(val)

                name = self.configTable.getFormatFieldName(ID)
                mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.formats[ID].desc, tags=tags,
                                          number=num)

        return mutation

    def _addInfoDataToMutation(self, mutation, record, index):
        """
        This method

        :param mutation:
        :param record:
        :param index:
        :return:
        """
        IDs = self.configTable.getInfoFieldIDs()

        for ID in IDs:
            val = ""
            dataType = self.vcf_reader.infos[ID].type

            num = self.vcf_reader.infos[ID].num
            tags = [TagConstants.INFO]

            if ID in record.INFO:
                isSplitTag = self._determineIsSplit(ID, num, "INFO")
                if isSplitTag:
                    val = record.INFO[ID][index]
                else:
                    val = record.INFO[ID]

                if isSplitTag:
                    tags += [TagConstants.SPLIT]
                else:
                    tags += [TagConstants.NOT_SPLIT]

            if (val is None) or (val == ""):
                if dataType == "Flag":
                    val = "False"
                else:
                    val = ""
            elif isinstance(val, list):
                val = string.join(["" if v is None else str(v) for v in val], ",")
            else:
                val = str(val)

            name = self.configTable.getInfoFieldName(ID)
            mutation.createAnnotation(name, val, "INPUT", dataType, self.vcf_reader.infos[ID].desc, tags=tags,
                                      number=num)

        return mutation

    def _determineIsSplit(self, ID, num, fieldType):
        """

        :param ID:
        :param num:
        :param fieldType:
        :return:
        """
        if num == -2:  # by the number of samples
            if fieldType == "FORMAT":
                if self.configTable.isFieldIDInSplitSet(fieldType, ID):
                    return True
            else:
                return False
        elif num == -1:  # by the number of alternates
            if self.configTable.isFieldIDInNotSplitSet(fieldType, ID):  # override the default using the config file
                return False
            else:
                return True
        elif num == 0:
            return False
        elif num is None:
            if self.configTable.isFieldIDInSplitSet(fieldType, ID):  # override the default using the config file
                return True
            else:
                return False

        return False

    def _determineAltSeen(self, sample_gt_str, index):
        """Look at the genotype string to see if the alternate is present.  GT of ./. is considered a 'yes'"""
        # TODO: Confirm that it is necessary to take into account the index.  Take into account the index.
        is_alt_seen = "True"
        if sample_gt_str is not None:
            # Split genotype field into number of entries as ploidy.  Using chars '/' or '|'
            gt_haploid_list = re.split('/|\|', sample_gt_str)
            if all([haploid != str(index) for haploid in gt_haploid_list]):
                is_alt_seen = "False"
        else:
            # GT is ./.
            is_alt_seen = "False"
        return is_alt_seen

    def createMutations(self):
        """ Creates a mutation for each mutation by each sample, regardless of allelic depth, etc.
            
            Annotations that this will generate (as source = "INPUT"):
                sampleName
                isCalled
                isAltAllele
                allelic_depth -- DP in the vcf

            TODO: Complete documentation
        """
        self.configTable = self.configTableBuilder.getConfigTable(filename=self.filename,
                                                                  configFilename=self.configFilename)

        for record in self.vcf_reader:
            for index in range(len(record.ALT)):
                mut = self._createMutation(record, index)

                if len(record.samples) <= 0:
                    yield mut
                else:
                    sampleRecList = record.samples
                    sample_names = [s.sample for s in sampleRecList]
                    is_tumor_normal_vcf = "NORMAL" in sample_names and len(sample_names) == 2
                    if is_tumor_normal_vcf:
                        logging.getLogger(__name__).info("Tumor-Normal VCF detected.  The Normal will assume GT= 0/0, "
                                                         "unless GT field specified otherwise.")

                    for sample in sampleRecList:
                        sampleMut = copy.deepcopy(mut)
                        sample_name = sample.sample
                        if is_tumor_normal_vcf and sample_name != "NORMAL":
                            sampleMut.createAnnotation("tumor_barcode", sample_name, "INPUT")
                        sampleMut.createAnnotation("sample_name", sample_name, "INPUT")

                        #TODO: Confirm that alt_allele_seen will be False in all cases of GT = ./.
                        genotype = "GT"
                        is_alt_seen = "True"
                        if genotype in sample.data._fields:
                            is_alt_seen = self._determineAltSeen(sample.data.GT, index + 1)

                        # HACK: If the sample name is NORMAL, there is more than one sample, and
                        # there is no GT field (or GT is ./.) then assume that this is alt_allele_seen of False
                        if is_tumor_normal_vcf and sample_name == "NORMAL" and (genotype not in sample.data._fields):
                            is_alt_seen = "False"
                        sampleMut["alt_allele_seen"] = is_alt_seen
                        sampleMut = self._addGenotypeDataToMutation(sampleMut, record, index)

                        yield sampleMut

    def _createMutation(self, record, index):
        """

        :param record:
        :param index:
        :return:
        """
        chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)

        ref = record.REF
        if ref == ".":
            ref = ""

        alt = record.ALT[index]
        if alt is None:
            alt = ""
        else:
            alt = str(alt)

        # Write end position as it would be in MAF format
        endPos = int(record.POS)
        if len(alt) < len(ref):  # deletion
            endPos += len(ref) - len(alt) - 1
        elif len(alt) > len(ref):  # insertion
            endPos += 1
        mut = MutationData(chrom, record.POS, endPos, ref, alt, "hg19")

        ID = record.ID
        if ID is None:
            ID = ""
        mut.createAnnotation("id", ID, "INPUT", tags=[TagConstants.ID])

        mut.createAnnotation("qual", str(record.QUAL), "INPUT", tags=[TagConstants.QUAL])
        for flt in self.vcf_reader.filters:  # for each filter in the header
            description = self.vcf_reader.filters[flt].desc  # parse the description
            if (len(record.FILTER) != 0) and \
                    (flt in record.FILTER):  # if the filter is mentioned for this variant, then it failed
                mut.createAnnotation(flt, "FAIL", "INPUT", annotationDescription=description,
                                     tags=[TagConstants.FILTER])
            else:
                mut.createAnnotation(flt, "PASS", "INPUT", annotationDescription=description,
                                     tags=[TagConstants.FILTER])
        mut.createAnnotation("alt_allele_seen", str(True), "INPUT")
        mut = self._addInfoDataToMutation(mut, record, index)
        return mut

    def reset(self):
        """ Resets the internal state, so that mutations can be generated. """
        self.vcf_reader = vcf.Reader(filename=self.filename, strict_whitespace=True)

    def _parseMiscellaneousMetadata(self):
        """


        :return:
        """
        comments = []
        keys = self.vcf_reader.metadata.keys()
        for key in keys:
            vals = self.vcf_reader.metadata[key]
            if isinstance(vals, list):
                for itm in vals:
                    if isinstance(itm, dict):
                        val = [string.join([str(k), str(v)], "=") for k, v in itm.iteritems()]
                        val = string.join(val, ",")
                        val = string.join(["<", ">"], val)
                    else:
                        val = str(itm)

                    comment = string.join([key, val], "=")
                    comments.append(comment)
            elif isinstance(vals, dict):
                val = [string.join([str(k), str(v)], "=") for k, v in vals.iteritems()]
                val = string.join(val, ",")
                val = string.join(["<", ">"], val)

                comment = string.join([key, val], "=")
                comments.append(comment)
            else:
                val = str(vals)
                comment = string.join([key, val], "=")
                comments.append(comment)
        return comments

    def _parseContigsMetadata(self):
        """


        :return:
        """
        comments = []
        keys = self.vcf_reader.contigs.keys()
        for key in keys:
            val = self.vcf_reader.contigs[key]
            ID = val.id
            length = str(val.length)
            val = string.join(["contig=<ID=", ID, ",length=", length, ">"], "")
            comments.append(val)
        return comments

    def _parseAltsMetadata(self):
        pass

    def getComments(self):
        """ Comments often need to be passed into the output.  Get the comments from the input file."""
        comments = self._parseMiscellaneousMetadata()
        comments += self._parseContigsMetadata()
        return comments

    def _addFormatFields2Metadata(self, metadata):
        """

        :param metadata:
        :return:
        """
        for ID in self.configTable.getFormatFieldIDs():
            name = self.configTable.getFormatFieldName(ID)
            num = self.vcf_reader.formats[ID].num
            tags = [TagConstants.FORMAT]
            isSplitTag = self._determineIsSplit(ID, num, "FORMAT")
            if isSplitTag:
                tags += [TagConstants.SPLIT]
            metadata[name] = Annotation("", "INPUT", self.vcf_reader.formats[ID].type, self.vcf_reader.formats[ID].desc,
                                        tags=tags, number=num)
        return metadata

    def _addInfoFields2Metadata(self, metadata):
        """

        :param metadata:
        :return:
        """
        for ID in self.configTable.getInfoFieldIDs():
            name = self.configTable.getInfoFieldName(ID)
            num = self.vcf_reader.infos[ID].num
            tags = [TagConstants.INFO]
            isSplitTag = self._determineIsSplit(ID, num, "INFO")
            if isSplitTag:
                tags += [TagConstants.SPLIT]
            metadata[name] = Annotation("", "INPUT", self.vcf_reader.infos[ID].type, self.vcf_reader.infos[ID].desc,
                                        tags=tags, number=num)
        return metadata

    def _addFilterFields2Metadata(self, metadata):
        """

        :param metadata:
        :return:
        """
        for filt in self.vcf_reader.filters:  # for each filter in the header
            metadata[filt] = Annotation("", "INPUT", "String", self.vcf_reader.filters[filt].desc,
                                        tags=[TagConstants.FILTER])
        return metadata

    def _createMetadata(self):
        """


        :return:
        """
        self.configTable = self.configTableBuilder.getConfigTable(filename=self.filename,
                                                                  configFilename=self.configFilename)

        metadata = Metadata()
        metadata = self._addFilterFields2Metadata(metadata)
        metadata = self._addFormatFields2Metadata(metadata)
        metaData = self._addInfoFields2Metadata(metadata)

        metaData["id"] = Annotation("", "INPUT", "String", "", [TagConstants.ID])
        metaData["qual"] = Annotation("", "INPUT", "String", "", [TagConstants.QUAL])

        return metadata

    def getMetadata(self):
        metadata = self._createMetadata()
        return metadata
