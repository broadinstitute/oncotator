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

import copy
import logging
import string
import vcf
from oncotator.utils.TagConstants import TagConstants
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.MutUtils import MutUtils
import operator


class IndexedVcfDatasource(Datasource):
    """
    A datasource derived from a VCF file.  Expects a bgzipped vcf using Tabix.

    Instructions on how to index file using Tabix prior Oncotator:

    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.


    The following VCF columns will be added to the output.
        REF, ALT and INFO


    Multiple annotation data for the same mutation will be delimited by "|".

    TODO: Support for FORMAT and sample columns in output

    """
    def __init__(self, src_file, title='', version=None, match_mode="exact"):
        # all of the info is coming from config file
        super(IndexedVcfDatasource, self).__init__(src_file, title=title, version=version)
        self.vcf_reader = vcf.Reader(filename=src_file, strict_whitespace=True)
        self.vcf_info_headers = self.vcf_reader.infos.keys()
        # Annotation name from INFO columns
        self.output_vcf_headers = dict([(ID, '_'.join([self.title, ID])) for ID in self.vcf_info_headers])
        # Type information
        self.output_vcf_types = dict([(ID, self.vcf_reader.infos[ID].type) for ID in self.vcf_info_headers])
        # Number information
        self.output_vcf_nums = dict([(ID, self.vcf_reader.infos[ID].num) for ID in self.vcf_info_headers])
        # Descriptions
        self.output_vcf_descs = dict([(ID, self.vcf_reader.infos[ID].desc) for ID in self.vcf_info_headers])
        self.match_mode = match_mode

    def _determine_tags(self):
        """
        Determine tag constants associated with annotation IDs.

        :return: map of tag IDs and respective tag constants
        """
        tagsMap = {}
        for ID in self.vcf_info_headers:
            num = self.output_vcf_nums[ID]
            if num == -1:  # the num field decides whether to split or not
                tagsMap[ID] = [TagConstants.INFO, TagConstants.SPLIT]
            else:
                tagsMap[ID] = [TagConstants.INFO, TagConstants.NOT_SPLIT]
        return tagsMap

    def _determine_info_annotation_value(self, vcf_record, ID, alt_index):
        """
        Determine the appropriate value that corresponds to a given ID.
        This method checks whether ID is present in the input Vcf record or not. If it is not found, then an appropriate
        missing value is computed.

        :param alt_index:
        :param vcf_record: input Vcf record
        :param ID: ID
        :return: value that corresponds to a given ID
        """
        if alt_index is None:
            val = self._determine_missing_value(ID)
            return val

        if ID in vcf_record.INFO:
            vals = vcf_record.INFO[ID]
            num = self.output_vcf_nums[ID]

            if isinstance(vals, list):
                if num == -1 and alt_index >= 0 and alt_index is not None:  # split by alternative allele
                    val = [vals[alt_index]]
                else:
                    val = vals
            else:  # force it to be a list
                val = [vals]
        else:
            val = self._determine_missing_value(ID)

        return val

    def _determine_missing_value(self, ID):
        """
        Determine the missing value that corresponds to a given ID.
        This method takes into account the number field and computes an appropriate missing value for a given ID.

        :param ID: ID
        :return: value that corresponds to a given ID
        """
        nsamples = len(self.vcf_reader.samples)
        num = self.output_vcf_nums[ID]
        if num == -2:
            return [""]*nsamples
        elif num == -1:
            return [""]
        elif num == 0:
            return [False]
        elif num is None:
            return [""]

        return [""]*num

    def _determine_matching_alt_indices(self, mut, record, build):
        """

        :param mut:
        :param record:
        :return:
        """
        indices = []
        if record.is_monomorphic:
            chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
            startPos = record.POS
            endPos = record.POS
            ref_allele = record.REF

            if self.match_mode == "exact":
                if mut.chr == chrom and mut.ref_allele == ref_allele:
                    indices = [-1]
            else:
                if mut.chr == chrom and int(mut.start) <= startPos and int(mut.end) >= endPos:
                    indices = [-1]
        else:
            # Iterate over all alternates in the record
            for index in xrange(0, len(record.ALT)):
                chrom = MutUtils.convertChromosomeStringToMutationDataFormat(record.CHROM)
                startPos = record.POS
                endPos = record.POS
                ref = str(record.REF)
                alt = str(record.ALT[index])
                ds_mut = MutUtils.initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build)

                if self.match_mode == "exact":
                    if mut.chr == ds_mut.chr and mut.ref_allele == ds_mut.ref_allele \
                        and mut.alt_allele == ds_mut.alt_allele and int(mut.start) == int(ds_mut.start) \
                        and int(mut.end) == int(ds_mut.end):
                        indices += [index]
                else:  # cases whether the match mode isn't exact
                    if mut.chr == ds_mut.chr and int(mut.start) == int(ds_mut.start) and int(mut.end) == int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) >= int(ds_mut.start) \
                        and int(mut.end) >= int(ds_mut.end) and int(mut.start) <= int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) <= int(ds_mut.start) and int(mut.end) >= int(ds_mut.end):
                        indices += [index]
                    elif mut.chr == ds_mut.chr and int(mut.start) <= int(ds_mut.start) \
                        and int(mut.end) <= int(ds_mut.end) and int(mut.end) >= int(ds_mut.start):
                        indices += [index]

        # if len(indices) == 0:
        #     indices = [None]

        return indices

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from Tabix indexed Vcf file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """

        # Get all the records corresponding to the given mutation
        mut_start = int(mutation.start)
        mut_end = int(mutation.end)
        vals = dict()
        record = None

        try:
            vcf_records = self.vcf_reader.fetch(mutation.chr, mut_start - 1, mut_end)  # query database for records
        except ValueError as ve:
            logging.getLogger(__name__).debug("Exception when looking for vcf records. Empty set of records being "
                                              "returned: " + repr(ve))
        else:
            # Process values.
            build = "hg19"
            for record in vcf_records:
                alt_indices = self._determine_matching_alt_indices(mutation, record, build)
                for ID in self.vcf_info_headers:
                    for alt_index in alt_indices:
                        val = self._determine_info_annotation_value(record, ID, alt_index)
                        if ID not in vals:  # dictionary of list of lists
                            vals[ID] = [val]
                        else:
                            vals[ID] += [val]

        if len(vals) == 0:
            for ID in self.vcf_info_headers:
                vals[ID] = [[""]]

        if record is None:
            msg = "Exception when looking for tsv records for chr%s:%s-%s. " \
                  "Empty set of records being returned." % (mutation.chr, mutation.start, mutation.end)
            logging.getLogger(__name__).debug(msg)

        tags = self._determine_tags()
        for ID in self.vcf_info_headers:
            if ID not in vals:  # this happens in cases where no matching records are found
                val = self._determine_missing_value(ID)
                vals[ID] = [val]

            if self.match_mode == "exact":
                # TODO: monomorphic records?
                val = []
                for v in vals[ID]:  # list of lists
                    v = string.join([str(v) if v is not None else "" for v in v], ",")
                    if v:
                        val += [v]
                val = string.join(val, "|")
            elif self.match_mode == "avg":
                if self.output_vcf_types[ID] in ("Integer", "Float"):  # integers and floats
                    num = self.output_vcf_nums[ID]
                    val = ""
                    if num in (None, -1, 0, 1,):
                        v = reduce(operator.add, vals[ID])  # flattens list of lists to a list
                        v = [v for v in v if v is not None and v != ""]  # discard non-numerics
                        if len(v) > 1:
                            val = str(float(sum(v))/len(v))
                        elif len(v) == 1:
                            val = str(float(v[0]))
                    else:
                        if len(vals[ID][0]) != num:  # addresses cases where no records were found
                            val = [""]
                        else:
                            nvals = len(vals[ID])  # number of lists
                            val = [[]]*num
                            for i in xrange(num):  # iterate over all numbers
                                for j in xrange(nvals):  # iterate over each list
                                    v = vals[ID][j][i]
                                    if v not in (None, "", "."):
                                        val[i] = val[i] + [v]
                                if len(val[i]) >= 1:
                                    val[i] = str(float(sum(val[i]))/len(val[i]))
                                else:
                                    val[i] = ""
                        val = string.join(val, ",")
                    self.output_vcf_types[ID] = "Float"
                elif self.output_vcf_types[ID] == "Character":
                    val = []
                    for v in vals[ID]:
                        val += [map(str, [string.join([v if v is not None else "" for v in v], ",")])]
                    val = string.join(val, "|")
                    self.output_vcf_types[ID] = "String"
                else:
                    val = []
                    for v in vals[ID]:
                        v = string.join(map(str, [v if v is not None else "" for v in v]), ",")
                        if v:
                            val += [v]
                    val = string.join(val, "|")
                    self.output_vcf_types[ID] = "String"
                    if self.output_vcf_nums[ID] == 0:
                        self.output_vcf_nums[ID] = None
            else:
                self.output_vcf_types[ID] = "String"  # everything is a String
                if self.output_vcf_nums[ID] == 0:
                    self.output_vcf_nums[ID] = None
                # values are delimited by pipe, and thus, the type is forced to be a String
                val = []
                for v in vals[ID]:  # list of lists
                    v = string.join([str(v) if v is not None else "" for v in v], ",")
                    if v:
                        val += [v]
                val = string.join(val, "|")

            mutation.createAnnotation(self.output_vcf_headers[ID], val, self.title, self.output_vcf_types[ID],
                                      self.output_vcf_descs[ID], tags=copy.copy(tags[ID]),
                                      number=self.output_vcf_nums[ID])

        return mutation