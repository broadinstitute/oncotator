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

import logging
import pysam
from oncotator.datasources.Datasource import Datasource
from oncotator.utils.TagConstants import TagConstants
import string
from oncotator.utils.MutUtils import MutUtils
from oncotator.MutationData import MutationData


class IndexedTsvDatasource(Datasource):
    """
    A datasource derived from a TSV file.  Expects a bgzipped vcf using Tabix.

    Instructions on how to index file using Tabix prior Oncotator:

    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.

    Multiple annotation data for the same mutation will be delimited by "|".
    """
    def __init__(self, src_file, title, version, colNames, indexColNames, annotationColNames, match_mode, colDataTypes):
        super(IndexedTsvDatasource, self).__init__(src_file, title=title, version=version)
        self.tsv_reader = pysam.Tabixfile(filename=src_file)  # initialize the tsv reader
        # Index the column names in tsv header
        self.tsv_headers = dict([(colName, index) for (index, colName) in enumerate(colNames)])
        # Annotation column names are the only columns that are used for annotating a mutation object.
        self.output_tsv_headers = dict([(colName, "_".join([self.title, colName])) for colName in annotationColNames])
        self.tsv_index = dict()
        self.tsv_index["chrom"] = self.tsv_headers[indexColNames[0]]
        self.tsv_index["start"] = self.tsv_headers[indexColNames[1]]
        self.tsv_index["end"] = self.tsv_headers[indexColNames[2]]
        if len(indexColNames) == 5:
            self.tsv_index["ref"] = self.tsv_headers[indexColNames[3]]
            self.tsv_index["alt"] = self.tsv_headers[indexColNames[4]]

        self.dataTypes = colDataTypes
        self.match_mode = match_mode

    def _is_matching(self, mut, tsv_record):

        chrom = tsv_record[self.tsv_index["chrom"]]
        startPos = tsv_record[self.tsv_index["start"]]
        endPos = tsv_record[self.tsv_index["end"]]
        build = "hg19"

        if self.match_mode == "exact":
            if "ref" in self.tsv_index and "alt" in self.tsv_index:  # ref and alt information is present
                ref = tsv_record[self.tsv_index["ref"]]
                alt = tsv_record[self.tsv_index["alt"]]
                if ref == "-" or alt == "-":  # addresses Mutation Annotation Format based tsv records
                    ds_mut = MutationData(chrom, startPos, endPos, ref, alt, build)
                else:  # addresses tsv records where the input isn't a Mutation Annotation Format file
                    ds_mut = MutUtils.initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build)

                if mut.chr == ds_mut.chr and mut.ref_allele == ds_mut.ref_allele \
                    and mut.alt_allele == ds_mut.alt_allele and int(mut.start) == int(ds_mut.start) \
                    and int(mut.end) == int(ds_mut.end):
                    return True
            else:  # do not use ref and alt information
                if mut.chr == chrom and int(mut.start) == int(startPos) and int(mut.end) == int(endPos):
                    return True
        else:
            if mut.chr == chrom and int(mut.start) == int(startPos) and int(mut.end) == int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) >= int(startPos) and int(mut.end) >= int(endPos) \
                and int(mut.start) <= int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) <= int(startPos) and int(mut.end) >= int(endPos):
                return True
            elif mut.chr == chrom and int(mut.start) <= int(startPos) and int(mut.end) <= int(endPos) \
                and int(mut.end) >= int(startPos):
                return True

        return False

    def _create_column_dict_from_tabix_index(self, mutation):
        mut_start = int(mutation.start)
        mut_end = int(mutation.end)
        chrom = mutation.chr
        vals = {}
        try:
            # tabix needs position - 1
            tsv_records = self.tsv_reader.fetch(chrom, mut_start - 1, mut_end, parser=pysam.asTuple())
            i = -1
            for i, tsv_record in enumerate(tsv_records):
                if not tsv_record:  # skip in case no records are found
                    continue

                logging.getLogger(__name__).debug("Got a record.")
                # Determine whether the new tsv record matches mutation or not
                if self._is_matching(mutation, tsv_record):
                    for colName in self.output_tsv_headers:
                        if colName.strip() == "":
                            continue
                        val = tsv_record[self.tsv_headers[colName]]
                        if colName not in vals:
                            vals[colName] = [val]
                        else:
                            vals[colName] += [val]
            logging.getLogger(__name__).debug("Processed %d records." % (i + 1))
        except ValueError as ve:
            msg = "Exception when looking for tsv records. Empty set of records being returned: " + repr(ve)
            logging.getLogger(__name__).debug(msg)
        return vals

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from a Tabix indexed TSV file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """

        vals = self._create_column_dict_from_tabix_index(mutation)

        for colName in self.output_tsv_headers:
            if colName in vals:  # this case happens when there are no matching records to be found
                val = string.join(vals[colName], "|")
            else:
                val = ""

            ds_type = self.dataTypes.get(colName, 'String')
            if self.match_mode == "exact":
                if "ref" not in self.tsv_index or "alt" not in self.tsv_index:
                    ds_type = "String"
            elif self.match_mode == "avg":
                if ds_type == "Integer" or ds_type == "Float":
                    ds_type = "Float"
                    if len(vals) != 0:
                        try:
                            v = [float(val) for val in vals[colName] if val not in ("", ".", "-",)]
                            if len(v) != 0:
                                val = str(sum(v)/len(v))
                            else:
                                val = ""
                        except ValueError:  # in cases where it's unsuccessful, default to overlap behavior
                            msg = "Exception when trying to cast %s value to float." % colName
                            logging.getLogger(__name__).warn(msg)
                            val = string.join(map(str, vals[colName]), "|")
                            ds_type = "String"
                else:
                    ds_type = "String"
            else:  # default: overlap
                ds_type = "String"

            if len(vals) == 0:
                msg = "Exception when looking for tsv records for chr%s:%s-%s. " \
                      "Empty set of records being returned." % (mutation.chr, mutation.start, mutation.end)
                logging.getLogger(__name__).debug(msg)

            mutation.createAnnotation(self.output_tsv_headers[colName], val, self.title, annotationDataType=ds_type,
                                      tags=[TagConstants.INFO, TagConstants.NOT_SPLIT])

        return mutation