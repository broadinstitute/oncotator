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
import csv
import string
import vcf
import logging
import collections
import itertools
import copy
import tempfile
import re
import os

from oncotator.utils.TsvFileSorter import TsvFileSorter
from oncotator.output.VcfOutputAnnotation import VcfOutputAnnotation
from oncotator.utils.TagConstants import TagConstants
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.SampleNameSelector import SampleNameSelector


class OutputDataManager:
    def __init__(self, configTable, muts, comments=None, md=None, path=os.getcwd()):
        """
        Initialize an instance of OutputDataManager.

        This class is used by the VcfOutputRenderer, but could be used elsewhere.

        :param muts: generator of MutationData objects
        :param path: location where to write temporary files
        :param configTable: tabular representation of the config
        :param comments: lines as key=value pairs (used in the header)
        :param md: meta-data
        """
        comments = [] if comments is None else comments
        md = [] if md is None else md

        self.logger = logging.getLogger(__name__)
        self.delimiter = "\t"
        self.lineterminator = "\n"
        self.comments = comments
        self.configTable = configTable
        self.metadata = md
        self.sampleNames = []
        self.mutation, self.mutations = self._fetchFirstMutation(muts)
        self.annotationTable, self.reverseAnnotationTable = self._createTables(self.metadata, self.mutation)
        self.chroms, self.sampleNames, self.filename = self._writeMuts2Tsv(self.mutations, path)

    def getSampleNames(self):
        return self.sampleNames

    def _fetchFirstMutation(self, muts):
        """
        Get the first mutation from the generator of mutations.

        :param muts: generator of mutations
        :return: first mutation
        """
        mut = None
        for mutation in muts:
            mut = copy.deepcopy(mutation)
            lst = [muts, (mut for mut in [mut])]
            muts = itertools.chain(*lst)
            break

        return mut, muts

    def getSortedTsvFilename(self, path):
        """


        :param path:
        :return:
        """
        chrom2HashCode = MutUtils.createChrom2HashCodeTable(self.chroms)
        tsvFileSorter = TsvFileSorter(self.filename)
        sortedTempTsvFile = tempfile.NamedTemporaryFile(dir=path, delete=False)
        func = lambda val: (chrom2HashCode[val["chr"]], int(val["start"]), val["alt_allele"])
        self.logger.debug("Sorting tmp tsv %s->%s", self.filename, sortedTempTsvFile.name)
        tsvFileSorter.sortFile(sortedTempTsvFile.name, func)
        os.remove(self.filename)

        return sortedTempTsvFile.name

    def _writeMuts2Tsv(self, muts, path):
        """
        Given a mutation generator, this methods writes a tab separated file for all mutations in the mutation
        generator. In addition, this method computes the appropriate sample name in scenarios where the mutation is
        missing sample name annotation. It also computes a list of all chromosomes and sample names contained within
        the generator.

        :param path: temporary filename
        :param muts: generator object with mutations
        """

        sampleNames = set()
        chroms = set()

        writer = None

        # create a temporary file to write tab-separated file
        tempTsvFile = tempfile.NamedTemporaryFile(dir=path, delete=False)
        self.logger.debug("Creating intermediate tsv file at %s" % tempTsvFile.name)

        mutAttributeNames = []
        sampleNameSelector = SampleNameSelector(self.mutation,
                                                configFile=self.configTable.getConfigFilename(),
                                                section="OTHER")

        with open(tempTsvFile.name, 'w') as fptr:
            ctr = 0
            sampleNameAnnotationName = sampleNameSelector.getOutputAnnotationName()
            sampleNameSource = sampleNameSelector.getAnnotationSource()

            for mut in muts:
                if len(mutAttributeNames) == 0:
                    mutAttributeNames = mut.getAttributeNames()

                sampleName = sampleNameSelector.getSampleName(mut)
                if sampleName is not None:
                    if mut.get(sampleNameAnnotationName, None) is None:
                        mut.createAnnotation(sampleNameAnnotationName, sampleName, sampleNameSource)
                    sampleNames.add(sampleName)

                # Parse chromosome
                chroms.add(mut.chr)

                updated_start, updated_ref_allele, updated_alt_allele = MutUtils.retrieveMutCoordinatesForRendering(mut)
                mut.ref_allele = updated_ref_allele
                mut.alt_allele = updated_alt_allele
                mut.start = updated_start

                if ctr == 0:
                    fieldnames2Render = MutUtils.getAllAttributeNames(mut)
                    if sampleNameAnnotationName is not None:
                        fieldnames2Render += [sampleNameAnnotationName]
                    for fieldname in fieldnames2Render:  # fieldnames that start "_" aren't rendered
                        if fieldname.startswith("_"):
                            fieldnames2Render.remove(fieldname)

                    writer = csv.DictWriter(fptr, fieldnames2Render, extrasaction='ignore', delimiter=self.delimiter,
                                            lineterminator=self.lineterminator)
                    writer.writeheader()

                writer.writerow(mut)

                ctr += 1
                if (ctr % 1000) == 0:
                    self.logger.info("Wrote " + str(ctr) + " mutations to tsv.")

        sampleNames = list(sampleNames)
        sampleNames.sort()
        chroms = list(chroms)

        return chroms, sampleNames, tempTsvFile.name

    def getHeader(self):
        """
        Returns correctly formatted header (meta-information and header lines).

        :return: header (meta-information and header lines)
        """
        return self._createHeader(self.comments, self.delimiter, self.lineterminator)

    def getOutputAnnotation(self, name):
        """
        Return VcfOutputAnnotation for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: instance of VcfOutputAnnotation
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is None:
            return VcfOutputAnnotation(name, "FORMAT", False, "INPUT", "String", "Unknown", None)
        return annotation

    def getFieldType(self, name):
        """
        Return field type (for example, INFO) for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: field type
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getFieldType()
        return "FORMAT"

    def getFieldID(self, name):
        """
        Return field ID for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: field ID
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getID()
        return name

    def getFieldDataType(self, name):
        """
        Return data type for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: data type
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getDataType()
        return "String"

    def getFieldNum(self, name):
        """
        Return number for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: number
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getNumber()
        return None

    def getFieldDesc(self, name):
        """
        Return description for the given field name.

        :param name: field name (or, unmapped field ID)
        :return: description
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getDescription()
        return "Unknown"

    def getFieldDatasource(self, name):
        """
        Returns whether the annotation associated with the given field name was generated by the program or the user.

        :param name: field name (or, unmapped field ID)
        :return: "INPUT", etc.
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.getDatasource()
        return "INPUT"

    def isFieldSplit(self, name):
        """
        Returns whether the given field name's value was split by the alternate allele or not.

        :param name: field name (or, unmapped field ID)
        :return: true or false
        """
        annotation = self.annotationTable.get(name, None)
        if annotation is not None:
            return annotation.isSplit()
        return False

    def getAnnotationNames(self, fieldType):
        """
        Returns a list of all field names associated with the given field type.

        :param fieldType: type of Vcf field (for example, INFO)
        :return: list of all field names associated with the field type
        """
        if fieldType in self.reverseAnnotationTable:
            return self.reverseAnnotationTable[fieldType]
        return []

    def _createHeader(self, comments=None, delimiter="\t", lineterminator="\n"):
        """
        Constructs VCF file header (meta-information and header lines).
        First, this method adds all meta-information lines as key=value pairs.
        Second, this method adds information lines describing the INFO, FILTER and FORMAT entries.
        Lastly, this method constructs the header line depending upon whether "FORMAT" tags are available or not.

        :param comments: lines as key=value pairs
        :param delimiter: special character (for example, tab) that separators words
        :param lineterminator: special character (for example, newline) signifying the end of line
        :return: VCF file header
        """
        comments = [] if comments is None else comments

        headers = ["##fileformat=VCFv4.1"]  # 'fileformat' is a required field; fixed since output vcf will be v4.1
        contigs = []
        alts = []

        for comment in comments:
            if comment:
                if comment.startswith("##ALT"):
                    headers += [comment]
                elif comment.startswith('##contig'):
                    headers += [comment]
                else:
                    pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
                    match = pattern.match(comment)
                    if match is not None:
                        if match.group("key") == "fileformat":
                            continue
                        headers += [comment]
                    else:
                        key = "comment"
                        val = string.replace(comment, "#", "")
                        val = string.replace(val, " ", "_")
                        if val.startswith("Oncotator"):
                            key = "oncotator_version"
                        comment = string.join([key, val], "=")
                        comment = string.join(["##", comment], "")
                        headers += [comment]

        annotations = self.annotationTable.values()
        gt = False
        for annotation in annotations:
            fieldType = annotation.getFieldType()

            ID = annotation.getID()
            if ID == "GT":
                gt = True

            description = annotation.getDescription()
            dataType = annotation.getDataType()
            num = annotation.getNumber()

            headers += [self._annotation2str(fieldType=fieldType, ID=ID, desc=description, dataType=dataType, num=num)]

        if not gt:
            headers += [self._annotation2str(fieldType="FORMAT", ID="GT", desc="Genotype", dataType="String", num=1)]

        # Add alts and contigs information
        headers += alts
        headers += contigs

        # Add header line
        if len(self.sampleNames) != 0:  # 'FORMAT' tags exist
            headers += [string.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], delimiter)]
            headers[len(headers)-1] += string.join(["", 'FORMAT'] + self.sampleNames, delimiter)
        else:
            headers += [string.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], delimiter)]

        headers = collections.OrderedDict.fromkeys(headers).keys()  # remove any duplicates and keep the original order
        headers = [header for header in headers if header not in ("", None,)]
        header = string.join(headers, lineterminator)
        return header

    def _getFieldnames(self, md, mut):
        """
        Given metadata and annotations from the first mutation, this function determines which fieldnames to render.

        :param md: metadata
        :param mut: first mutation
        :return: fieldnames
        """
        fieldnames = []
        if mut is not None:
            if md is not None and len(md) != 0:
                fieldnames = md.keys()
            fieldnames = set(fieldnames)  # convert list to set
            fieldnames = fieldnames.union(mut.keys())
            fieldnames = fieldnames.difference(['chr', 'start', 'end', 'ref_allele', 'alt_allele', 'alt_allele_seen',
                                                'build'])
            fieldnames = list(fieldnames)  # convert back to list
        else:
            if md is not None and len(md) != 0:
                fieldnames = md.keys()

        return fieldnames

    def _correctFieldName(self, fieldName):
        """
        Replaces unwanted characters in the field name

        :param fieldName:
        :return: corrected field name
        """
        fieldName = MutUtils.replaceChrs(fieldName, "=; :", "~|_>")  # Replace whitespace and other characters
        if fieldName.endswith("__FORMAT__"):  # Drop "__FORMAT__" from the end
            fieldName = fieldName[0:len(fieldName)-len("__FORMAT__")]
        elif fieldName.endswith("__INFO__"):  # Drop "__INFO__" from the end
            fieldName = fieldName[0:len(fieldName)-len("__INFO__")]
        return fieldName

    def _createTables(self, md, mut):
        """
        Constructs two tables, one that maps each field name to an VcfOutputAnnotation object and the other that lists
        all field names corresponding to a field type.
        First, this method computes a set of field names that are unique to both, meta-data and the first mutation.
        Second, this method resolves properties such as field type, ID and data type for every field name.

        :param md: meta-data
        :param mut: first mutation
        :return: two tables, one that maps each field name to an VcfOutputAnnotation object and the other that lists all
        field names corresponding to a field type
        """
        fieldnames = self._getFieldnames(md, mut)

        table = dict()
        revTable = dict()
        for fieldName in fieldnames:
            if fieldName.startswith("_"):  # skip rendering any fieldnames that start with a "_"
                continue

            if fieldName in md:
                annotation = md[fieldName]
            elif fieldName in mut:
                annotation = mut.getAnnotation(fieldName)

            num = annotation.getNumber()
            tags = annotation.getTags()
            dataType = annotation.getDataType()
            desc = annotation.getDescription()
            src = annotation.getDatasource()

            name = self._correctFieldName(fieldName)

            fieldType = self._resolveFieldType(name, tags)
            ID = self._resolveFieldID(fieldType, name)
            dataType = self._resolveFieldDataType(fieldType, ID, dataType)
            desc = self._resolveFieldDescription(fieldType, name, desc)

            isSplit = False
            if fieldType in ("INFO", "FORMAT",):  # determine whether to split or not for only INFO and FORMAT fields
                isSplit = self._determineIsSplit(fieldName, num, fieldType, tags)
            table[fieldName] = VcfOutputAnnotation(ID, fieldType, isSplit, src, dataType, desc, num)
            if fieldType not in revTable:
                revTable[fieldType] = [fieldName]
            else:
                revTable[fieldType] += [fieldName]

        return table, revTable

    def _determineIsSplit(self, name, num, field_type, tags=None):
        """
        Determines whether a given name's value was split by the alternate allele or not.
        This method implements the following decision tree for the number and field type corresponding to the ID:
        (1) Number is -2 (thus, G): by default it is assumed that the field name's value was NOT split by the alternate
        allele but is when the field type is FORMAT and the field name appears in the SPLIT_TAGS section of the config
        file.
        (2) Number is -1 (thus, A): by default it is assumed that the field name's value was split by the alternate allele
        but is NOT when the field name for a given field type appears in the NOT_SPLIT_TAGS section of the config file.
        (3) Number is None (thus, .): by default it is assumed that the field name's value was NOT split by the alternate
        allele but is when the field name for a given field type appears in the NOT_SPLIT_TAGS section of the config file.
        (4) Number is an integer: the default it is assumed that the field name's values was NOT split by the alternate
        allele but is when the field name for a given field type appears in the NOT_SPLIT_TAGS section of the config file.

        :param name: field name
        :param num: integer that describes the number of values that can be included with the field
        :param field_type: type of Vcf field (for example, INFO)
        :return: whether the field ID's value was split by the alternate or not
        """
        tags = [] if tags is None else tags

        if field_type != "FORMAT" and num == -2:
            logging.getLogger(__name__).warn("Non-FORMAT field type (%s) with Number=G name: %s  A datasource or VCF file may be misconfigured." % (field_type, name))

        # Whether the annotation indicates split, not split, or unspecified
        is_annotation_split = (TagConstants.SPLIT in tags)
        is_annotation_not_split = (TagConstants.NOT_SPLIT in tags)
        is_annotation_split_unknown = (not is_annotation_split) and (not is_annotation_not_split)

        # Whether the config file indicates split, not split, or unspecified
        is_config_split = self.configTable.isFieldNameInSplitSet(field_type, name)
        is_config_not_split = self.configTable.isFieldNameInNotSplitSet(field_type, name)
        is_config_unknown = (not is_config_not_split) and (not is_config_split)

        if num == -2:  # by the number of samples (Number=G)
            # Always return false, but show log message if config file or annotation indicates it should be true
            if is_config_split or is_annotation_split:
                logging.getLogger(__name__).warn("Annotation or config file specifying is split for field type (%s) with Number=G name: %s  A datasource or VCF file may be misconfigured." % (field_type, name))
            return False

        elif num == -1:  # by the number of alternates (Number=A)
            # Always return true, but show log message if config file or annotation indicates it should be false
            if is_config_split or is_annotation_split:
                logging.getLogger(__name__).warn("Annotation or config file specifying is not split for field type (%s) with Number=A name: %s  A datasource or VCF file may be misconfigured." % (field_type, name))
            return True

        elif num is None or num > -1:  # number is unknown or greater than -1 (Use config file, then annotation to decide... otherwise, True for num=="." and False for num >=0)
            if is_config_split:
                return True
            elif is_config_not_split:
                return False
            elif is_annotation_split:
                return True
            elif is_annotation_not_split:
                return False

            # both unknown
            else:
                if num is None:
                    return True
                else:
                    return False

        # Should never get here.  Throw warning and return True
        logging.getLogger(__name__).warn("%s %s num: %d  -- number is unrecognized.  Assuming is split." % (field_type, name, num))
        return True

    def _resolveFieldType(self, name, tags):
        """
        Determine the type of Vcf field corresponding to the given name.
        First, this method inspects the tags section for known keys and attempts to determine the type of field for the
        given name.
        Second, this method inspects the config file and attempts to determine the type of field for the given name.
        Third, this method inspects the reserved key words and attempts to determine the type of field for the given
        name.
        Fourth, this method defaults the field type to the INFO field.

        :param name: unmapped field ID
        :param tags: list of keys associated passed in the mut object
        :return: type of Vcf field (for example, INFO)
        """
        m = {TagConstants.INFO: "INFO", TagConstants.FORMAT: "FORMAT", TagConstants.FILTER: "FILTER",
             TagConstants.ID: "ID", TagConstants.QUAL: "QUAL"}

        for tag in tags:
            if tag in m.keys():
                return m[tag]

        if name in self.configTable.getInfoFieldNames():
            return "INFO"
        elif name in self.configTable.getFormatFieldNames():
            return "FORMAT"
        elif name in self.configTable.getFilterFieldNames():
            return "FILTER"
        elif name in self.configTable.getOtherFieldNames():
            return self.configTable.getOtherFieldID(name)

        if name.upper() in vcf.parser.RESERVED_INFO:
            return "INFO"
        elif name.upper() in vcf.parser.RESERVED_FORMAT:
            return "FORMAT"

        return "INFO"

    def _resolveFieldID(self, fieldType, name):
        """
        Determines the field ID corresponding to the name and field type.
        This method retrieves from the config file the field ID corresponding to both, the name and the field type.
        In cases where the name is not found, this method returns the name itself as the field ID.

        :param fieldType: type of Vcf field (for example, INFO)
        :param name: unmapped field ID
        :return: mapped field ID
        """
        if fieldType == "FORMAT":
            return self.configTable.getFormatFieldID(name)
        elif fieldType == "INFO":
            return self.configTable.getInfoFieldID(name)
        elif fieldType == "FILTER":
            return self.configTable.getFormatFieldID(name)

        return name

    def _resolveFieldDataType(self, fieldType, ID, dataType):
        """
        Determines whether the value corresponding to a reserved field ID is an integer, a float, a character or a
        string.

        :param fieldType: type of Vcf field (for example, INFO)
        :param ID: field ID
        :param dataType: type (integer/float/character/string) passed in the mut object
        :return: corrected type (integer/float/character/string)
        """
        if dataType == "String":
            if fieldType == "FILTER":
                return vcf.parser.RESERVED_INFO.get(ID, dataType)
            elif fieldType == "FORMAT":
                return vcf.parser.RESERVED_FORMAT.get(ID, dataType)

        return dataType

    def _resolveFieldDescription(self, fieldType, name, description):
        """
        Determines description corresponding to a given field ID.
        This method overwrites the default description passed in the mut object with corresponding description in
        the config file when present.

        :param fieldType: type of Vcf field (for example, INFO)
        :param name: field name
        :param description: description passed in mut object
        :return: description (default: 'Unknown'), must be surrounded by double-quotes
        """
        if description is None or description.lower() == "unknown":
            if fieldType == "FILTER" and name in self.configTable.getFilterFieldNames():
                return self.configTable.getFilterFieldNameDescription(name)
            elif fieldType == "FORMAT" and name in self.configTable.getFormatFieldNames():
                return self.configTable.getFormatFieldNameDescription(name)
            elif fieldType == "INFO" and name in self.configTable.getInfoFieldNames():
                return self.configTable.getInfoFieldNameDescription(name)

        if description is None or description == "":
            return "Unknown"

        return description

    def _annotation2str(self, fieldType, ID, desc="Unknown", dataType="String", num=None):
        """
        Converts FORMAT, INFO and FILTER to their respective meta-information descriptions.

        :param fieldType: type of Vcf field (for example, INFO)
        :param ID: field ID
        :param desc: description
        :param dataType: identifies whether the value corresponding to the ID is an integer, a float, a character or a
        string
        :param num: integer that describes the number of values that can be included with the field
        :return: correct meta-information description
        """
        num = "." if num is None else num

        if fieldType in ("FORMAT", "INFO",):
            return "##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % (fieldType, ID, num, dataType, desc)
        elif fieldType in ("FILTER",):
            return "##%s=<ID=%s,Description=\"%s\">" % (fieldType, ID, desc)

        return ""