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
        tsvFileSorter.sortFile(sortedTempTsvFile.name, func)
        os.remove(self.filename)

        return sortedTempTsvFile.name

    def _writeMuts2Tsv(self, muts, path):
        """
        Given a mutation generator, this methods writes a tab separated file for all mutations in the mutation
        generator. In addition, this method computes the appropriate sample name in scenarios where the mutation is
        missing sample name annotation. It also computes a list of all chromosomes and sample names contained within
        the generator.

        :param filename: temporary filename
        :param muts: generator object with mutations
        """

        sampleNames = set()
        chroms = set()

        writer = None

        # create a temporary file to write tab-separated file
        tempTsvFile = tempfile.NamedTemporaryFile(dir=path, delete=False)
        self.logger.info("Creating intermediate tsv file...")

        sampleNameAnnotationNames = self.getAnnotationNames("SAMPLE_NAME")
        tumorSampleNameAnnotationNames = self.getAnnotationNames("SAMPLE_TUMOR_NAME")
        normalSampleNameAnnotationNames = self.getAnnotationNames("SAMPLE_NORMAL_NAME")

        mutAttributeNames = []

        with open(tempTsvFile.name, 'w') as fptr:
            ctr = 0
            for mut in muts:

                sampleName = None
                sampleNameAnnotationName = None

                if len(mutAttributeNames) == 0:
                    mutAttributeNames = mut.getAttributeNames()

                # Sample name annotation is present
                if len(sampleNameAnnotationNames) != 0:
                    sampleNameAnnotationName = sampleNameAnnotationNames[0]
                    sampleName = mut[sampleNameAnnotationName]
                # Both, tumor and normal sample name annotations are present
                elif len(tumorSampleNameAnnotationNames) != 0 and len(normalSampleNameAnnotationNames) != 0:
                    tumorSampleNameAnnotationName = tumorSampleNameAnnotationNames[0]
                    normalSampleNameAnnotationName = normalSampleNameAnnotationNames[0]
                    sampleName = string.join([mut[normalSampleNameAnnotationName],
                                              mut[tumorSampleNameAnnotationName]], sep="-")
                    sampleNameAnnotationName = MutUtils.SAMPLE_NAME_ANNOTATION_NAME
                    mut.createAnnotation(sampleNameAnnotationName, sampleName, "INPUT")
                    if ctr == 0:
                        self.logger.info("Sample name is the concatenation of %s and %s columns."
                                         % (normalSampleNameAnnotationName, tumorSampleNameAnnotationName))
                # Only tumor sample name is present
                elif len(tumorSampleNameAnnotationNames) != 0:
                    tumorSampleNameAnnotationName = tumorSampleNameAnnotationNames[0]
                    sampleName = mut[tumorSampleNameAnnotationName]
                    sampleNameAnnotationName = MutUtils.SAMPLE_NAME_ANNOTATION_NAME
                    mut.createAnnotation(sampleNameAnnotationName, sampleName, "INPUT")
                    if ctr == 0:
                        self.logger.info("Sample name is %s column." % tumorSampleNameAnnotationName)
                # Only normal sample name is present
                elif len(normalSampleNameAnnotationNames) != 0:
                    normalSampleNameAnnotationName = normalSampleNameAnnotationNames[0]
                    sampleName = mut[normalSampleNameAnnotationName]
                    sampleNameAnnotationName = MutUtils.SAMPLE_NAME_ANNOTATION_NAME
                    mut.createAnnotation(sampleNameAnnotationName, sampleName, "INPUT")
                    if ctr == 0:
                        self.logger.info("Sample name is %s column." % normalSampleNameAnnotationName)

                if sampleName is not None:
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

    def _determineIsSplit(self, name, num, fieldType, tags=None):
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
        :param fieldType: type of Vcf field (for example, INFO)
        :return: whether the field ID's value was split by the alternate or not
        """
        tags = [] if tags is None else tags

        if num == -2:  # by the number of samples
            if fieldType == "FORMAT":
                if TagConstants.SPLIT in tags:  # override the default using the tags section
                    return True
                elif self.configTable.isFieldNameInSplitSet(fieldType, name):  # override the default using the config file
                    return True
                else:
                    return False
        elif num == -1:  # by the number of alternates
            if TagConstants.NOT_SPLIT in tags:  # override the default using the tags section
                return False
            elif self.configTable.isFieldNameInNotSplitSet(fieldType, name):
                return False
            else:
                return True
        elif num is None:  # number is unknown
            if TagConstants.SPLIT in tags:  # override the default using the tags section
                return True
            elif self.configTable.isFieldNameInSplitSet(fieldType, name):  # override the default using the config file
                return True
            else:
                return False
        else:
            if TagConstants.NOT_SPLIT in tags:  # override the default using the tags section
                return True
            elif self.configTable.isFieldNameInSplitSet(fieldType, name):  # override the default using the config file
                return True
            else:
                return False

        return False

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