# LICENSE_GOES_HERE
import logging
from oncotator.utils.MissingRequiredAnnotationException import MissingRequiredAnnotationException


"""
Created on Nov 13, 2012

@author: lichtens
"""
from oncotator.MutationData import MutationData
from oncotator.utils.MutationValidationFailureException import MutationValidationFailureException
import re
import collections
import string
import shutil


class MutUtils(object):
    """
    Static class containing utility functions for Mutations. 
    """
    proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)[_]*[A-Z]{0,1}([0-9]*)")
    SAMPLE_NAME_ANNOTATION_NAME = "sample_name"
    PRECEDING_BASES_ANNOTATION_NAME = "_preceding_bases"

    def __init__(self, params):
        """
        Constructor -- should never be called.
        """
        pass

    @staticmethod
    def isSNP(m):
        if len(m.ref_allele) > 1:
            return False
        if m.alt_allele not in ["A", "C", "T", "G"]:
            return False
        return True

    @staticmethod
    def isDeletion(m):
        if len(m.ref_allele) > len(m.alt_allele):
            return True
        return False

    @staticmethod
    def isInsertion(m):
        if len(m.ref_allele) < len(m.alt_allele):
            return True
        return False

    @staticmethod
    def determineVariantType(m):
        if MutUtils.isSNP(m):
            return "snp"
        elif MutUtils.isDeletion(m):
            return "del"
        elif MutUtils.isInsertion(m):
            return "ins"
        return "unknown"

    @staticmethod
    def initializeMutFromAttributes(chrom, startPos, endPos, ref, alt, build):
        mut = MutationData(chrom, startPos, endPos, ref, alt, build)
        varType = MutUtils.determineVariantType(mut)

        if varType == "snp":  # Snps
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue="")
        if varType == "del":  # deletion
            preceding_bases, updated_ref_allele, updated_start, updated_end =\
                MutUtils.retrievePrecedingBasesForDeletions(mut)
            mut.ref_allele = updated_ref_allele
            mut["ref_allele"] = updated_ref_allele
            mut.alt_allele = "-"
            mut["alt_allele"] = "-"
            mut.start = updated_start
            mut["start"] = updated_start
            mut.end = updated_end
            mut["end"] = updated_end
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME,
                                 annotationValue=preceding_bases)
        elif varType == "ins":  # insertion
            preceding_bases, updated_alt_allele, updated_start, updated_end = \
                MutUtils.retrievePrecedingBasesForInsertions(mut)
            mut.ref_allele = "-"
            mut["ref_allele"] = "-"
            mut.alt_allele = updated_alt_allele
            mut["alt_allele"] = updated_alt_allele
            mut.start = updated_start
            mut["start"] = updated_start
            mut.end = updated_end
            mut["end"] = updated_end
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME,
                                 annotationValue=preceding_bases)

        return mut

    @staticmethod
    def retrievePrecedingBasesForInsertions(m):
        ref_allele = m.ref_allele
        alt_allele = m.alt_allele
        start = int(m.start)

        preceding_bases = ref_allele
        updated_alt_allele = alt_allele[len(preceding_bases):]
        updated_start = start + len(preceding_bases) - 1
        updated_end = start + len(preceding_bases)

        return preceding_bases, updated_alt_allele, updated_start, updated_end

    @staticmethod
    def retrievePrecedingBasesForDeletions(m):
        ref_allele = m.ref_allele
        alt_allele = m.alt_allele
        start = int(m.start)

        preceding_bases = alt_allele
        updated_ref_allele = ref_allele[len(preceding_bases):]
        updated_start = start + len(preceding_bases)
        updated_end = updated_start + len(updated_ref_allele) - 1

        return preceding_bases, updated_ref_allele, updated_start, updated_end

    @staticmethod
    def removeDir(currentDir):
        shutil.rmtree(path=currentDir, ignore_errors=True)

    @staticmethod
    def createChrom2HashCodeTable(chroms):
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

    @staticmethod
    def replaceChrs(text, frm, to):
        tbl = string.maketrans(frm, to)
        return text.translate(tbl)

    @staticmethod
    def getAllAttributeNames(mut):
        """
        :param mut: mutation object
        :return: set of attribute names that are encapsulated in the mutation
        """
        attrs = []
        if mut is not None:
            attrs = mut.keys() + mut.getAttributeNames()
            return collections.OrderedDict.fromkeys(attrs).keys()
        return attrs

    @staticmethod
    def str2bool(v):
        """ Given an input string, v, returns whether the input string is a boolean 
            or not. """
        return v.lower() in ("yes", "true", "t", "1", "y")

    @staticmethod
    def enum(*sequential, **named):
        enums = dict(zip(sequential, range(len(sequential))), **named)
        reverse = dict((value, key) for key, value in enums.iteritems())
        enums['reverse_mapping'] = reverse
        return type('Enum', (), enums)

    @staticmethod
    def getTokens(text, delimiter="\t", lineterminator="\n"):
        lines = text.split(lineterminator)
        if not lines[len(lines)-1]:
            lines = lines[0:len(lines)-1]
        toks = []
        for line in lines:
            toks += line.split(delimiter)
        return toks

    @staticmethod
    def convertChromosomeStringToMutationDataFormat(chrom, build="hg19"):
        """ This method only covers strings.  This method does not check that the chromosome exists for a given genome build. 
        
        The standard for MutationData is "1", "X", "M" -- not: "chrX"
        
        Specifically:
        Converts MT to M
        Strips "chr"
        If the chromosome is 23 or 24 and the build begins with "hg", then convert to "X" and "Y", respectively.

        If a chromosome contains "<" or  ">", then replace those with "".

        TODO:  This will be changed to a framework that better supports build-specific classes that do this conversion.
        """
        result = chrom
        if chrom.startswith('chr'):
            result = result.replace('chr', '')
        if chrom == "MT":
            result = "M"
        if chrom == "'M_rCRS'":
            result = "M"

        if build.startswith("hg") and (chrom == "23" or chrom == "24"):
            if chrom == "23":
                result = "X"
            if chrom == "24":
                result = "Y"

        result = result.replace('<', '')
        result = result.replace('>', '')

        return result
        
    @staticmethod
    def getAnnotationsByDatasource(mutation, dsTitle):
        """ Given a datasource title, return all of the annotation names with
            that datasource. """
        result = []
        ks = mutation.keys()
        for k in ks:
            if mutation.getAnnotation(k).getDatasource() == dsTitle:
                result.append(k)
        return result
    
    @staticmethod
    def getUnknownAnnotations(mutation):   
        """ Returns annotations where source is unknown, None, or blank. """  
        result = MutUtils.getAnnotationsByDatasource(mutation, "Unknown")
        result.extend(MutUtils.getAnnotationsByDatasource(mutation, ""))
        result.extend(MutUtils.getAnnotationsByDatasource(mutation, None))
        return result
    
    @staticmethod
    def prettyPrint(mutation):
        for annotation in mutation: 
            print mutation[annotation]

    @staticmethod
    def validateMutation(mutation):
        """ Does some basic sanity checks that the given mutationData is coherent.
        
        Checks:
            All attributes are present and not None
            No attribute is blank.
            Make sure that the dict interface matches the attribute.
            
            TODO: Validate chromosome value given the genome build
            
        Throws MutationValidationFailureException.  Otherwise, returns True
        """

        attribute = "chr"
        noneAttributes = []
        for attribute in MutationData.attributes:
            if mutation[attribute] is None:
                noneAttributes.append(attribute)
        if len(noneAttributes) > 0:
            raise MutationValidationFailureException("None values found for attributes: " + str(noneAttributes))

        blankAttributes = []
        for attribute in MutationData.attributes:
            if mutation[attribute] is '':
                blankAttributes.append(attribute)
        if len(blankAttributes) > 0:
            raise MutationValidationFailureException("Blank values found for attributes: " + str(blankAttributes))

        noMatchAttributes = []
        for attribute in MutationData.attributes:
            if mutation[attribute] != mutation.__dict__[attribute]:
                noMatchAttributes.append(attribute)
        if len(noMatchAttributes) > 0:
            raise MutationValidationFailureException("Attribute did not match dictionary value for " +
                                                     str(noMatchAttributes) + " (dict, attribute): " +
                                                     mutation[attribute] + ", " + mutation.__dict__[attribute])
        
        if mutation.chr.startswith("chr"):
            raise MutationValidationFailureException("Chromosome value started with chr: " + str(mutation.chr))
        
        if mutation.chr == "MT":
            raise MutationValidationFailureException("Mitochondria must be M, not MT.")
        # TODO: Check for valid chromosome values given the genome build
        
        return True

    @staticmethod
    def retrieveMissingAnnotations(m, annotationNames):
        """
        Get a subset of the given annotation names that are not present in the mutation.

        The returned list is sorted.

        :param m: MutationData
        :param annotationNames: List of string annotations that should be searched for in mut.
        :return missingAnnotations: List of annotations in annotationNames that were not present in the mutation.
        """
        annotationNamesInMut = set(m.keys())
        result = list(set(annotationNames).difference(annotationNamesInMut))
        return sorted(result)

    @staticmethod
    def extractProteinPosition(proteinChange):
        """
        Extract start and end in a protein change string, such as:
        """
        result = ['','']

        m = MutUtils.proteinRegexp.search(proteinChange)
        if m is not None:
            result[0] = m.group(1)
            if (m.group(2) is None) or (m.group(2).strip() == ""):
                result[1] = result[0]
            else:
                result[1] = m.group(2)
        return result

    @staticmethod
    def create_variant_key_by_mutation(m, other_info=""):
        return MutUtils.create_variant_key(m.chr, m.start, m.end, m.ref_allele, m.alt_allele, other_info)

    @staticmethod
    def create_variant_key(chrom, start, end, ref_allele, alt_allele, other_info=""):
        return "%s_%s_%s_%s_%s_%s" % (chrom, start, end, ref_allele, alt_allele, other_info)

    @staticmethod
    def createFieldsMapping(headers, annotations, alternativeDictionary, isRenderInternalFields=True,
                            exposedFields=None, prepend="i_"):
        """ Creates a dictionary of the output maf file headers to the annotations.
        :param headers:  optional and required fields
        :param annotations: annotations seen on a mutation
        :param alternativeDictionary: mapping from headers to acceptable annotation names.
        :param isRenderInternalFields: should the resulting dictionary include
                unused annotations (as [prepend]annotation_name)?
        :param exposedFields: These are fields that, if seen, should not be treated as internal.
        :param prepend: prepend to use for internal annotations (i.e. ones that are not seen in the required or optional fields)


        Output:
        Final headers that should be used
         """

        if prepend is None:
            prepend = ""

        if exposedFields is None:
            exposedFields = set()

        result = dict()
        for h in headers:
            hdr = h.lower()
            result[h] = h
            if h not in annotations:
                if hdr in alternativeDictionary.keys():
                    alternatives = alternativeDictionary[hdr]
                    for alternative in alternatives:
                        if alternative in annotations:
                            result[h] = alternative
                            break

        if isRenderInternalFields:
            sAnnotations = set(annotations)
            internalFields = sAnnotations.difference(result.values())
            for i in internalFields:
                if not i.startswith('_') and i is not "transcripts":
                    if prepend.strip() == "" or i.startswith(prepend) or i in exposedFields:
                        result[i] = i
                    else:
                        result[prepend + i] = i

        return result

    @staticmethod
    def retrieveMutCoordinatesForRendering(mut):
        updated_start = mut.start
        updated_ref_allele = mut.ref_allele
        updated_alt_allele = mut.alt_allele
        if mut.ref_allele == "-":  # detects insertions in cases where the input is a maf
            if MutUtils.PRECEDING_BASES_ANNOTATION_NAME in mut:
                updated_ref_allele, updated_alt_allele, updated_start = \
                    MutUtils.retrievePrecedingBaseFromAnnotationForInsertions(mut)
            else:
                updated_ref_allele, updated_alt_allele, updated_start = \
                    MutUtils.retrievePrecedingBaseFromReference(mut)
        elif mut.alt_allele == "-":  # detects deletions in cases where the input is a maf
            if MutUtils.PRECEDING_BASES_ANNOTATION_NAME in mut:
                updated_ref_allele, updated_alt_allele, updated_start = \
                    MutUtils.retrievePrecedingBaseFromAnnotationForDeletions(mut)
            else:
                updated_ref_allele, updated_alt_allele, updated_start = \
                    MutUtils.retrievePrecedingBaseFromReference(mut)
        elif mut.ref_allele == mut.alt_allele:  # detects monomorphic SNPs
            updated_alt_allele = ""

        return updated_start, updated_ref_allele, updated_alt_allele

    @staticmethod
    def retrievePrecedingBaseFromAnnotationForDeletions(mut):
        preceding_bases = mut[MutUtils.PRECEDING_BASES_ANNOTATION_NAME]
        alt_allele = preceding_bases
        ref_allele = preceding_bases + mut.ref_allele
        updated_start = mut.start - len(preceding_bases)
        return ref_allele, alt_allele, updated_start

    @staticmethod
    def retrievePrecedingBaseFromAnnotationForInsertions(mut):
        preceding_bases = mut[MutUtils.PRECEDING_BASES_ANNOTATION_NAME]
        alt_allele = preceding_bases + mut.alt_allele
        ref_allele = preceding_bases
        updated_start = mut.start - len(preceding_bases) + 1

        return ref_allele, alt_allele, updated_start

    @staticmethod
    def retrievePrecedingBaseFromReference(m):
        updated_start = m.start
        ref_allele = m.ref_allele
        if ref_allele == "-":
            ref_allele = "."

        alt_allele = m.alt_allele
        if alt_allele == "-":
            alt_allele = "."

        if "ref_context" in m:
            ref_context = m['ref_context']
        else:
            raise MissingRequiredAnnotationException("Missing ref_context annotation in mutation. "
                                                     "Please add ref_[genome_build] to your dbDir.")

        # Insert only
        if ref_allele == ".":
            if ref_context == "" or (len(ref_context) < 11):
                logging.getLogger(__name__).warn("WARNING: Could not find ref_context when required.  "
                                                 "Unable to render: " + m.chr + ":" + m.start + "-" + m.end)
                return None
            ref_allele = ref_context[10].upper()
            alt_allele = ref_allele + alt_allele

        # Deletion only
        if alt_allele == ".":
            if ref_context == "" or (len(ref_context) < 11):
                logging.getLogger(__name__).warn("WARNING: Could not find ref_context when required.  "
                                                 "Unable to render: " + m.chr + ":" + m.start + "-" + m.end)
                return None
            ref_allele = ref_context[9].upper() + ref_allele
            alt_allele = ref_context[9].upper()
            updated_start = str(int(m.start) - 1)

        return ref_allele, alt_allele, updated_start
