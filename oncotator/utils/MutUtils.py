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


"""
Created on Nov 13, 2012

@author: lichtens
"""
from oncotator.MutationData import MutationData
from oncotator.utils.MutationValidationFailureException import MutationValidationFailureException
import re
class MutUtils(object):
    """
    Static class containing utility functions for Mutations. 
    """

    proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)[_]*[A-Z]{0,1}([0-9]*)")

    def __init__(self,params):
        """
        Constructor -- should never be called.
        """
        pass
    
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

        if build.startswith("hg") and (chrom == "23" or chrom =="24"):
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
        return prettyStr
    
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
        
        for attribute in MutationData.attributes:
            noneAttributes = []
            if mutation[attribute] is None:
                noneAttributes.append(attribute)
        if len(noneAttributes) > 0:
            raise MutationValidationFailureException("None values found for attributes: " + str(noneAttributes))

        for attribute in MutationData.attributes:
            blankAttributes= []
            if mutation[attribute] is '':
                blankAttributes.append(attribute)
        if len(blankAttributes) > 0:
            raise MutationValidationFailureException("Blank values found for attributes: " + str(blankAttributes))
        
        for attribute in MutationData.attributes:
            noMatchAttributes= []
            if mutation[attribute] <> mutation.__dict__[attribute]:
                noMatchAttributes.append(attribute)
        if len(noMatchAttributes) > 0:
            raise MutationValidationFailureException("Attribute did not match dictionary value for " + str(noMatchAttributes) + " (dict, attribute): " + mutation[attribute] + ", " + mutation.__dict__[attribute])
        
        if mutation.chr.startswith("chr"):
            raise MutationValidationFailureException("Chromosome value started with chr: " + str(mutation.chr))
        
        if (mutation.chr == "MT"):
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
    def createFieldsMapping(headers, annotations, alternativeDictionary, isRenderInternalFields=True, exposedFields=set()):
        """ Creates a dictionary of the output maf file headers to the annotations.
        Input:
            headers -- optional and required fields
            annotations -- annotations seen on a mutation
            alternativeDictionary -- mapping from headers to acceptable annotation names.
            isRenderInternalFields -- should the resulting dictionary include
                unused annotations (as i_annotation_name)?
            exposedFields -- These are fields that, if seen, should not be treated as internal.

        Output:
        Final headers that should be used
         """
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

        if isRenderInternalFields:
            sAnnotations = set(annotations)
            internalFields = sAnnotations.difference(result.values())
            for i in internalFields:
                if not i.startswith('_') and i is not "transcripts":
                    if i.startswith("i_") or i in exposedFields:
                        result[i] = i
                    else:
                        result["i_" + i] = i

        return result