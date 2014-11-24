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


"""
Created on Oct 22, 2012

@author: lichtens
"""
from DuplicateAnnotationException import DuplicateAnnotationException
import collections
import logging
from Annotation import Annotation
from collections import OrderedDict
#from multiprocessing import Lock


class MutationData(collections.MutableMapping):
    """
    Intermediate class for storing a mutation.
    
    Usage notes:
    
        MutationData classes have set attributes:
            chr -- chromosome without "chr".  For example, "1" or "X" or "GL000194.1" or "MT" (without quotes)
            start -- start position (Assumes chromosome starts with position 1)
            end -- end position (Assumes chromosome starts with position 1)
            ref_allele -- reference at the specified positions
            alt_allele -- observed allele
            build -- genome build (e.g. mm9, hg18, hg19, or canFam2)
                        
            All are stored as strings and can be treated as annotations as well (see below)
            
        All other attributes are considered annotations and can only be accessed via the dictionary interface.
            For example, given MutationData m:
                print m['some_annotation']
                print m.chr
                
        Notes:
            len(m) will return the number of annotations not including the attributes above.
            Annotations can be accessed as if in a dictionary.  For example, tmp = m['my_annotation1']
            Annotations can be set as if the mutation was a dictionary, but this is not recommended.  Doing this will cause the annotation to only have a value -- no source, type, etc.
            
        MutationData is pickle-able.

    """
    
    """ internal annotations that will show as both annotations and attributes.   If this changes, updates should probably be made to the maflite config."""
    attributes = {"chr", "start", "end", "ref_allele", "alt_allele", "build"}

    def __init__(self, chr="", start="", end="", ref_allele="", alt_allele="", build=""):
        """
        Constructor
        """
        self.__dict__.update(locals())
        self.annotations = dict()

        for k in MutationData.attributes:
            self.annotations[k] = locals()[k]

#        self.lock = Lock()
        
    def createAnnotation(self, annotationName, annotationValue, annotationSource="Unknown", annotationDataType="String", annotationDescription="", newRequired=True, tags=None, number=None):
        """
        newRequired implies that this cannot update an existing value.  If a value exists, throw an exception.
        
        This method must be called to add an annotation to a mutation.  Do not use: mut['new_annotation_name'] = 'annotation_value'
        
        """
        tags = [] if tags is None else tags

#        self.lock.acquire()
        if newRequired and (annotationName in self.annotations.keys()) and (annotationName not in MutationData.attributes):
#            self.lock.release()
            if annotationValue == self.annotations[annotationName].value:
                logging.getLogger(__name__).warn("Attempting to create an annotation multiple times, but with the same value: " + str(annotationName) +  "  :  " + str(annotationValue))
            else:
                raise DuplicateAnnotationException('Attempting to create an annotation multiple times (' + annotationName + ') with old, new values of (' + str(self.annotations[annotationName].value) + ", " + str(annotationValue) + ")")
        if annotationName in MutationData.attributes:
            # FYI ... logging.getLogger(__name__).debug("Attempting to create an attribute with createAnnotation.  Should be using instance attribute setting.  x." + str(annotationName) + " = " + str(annotationValue) + " ... Ignoring annotationSource, but setting attribute.")
            self[annotationName] = annotationValue
        else:
            self.annotations[annotationName] = Annotation(annotationValue, annotationSource, annotationDataType, annotationDescription, tags=tags, number=number)
#        self.lock.release()
        
    def getAnnotation(self, annotationName):
        """ Returns the Annotation instance, rather than just the value """
        if annotationName in MutationData.attributes:
            return Annotation(self.__dict__[annotationName], "__ATTR__")
        return self.annotations[annotationName]

    def getAnnotations(self):
        """Returns a list of Annotation instances."""
        return self.annotations.values()

    def addAnnotations(self, annot_dict):
        """
        :param annot_dict: name:Annotation dictionary
        :return:
        """
        self.annotations.update(annot_dict)

    def addTagToAnnotation(self, annotationName, tag):
        """ Attach tag to a given annotation """
        if not (annotationName in MutationData.attributes):
            self.annotations[annotationName].addTag(tag)

    def getAttributeNames(self):
        return list(self.attributes)

    def print_annotations(self):
        """ Pretty print datasources, annotations, and annotations value."""
        outputs = list()
        for annotation in self.annotations:
            if annotation in ['build', 'chr', 'start', 'end', 'ref_allele', 'alt_allele']:
                ds_name = 'INPUT'
                ann_value = self[annotation]
            else:
                ds_name = self.annotations[annotation].getDatasource()
                ann_value = self.annotations[annotation].getValue()
            output = (ds_name, annotation, ann_value)
            outputs.append(output)
        outputs.sort(key=lambda x: x[1])
        outputs.sort(key=lambda x: x[0])
        outputs.sort(key=lambda x: x[0] == 'INPUT', reverse=True)
        print "Datasource -- Annotation -- Annotation Value"
        for output in outputs:
            print "%s -- %s -- %s" % output

    def attributesEqual(self, other):
        """
        are the attributes of this MutationData equal to the Attributes of other
        :param other: a different MutationData
        :return: are the attributes of this MutationData equal to the Attributes of other
        """
        if self.attributes == other.attributes:
            for a in self.attributes:
                if self[a] != other[a]:
                    return False
            return True
        return False

    def positionStr(self):
        return "%s:%s-%s %s:%s" % (self.chr, self.start, self.end, self.ref_allele, self.alt_allele)

    def __setitem__(self, key, value):
        
        if key in MutationData.attributes:
            self.__dict__[key] = value
        else:
#            self.lock.acquire()
            if key not in self.annotations.keys():
                logging.getLogger(__name__).warn("Attempting to create an annotation using dictionary output.  Cannot determine annotation source, but creating it anyway.")
                self.annotations[key] = Annotation(value)
            else:
                self.annotations[key].value = value
#            self.lock.release()
        
    def __delitem__(self, key):
        if key in MutationData.attributes:
            logging.getLogger(__name__).warn('Ignoring attempt to delete a MutationData attribute: ' + key)
        else:
            del(self.annotations[key])
        
    def __getitem__(self, key):
        if key in MutationData.attributes:
            return self.__dict__[key]
        return self.annotations[key].value
    
    def __contains__(self, key):
        return key in self.annotations.keys()
    
    def __len__(self):
        return len(self.annotations)
    
    def __iter__(self):
        return iter(self.annotations)
    
    def __str__(self):
        return str(self.annotations)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            try:
                for k in self.keys():
                    if self.getAnnotation(k) != other.getAnnotation(k):
                        return False
            except KeyError:
                return False
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

