# LICENSE_GOES_HERE


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
            if annotation in ['build', 'chr' , 'start', 'end', 'ref_allele', 'alt_allele']:
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
