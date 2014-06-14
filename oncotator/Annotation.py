"""
Created on Nov 20, 2012

@author: lichtens
"""
# LICENSE_GOES_HERE


class Annotation(object):
    """
    Simple class for storing annotations (particularly value and datasource)
    value should be a string
    datasource should be a string
    """

    def __init__(self, value, datasourceName="Unknown", dataType="String", description="", tags=list(), number=None):
        """
        
        """
        self.value = value
        self.datasourceName = datasourceName
        self.dataType = dataType
        self.description = description
        self.tags = tags
        self.number = number

    def setValue(self, val):
        self.value = val
    
    def setDatasource(self, ds):
        self.datasourceName = ds
    
    def setDataType(self, dt):
        self.dataType = dt
    
    def setDescription(self, desc):
        self.description = desc
    
    def setNumber(self, num):
        self.number = num
    
    def getValue(self):
        return self.value
    
    def getDatasource(self):
        return self.datasourceName

    def getDataType(self):
        return self.dataType

    def getDescription(self):
        return self.description
    
    def getNumber(self):
        return self.number
    
    def addTag(self, tag):
        self.tags.append(tag)
        
    def getTags(self):
        return self.tags

    def isEqual(self, annot):
        if self.value != annot.getValue():
            return False
        if self.datasourceName != annot.getDatasource():
            return False
        if self.dataType != annot.getDataType():
            return False
        if self.description != annot.getDescription():
            return False
        if self.number != annot.getNumber():
            return False
        if len(set(self.tags).symmetric_difference(set(annot.getTags()))) != 0:
            return False
        return True