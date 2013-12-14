from oncotator.Annotation import Annotation


class VcfOutputAnnotation(Annotation):

    def __init__(self, ID, fieldType="FORMAT", split=False, datasourceName="INPUT", dataType="String",
                 description="Unknown", number=None):
        Annotation.__init__(self, None, datasourceName, dataType, description, [], number)
        self.split = split
        self.fieldType = fieldType
        self.ID = ID

    def setID(self, ID):
        self.ID = ID

    def getID(self):
        return self.ID

    def setIsSplit(self, split):
        self.split = split

    def isSplit(self):
        return self.split

    def setFieldType(self, ft):
        self.fieldType = ft

    def getFieldType(self):
        return self.fieldType