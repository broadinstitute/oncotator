from oncotator.Annotation import Annotation


class VcfOutputAnnotation(Annotation):

    def __init__(self, ID, fieldType="FORMAT", split=False, datasourceName="INPUT", dataType="String",
                 description="Unknown", number=None):
        """


        :param ID: unique identifier
        :param fieldType:
        :param split:
        :param datasourceName:
        :param dataType:
        :param description:
        :param number:
        """
        Annotation.__init__(self, None, datasourceName, dataType, description, [], number)
        self.split = split
        self.fieldType = fieldType
        self.ID = ID

    def setID(self, ID):
        """

        :param ID:
        """
        self.ID = ID

    def getID(self):
        """


        :return:
        """
        return self.ID

    def setIsSplit(self, split):
        """

        :param split:
        """
        self.split = split

    def isSplit(self):
        """


        :return:
        """
        return self.split

    def setFieldType(self, ft):
        """

        :param ft:
        """
        self.fieldType = ft

    def getFieldType(self):
        """


        :return: field type
        """
        return self.fieldType