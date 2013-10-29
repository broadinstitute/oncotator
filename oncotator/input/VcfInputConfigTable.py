class VcfInputConfigTable():
    _infoFieldIDs = dict()
    _formatFieldIDs = dict()
    _splitSet = dict()
    _notSplitSet = dict()

    def addInfoFieldID(self, ID, name):
        """
        Add info field ID and corresponding name.

        :param ID:
        :param name:
        """
        self._infoFieldIDs[ID] = name

    def addFormatFieldID(self, ID, name):
        """
        Add format field ID and corresponding name.

        :param ID: format field ID
        :param name: format field name
        """
        self._formatFieldIDs[ID] = name

    def removeInfoFieldID(self, ID):
        """

        :param ID:
        """
        if ID in self._infoFieldIDs:
            del self._infoFieldIDs[ID]

    def removeFormatFieldID(self, ID):
        """

        :param ID:
        """
        if ID in self._formatFieldIDs:
            del self._formatFieldIDs[ID]

    def getInfoFieldIDs(self):
        """


        :return:
        """
        return self._infoFieldIDs.keys()

    def getFormatFieldIDs(self):
        """


        :return:
        """
        return self._formatFieldIDs.keys()

    def getInfoFieldName(self, ID):
        """

        :param ID:
        :return:
        """
        if ID in self._infoFieldIDs:
            return self._infoFieldIDs[ID]
        return ID

    def getFormatFieldName(self, ID):
        """

        :param ID:
        :return:
        """
        if ID in self._formatFieldIDs:
            return self._formatFieldIDs[ID]
        return ID

    def addFieldIDsToNotSplitSet(self, fieldType, IDs):
        """

        :param fieldType:
        :param IDs:
        """
        self._notSplitSet[fieldType] = IDs

    def addFieldIDsToSplitSet(self, fieldType, IDs):
        """

        :param fieldType:
        :param IDs:
        """
        self._splitSet[fieldType] = IDs

    def isFieldIDInSplitSet(self, fieldType, ID):
        """

        :param fieldType:
        :param ID:
        :return:
        """
        if fieldType in self._splitSet and ID in self._splitSet[fieldType]:
            return True
        return fieldType in self._splitSet and ID in self._splitSet[fieldType]

    def isFieldIDInNotSplitSet(self, fieldType, ID):
        """

        :param fieldType:
        :param ID:
        :return:
        """
        if fieldType in self._notSplitSet and ID in self._notSplitSet[fieldType]:
            return True
        return False