# LICENSE_GOES_HERE


class VcfInputConfigTable():
    """
    Container class that handles parsed input VCF configuration data.
    """

    _infoFieldIDs = dict()
    _formatFieldIDs = dict()

    _splitSet = dict()
    _notSplitSet = dict()

    def addInfoFieldID(self, ID, name):
        self._infoFieldIDs[ID] = name

    def addFormatFieldID(self, ID, name):
        self._formatFieldIDs[ID] = name

    def removeInfoFieldID(self, ID):
        if ID in self._infoFieldIDs:
            del self._infoFieldIDs[ID]

    def removeFormatFieldID(self, ID):
        if ID in self._formatFieldIDs:
            del self._formatFieldIDs[ID]

    def getInfoFieldIDs(self):
        return self._infoFieldIDs.keys()

    def getFormatFieldIDs(self):
        return self._formatFieldIDs.keys()

    def getInfoFieldName(self, ID):
        if ID in self._infoFieldIDs:
            return self._infoFieldIDs[ID]
        return ID

    def getFormatFieldName(self, ID):
        if ID in self._formatFieldIDs:
            return self._formatFieldIDs[ID]
        return ID

    def addFieldIDsToNotSplitSet(self, fieldType, IDs):
        self._notSplitSet[fieldType] = IDs

    def addFieldIDsToSplitSet(self, fieldType, IDs):
        self._splitSet[fieldType] = IDs

    def isFieldIDInSplitSet(self, fieldType, ID):
        if fieldType in self._splitSet and ID in self._splitSet[fieldType]:
            return True
        return fieldType in self._splitSet and ID in self._splitSet[fieldType]

    def isFieldIDInNotSplitSet(self, fieldType, ID):
        if fieldType in self._notSplitSet and ID in self._notSplitSet[fieldType]:
            return True
        return False