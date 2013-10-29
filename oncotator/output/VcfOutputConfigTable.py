
class VcfOutputConfigTable():
    _infoFieldNames = dict()  # info (name, ID) pairs
    _formatFieldNames = dict()  # format (name, ID) pairs
    _filterFieldNames = dict()
    _otrFieldNames = dict()  # ID, QUAL and FILTER (name, ID) pairs

    _infoFieldNamesDescriptions = dict()  # info (name, description) pairs
    _formatFieldNamesDescriptions = dict()  # format (name, description) pairs
    _filterFieldNamesDescriptions = dict()  # filter (name, description) pairs

    _splitSet = dict()  #
    _notSplitSet = dict()  #

    def addInfoFieldName(self, name, ID):
        self._infoFieldNames[name] = ID

    def addFormatFieldName(self, name, ID):
        self._formatFieldNames[name] = ID

    def addOtherFieldName(self, name, ID):
        self._otrFieldNames[name] = ID
        if ID == "FILTER":
            self._filterFieldNames[name] = name

    def removeInfoFieldName(self, name):
        if name in self._infoFieldNames:
            del self._infoFieldNames[name]

    def removeFormatFieldName(self, name):
        if name in self._formatFieldNames:
            del self._formatFieldNames[name]

    def removeOtherFieldName(self, name):
        if name in self._otrFieldNames:
            del self._otrFieldNames[name]

    def addInfoFieldNameDescription(self, name, description):
        self._infoFieldNamesDescriptions[name] = description

    def addFormatFieldNameDescription(self, name, description):
        self._formatFieldNamesDescriptions[name] = description

    def addFilterFieldNameDescription(self, name, description):
        self._filterFieldNamesDescriptions[name] = description

    def getInfoFieldNames(self):
        return self._infoFieldNames.keys()

    def getFormatFieldNames(self):
        return self._formatFieldNames.keys()

    def getFilterFieldNames(self):
        return self._filterFieldNames.keys()

    def getOtherFieldNames(self):
        return self._otrFieldNames.keys()

    def getInfoFieldID(self, name):
        if name in self._infoFieldNames:
            return self._infoFieldNames[name]
        return name

    def getFormatFieldID(self, name):
        if name in self._formatFieldNames:
            return self._formatFieldNames[name]
        return name

    def getOtherFieldID(self, name):
        if name in self._otrFieldNames:
            return self._otrFieldNames[name]
        return name

    def getInfoFieldNameDescription(self, name):
        if name in self._infoFieldNamesDescriptions:
            return self._infoFieldNamesDescriptions[name]
        return "Unknown"

    def getFormatFieldNameDescription(self, name):
        if name in self._formatFieldNamesDescriptions:
            return self._formatFieldNamesDescriptions[name]
        return "Unknown"

    def getFilterFieldNameDescription(self, name):
        if name in self._filterFieldNamesDescriptions:
            return self._filterFieldNamesDescriptions[name]
        return "Unknown"

    def addFieldNamesToNotSplitSet(self, fieldType, names):
        self._notSplitSet[fieldType] = names

    def addFieldNamesToSplitSet(self, fieldType, names):
        self._splitSet[fieldType] = names

    def isFieldNameInSplitSet(self, fieldType, name):
        if fieldType in self._splitSet and name in self._splitSet[fieldType]:
            return True
        return fieldType in self._splitSet and name in self._splitSet[fieldType]

    def isFieldNameInNotSplitSet(self, fieldType, name):
        if fieldType in self._notSplitSet and name in self._notSplitSet[fieldType]:
            return True
        return False