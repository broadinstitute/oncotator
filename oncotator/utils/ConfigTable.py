
class ConfigTable(object):
    _filtDesc = None
    _info = None
    _infoDesc = None
    _fmt = None
    _fmtDesc = None
    _split = None
    _notSplit = None

    def __init__(self):
        self._filtDesc = dict()
        self._info = dict()
        self._infoDesc = dict()
        self._fmt = dict()
        self._fmtDesc = dict()
        self._split = dict()
        self._notSplit = dict()

    def addInfoFieldID(self, ID, name):
        if ID not in self._info:
            self._info[ID] = name

    def addInfoFieldIDDesc(self, ID, desc):
        if ID not in self._info:
            self._infoDesc[ID] = desc

    def removeInfoFieldID(self, ID):
        if ID in self._info:
            del self._info[ID]
        if ID in self._infoDesc:
            del self._infoDesc[ID]

    def addFormatFieldID(self, ID, name):
        if ID not in self._fmt:
            self._fmt[ID] = name

    def addFormatFieldIDDesc(self, ID, desc):
        if ID not in self._info:
            self._fmtDesc[ID] = desc

    def removeFormatFieldID(self, ID):
        if ID in self._fmt:
            del self._fmt[ID]

    def splitFieldID(self, ID, fieldType):
        if fieldType not in self._split:
            return False
        if fieldType in self._split:
            if ID in self._split[fieldType]:
                return True
        return False

    def addFieldIDToSplit(self, fieldType, ID):
        if fieldType not in self._split:
            self._split[fieldType] = set(ID)
        self._split[fieldType].add(ID)

    def addFieldIDToNotSplit(self, fieldType, ID):
        if fieldType not in self._notSplit:
            self._notSplit[fieldType] = set(ID)
        self._notSplit[fieldType].add(ID)

    def getFormatFieldIDs(self):
        return self._fmt.keys()

    def getFormatFieldName(self, ID):
        if ID not in self._fmt:
            return ID
        return self._fmt[ID]

    def getInfoFieldIDs(self):
        return self._info.keys()

    def getInfoFieldName(self, ID):
        if ID not in self._info:
            return ID
        return self._info[ID]