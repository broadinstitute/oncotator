
class ConfigTable(object):
    _filt = None
    _filtDesc = None
    _info = None
    _infoDesc = None
    _fmt = None
    _fmtDesc = None
    _other = None
    _split = None
    _notSplit = None

    def __init__(self):
        self._filt = dict()
        self._filtDesc = dict()
        self._info = dict()
        self._infoDesc = dict()
        self._fmt = dict()
        self._fmtDesc = dict()
        self._other = dict()
        self._split = dict()
        self._notSplit = dict()

    def addInfoFieldID(self, ID, name):
        self._info[ID] = name

    def addFormatFieldID(self, ID, name):
        self._fmt[ID] = name

    def addFilterFieldID(self, ID, name):
        self._filt[ID] = name

    def addOtherFieldID(self, ID, name):
        self._other[ID] = name

    def addInfoFieldIDDesc(self, ID, desc):
        self._infoDesc[ID] = desc

    def addFormatFieldIDDesc(self, ID, desc):
        self._fmtDesc[ID] = desc

    def addFilterFieldIDDesc(self, ID, desc):
        self._filtDesc[ID] = desc

    def removeInfoFieldID(self, ID):
        if ID in self._info:
            del self._info[ID]
        if ID in self._infoDesc:
            del self._infoDesc[ID]

    def removeFormatFieldID(self, ID):
        if ID in self._fmt:
            del self._fmt[ID]
        if ID in self._fmtDesc:
            del self._fmtDesc[ID]

    def removeFilterFieldID(self, ID):
        if ID in self._filt:
            del self._filt[ID]
        if ID in self._filtDesc:
            del self._filtDesc[ID]

    def getInfoFieldIDs(self):
        return self._info.keys()

    def getFormatFieldIDs(self):
        return self._fmt.keys()

    def getFilterFieldIDs(self):
        return self._filt.keys()

    def getOtherFieldIDs(self):
        return self._other.keys()

    def getFormatFieldName(self, ID):
        if ID not in self._fmt:
            return ID
        return self._fmt[ID]

    def getInfoFieldName(self, ID):
        if ID not in self._info:
            return ID
        return self._info[ID]

    def getOtherFieldName(self, ID):
        if ID not in self._other:
            return ID
        return self._other[ID]

    def isFieldSplit(self, ID, fieldType):
        if fieldType not in self._split:
            return False
        if fieldType in self._split:
            if ID in self._split[fieldType]:
                return True
        return False

    def isFieldNotSplit(self, ID, fieldType):
        if fieldType not in self._notSplit:
            return False
        if fieldType in self._notSplit:
            if ID in self._notSplit[fieldType]:
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

    def getInfoFieldIDDesc(self, ID):
        if ID in self._infoDesc:
            return self._infoDesc[ID]
        return "Unknown"

    def getFormatFieldIDDesc(self, ID):
        if ID in self._fmtDesc:
            return self._fmtDesc[ID]
        return "Unknown"

    def getInfoFilterIDDesc(self, ID):
        if ID in self._filtDesc:
            return self._filtDesc[ID]
        return "Unknown"