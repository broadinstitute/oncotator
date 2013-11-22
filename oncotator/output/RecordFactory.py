import collections
import logging
import string
import vcf
from oncotator.utils.MutUtils import MutUtils


class RecordFactory:
    def __init__(self, chrom, pos, ref, sampleNames):
        """

        :param chrom:
        :param pos:
        :param ref:
        :param sampleNames:
        """
        self._logger = logging.getLogger(__name__)
        self._chrom = chrom
        self._pos = pos
        self._ID = None  # semi-colon separated list of unique identifiers where available
        self._refAllele = ref  # reference base(s)
        self._alts = []  # comma separated list of alternate non-reference alleles called on at least one of the samples
        self._qual = None  # phred-scaled quality score for the assertion made in ALT
        self._filt = []  # PASS if this position has passed all filters
        self._info = collections.OrderedDict()  # additional information
        self._infoFieldProperty = collections.OrderedDict()
        self._fmt = [None]*len(sampleNames)
        self._fmtIDs = []
        self._fmtFieldProperty = collections.OrderedDict()
        self._sampleNames = sampleNames
        self._sampleNameIndexes = dict([(x, i) for (i, x) in enumerate(sampleNames)])
        self._fieldProperty = collections.namedtuple(typename="Property", field_names=["num", "dataType", "isSplit"])

    def _map(self, func, iterable, bad=(".", "",)):
        """

        :param func:
        :param iterable:
        :param bad:
        :return:
        """
        return [func(v) if v not in bad else None for v in iterable]

    def _replaceChrs(self, text, frm, to):
        """

        :param text:
        :param frm:
        :param to:
        :return:
        """
        tbl = string.maketrans(frm, to)
        return text.translate(tbl)

    def _resolveInfo(self):
        """


        :return:
        """
        IDs = self._info.keys()
        info = collections.OrderedDict()

        nalts = len(self._alts)
        nsamples = len(self._sampleNames)
        for ID in IDs:
            val = self._info[ID]
            prop = self._infoFieldProperty[ID]

            if prop.num == -2:
                if len(val) == nsamples and len(filter(None, val)) != 0:
                    info[ID] = val
            elif prop.num == -1:
                if len(val) == nalts and len(filter(None, val)) != 0:
                    info[ID] = val
            elif prop.num == 0:
                if val:
                    info[ID] = val
            elif prop.num is None:
                if prop.isSplit:
                    if len(val) == nalts and len(filter(None, val)) != 0:
                        info[ID] = val
                elif len(filter(None, val)) != 0:
                    info[ID] = val
            else:
                if len(val) == prop.num and len(filter(None, val)) != 0:
                    info[ID] = val

        return info

    def _resolveSamples(self, record):
        """

        :param record:
        :return:
        """
        nalts = len(self._alts)
        samples = [None]*len(self._sampleNames)

        IDs = [None]*len(self._fmtIDs)
        dataTypes = [None]*len(self._fmtFieldProperty)
        nums = [None]*len(self._fmtFieldProperty)

        for sampleName in self._sampleNames:
            sampleData = [None]*len(self._fmtIDs)
            sampleNameIndex = self._sampleNameIndexes[sampleName]
            data = self._fmt[sampleNameIndex]
            for i in xrange(len(self._fmtIDs)):
                ID = self._fmtIDs[i]
                IDs[i] = self._fmtIDs[i]
                prop = self._fmtFieldProperty[ID]
                dataTypes[i] = prop.dataType
                nums[i] = prop.num if prop.num is not None else "."
                val = [None]

                if (data is not None) and (ID in data):
                    val = data[ID]
                    if prop.num == -2:
                        pass
                    elif prop.num == -1:
                        if len(val) != nalts:
                            val = nalts*[None]
                    elif prop.num == 0:
                        pass
                    elif prop.num is None:
                        if prop.isSplit:
                            if len(val) != nalts:
                                val = nalts*[None]
                    else:
                        if len(val) != prop.num:
                            val = abs(prop.num)*[None]

                if ID == "GT":
                    sampleData[i] = val[0]
                else:
                    sampleData[i] = val

            calldata = vcf.model.make_calldata_tuple(IDs)
            calldata._types = dataTypes
            calldata._nums = nums
            samples[self._sampleNameIndexes[sampleName]] = calldata(*sampleData)

        for sampleName in self._sampleNames:
            samples[self._sampleNameIndexes[sampleName]] = \
                vcf.model._Call(record, sampleName, samples[self._sampleNameIndexes[sampleName]])
        return samples

    def createRecord(self):
        """


        :return:
        """
        chrom = self._chrom
        pos = self._pos
        refAllele = self._refAllele
        alts = self._alts
        if len(alts) == 0:
            alts = ["."]

        ID = self._ID
        if ID is not None:
            ID = string.join(ID, ";")

        qual = self._qual
        if qual is None:
            self._logger.warn("Variant at chromosome %s and position %s is missing phred-scaled quality score "
                              "(typically: annotation 'qual')." % (chrom, pos))

        filt = self._filt
        info = self._resolveInfo()
        fmt = string.join(self._fmtIDs, ":")
        record = vcf.model._Record(chrom, pos, ID, refAllele, alts, qual, filt, info, fmt, self._sampleNameIndexes)
        record.samples = self._resolveSamples(record)
        return record

    def _correct(self, iterable, bad=(".", "",)):
        """

        :param iterable:
        :param bad:
        :return:
        """
        return [v if v not in bad else None for v in iterable]

    def _fixVal(self, val, isSplit):
        """

        :param val:
        :param isSplit:
        :return:
        """
        if isSplit:
            val = self._replaceChrs(val, ",=;\n\t ", "|~|#__")  # exclude ":"
        else:
            val = self._replaceChrs(val, "=;\n\t :", "~|#__>")  # exclude ":" and ","

        if not isSplit:
            val = self._correct(val.split(","))
        else:
            val = self._correct([val])

        return val

    def _appendVal2FixedNumField(self, data, ID, num, isSplit, val):
        """

        :param data:
        :param ID:
        :param num:
        :param isSplit:
        :param val:
        """
        if ID not in data:
            data[ID] = self._fixVal(val, isSplit)
        elif isSplit and num > 1:
            vals = data[ID]
            if len(vals) < num:
                vals += self._fixVal(val, isSplit)
                data[ID] = vals

    def addInfo(self, sampleName, ID, num=".", dataType="String", val=None, isSplit=True):
        """

        :param sampleName:
        :param ID:
        :param num:
        :param dataType:
        :param val:
        :param isSplit:
        """
        if num == -2:  # num is the number of samples
            nsamples = len(self._sampleNames)
            if sampleName in self._sampleNames:
                if ID not in self._info:
                    self._info[ID] = [None]*nsamples
                sampleNameIndex = self._sampleNameIndexes[sampleName]
                val = self._fixVal(val, isSplit)
                self._info[ID][sampleNameIndex] = val[0]
        elif num == -1:  # num is the number of alternate alleles
            nalts = len(self._alts)
            self._appendVal2FixedNumField(self._info, ID, nalts, isSplit, val)
        elif num == 0:  # num is either true or false
            if ID not in self._info:
                val = self._map(MutUtils.str2bool, self._fixVal(val, isSplit))  # convert the value to boolean
                self._info[ID] = val[0]
        elif num is None:  # num is unknown
            nalts = len(self._alts)
            self._appendVal2FixedNumField(self._info, ID, nalts, isSplit, val)
        else:
            self._appendVal2FixedNumField(self._info, ID, num, isSplit, val)

        if ID not in self._infoFieldProperty:
            self._infoFieldProperty[ID] = self._fieldProperty(num, dataType, isSplit)

    def addFormat(self, sampleName, ID, num=".", dataType="String", val=None, isSplit=True):
        """

        :param sampleName:
        :param ID:
        :param num:
        :param dataType:
        :param val:
        :param isSplit:
        """
        if sampleName in self._sampleNames and num != 0:
            sampleNameIndex = self._sampleNameIndexes[sampleName]
            if self._fmt[sampleNameIndex] is None:
                self._fmt[sampleNameIndex] = collections.OrderedDict()
                self._fmtIDs = ["GT"]
                self._fmtFieldProperty["GT"] = self._fieldProperty(1, "String", False)

            if num == -2:  # num is the number of samples
                nalts = len(self._alts)
                self._appendVal2FixedNumField(self._fmt[sampleNameIndex], ID, nalts, isSplit, val)
            elif num == -1:  # num is the number of alternate alleles
                nalts = len(self._alts)
                self._appendVal2FixedNumField(self._fmt[sampleNameIndex], ID, nalts, isSplit, val)
            elif num == 0:  # num is either true or false
                pass
            elif num is None:  # num is unknown
                nalts = len(self._alts)
                self._appendVal2FixedNumField(self._fmt[sampleNameIndex], ID, nalts, isSplit, val)
            else:
                self._appendVal2FixedNumField(self._fmt[sampleNameIndex], ID, num, isSplit, val)

            if ID not in self._fmtIDs:
                self._fmtIDs += [ID]

            if ID not in self._fmtFieldProperty:
                self._fmtFieldProperty[ID] = self._fieldProperty(num, dataType, isSplit)

    def addQual(self, qual):
        try:
            self._qual = int(qual)
        except ValueError:
            try:
                self._qual = float(qual)
            except ValueError:
                self._qual = None

    def addID(self, ID):
        if ID not in (".", "",):
            if self._ID is None:
                self._ID = [ID]
            else:
                if ID not in self._ID:
                    self._ID += [ID]

    def addAlt(self, alt):
        if alt in ("",):
            alt = "."
        if alt not in self._alts:
            self._alts += [alt]

    def addFilter(self, filt, val):
        if val not in ("PASS", ".",):
            if filt not in self._filt:
                self._filt += [filt]

    def setChrom(self, chrom):
        self._chrom = chrom

    def setPos(self, pos):
        self._pos = pos

    def setReferenceAllele(self, ref):
        self._refAllele = ref

    def setSampleNames(self, sampleNames):
        self._sampleNames = sampleNames
        self._sampleNameIndex = dict([(x, i) for (i, x) in enumerate(sampleNames)])