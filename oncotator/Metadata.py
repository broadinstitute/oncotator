# LICENSE_GOES_HERE

from collections import OrderedDict
from oncotator.Annotation import Annotation
from oncotator.MutationData import MutationData

__author__ = 'lichtens'


class Metadata(OrderedDict):
    """OrderedDict of key:Annotation pairs"""

    def asDict(self):
        result = OrderedDict()
        ks = self.keys()
        for k in ks:
            if isinstance(self[k], Annotation):
                result[k] = self[k].value
            else:
                result[k] = self[k]
        return result