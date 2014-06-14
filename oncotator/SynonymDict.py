# LICENSE_GOES_HERE


"""
Created on Oct 22, 2012

@author: lichtens
"""
import collections
class SynonymDict(collections.MutableMapping):
    """
    If an item is added and there is no key for a synonym nor a real key, then both are created.
    """


    def __init__(self):
        """
        Constructor
        
        real keys are present in the synonym map as well.
        """
        self._synonymMap = dict()
        self.data = dict()
        
    def addSynonym(self, synKey, realKey):
        self._synonymMap[synKey] = realKey
        
    def __setitem__(self, key, value):
        
        if (key not in self._synonymMap.keys()) and (key not in self.data.keys()) :
            self._synonymMap[key] = key
        realKey = self._synonymMap[key]
        self.data[realKey] = value
        
    def __delitem__(self, key):
        del(self._synonymMap[key])
        
    def __getitem__(self, key):
        realKey = self._synonymMap[key]
        return self.data[realKey]
    
    def __contains__(self, key):
        return key in self._synonymMap.keys()

    def __len__(self):
        return len(self._synonymMap)
    
    def __iter__(self):
        return iter(self._synonymMap)
