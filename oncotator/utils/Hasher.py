import hashlib


class Hasher(object):

    def __init__(self):
        self._md5hash = hashlib.md5()

    @staticmethod
    def md5_hash(self, input_string):
        md5hash = hashlib.md5()
        md5hash.update(input_string)
        return md5hash.hexdigest()

    def update(self, data):
        self._md5hash.update(data)

    def reset(self):
        self._md5hash = hashlib.md5()

    def hexdigest(self):
        return self._md5hash.hexdigest()


