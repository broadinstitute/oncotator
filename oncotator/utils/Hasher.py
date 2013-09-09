import hashlib


class Hasher(object):

    def __init__(self):
        pass

    @staticmethod
    def md5_hash(self, input_string):
        md5hash = hashlib.md5()
        md5hash.update(input_string)
        return md5hash.hexdigest()
