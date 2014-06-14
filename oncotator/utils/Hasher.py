# LICENSE_GOES_HERE
import hashlib
import logging


class Hasher(object):

    def __init__(self):
        self._md5hash = hashlib.md5()

    @staticmethod
    def md5_hash(input_string):
        md5hash = hashlib.md5()
        md5hash.update(input_string)
        return md5hash.hexdigest()

    def update(self, data):
        self._md5hash.update(data)

    def reset(self):
        self._md5hash = hashlib.md5()

    def hexdigest(self):
        return self._md5hash.hexdigest()

    def _read_md5_file(self, filepath):
        f1 = open(filepath, 'r')
        buf = f1.read(4096)
        f1.close()
        return buf

    def _create_md5_file(self, filepath):
        file_md5_hash = hashlib.md5()
        f1 = file(filepath, 'r')
        while 1:
            # Read file in as little chunks
            buf = f1.read(4096)
            if not buf : break
            file_md5_hash.update(hashlib.md5(buf).hexdigest())
        f1.close()
        outFP = file(filepath + ".md5", 'w')
        outFP.write(file_md5_hash.hexdigest())
        outFP.close()


    def _write_hex_to_md5_file(self, filepath, hex):
        f1 = file(filepath, 'w')
        f1.write(hex)
        f1.close()


    def _compute_md5(self, filepath):
        file_md5_hash = hashlib.md5()
        f1 = file(filepath, 'r')
        while 1:
            # Read file in as little chunks
            buf = f1.read(4096)
            if not buf : break
            file_md5_hash.update(hashlib.md5(buf).hexdigest())
        f1.close()
        return file_md5_hash.hexdigest()

    def create_hashcode_for_dir(self, directory):
        """Create a single md5 representing every file in every subdirectory for this directory. """
        import hashlib, os
        md5hash = hashlib.md5()

        if not os.path.exists(directory):
            logging.getLogger(__name__).error(directory + " not found.")
            raise ValueError(directory + " not found.")
        md5s = []
        for root, dirs, files in os.walk(directory):

            for names in files:
                if names.endswith(".old") or names.endswith("bkup") or names.endswith(".md5") or names.startswith(".") or names.endswith("~"):
                    logging.getLogger(__name__).info("Not creating new hash for " + names)
                    continue

                file_to_hash = names
                if os.path.exists(os.path.join(root, names + ".md5")):
                    logging.getLogger(__name__).info("using pre-existing hash for " + names)
                    file_to_hash = names + ".md5"

                filepath = os.path.join(root, file_to_hash)

                if not filepath.endswith(".md5"):
                    logging.getLogger(__name__).info("Computing new hash for " + filepath)
                    md5 = self._compute_md5(filepath)
                else:
                    md5 = self._read_md5_file(filepath)

                md5s.append(md5)

        # We have to sort the md5s, since the order can change when walking the directory.
        sorted_md5s = sorted(md5s)
        for md5 in sorted_md5s:
            md5hash.update(md5)
        return md5hash.hexdigest()

