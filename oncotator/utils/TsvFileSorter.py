# LICENSE_GOES_HERE


import collections
import heapq
import itertools
from oncotator.utils.GenericTsvReader import GenericTsvReader
import operator
import string
from oncotator.utils.CallbackException import CallbackException
from oncotator.utils.MutUtils import MutUtils


class TsvFileSorter(object):

    def __init__(self, filename, commentPrepend="#", delimiter="\t", lineterminator="\n"):
        """
        A class for sorting a delimited file. A merge sort procedure is used to sort the file.

        :type self: object
        :param filename: file that needs to be sorted
        :param commentPrepend: character that indicates beginning of a comment line
        :param delimiter: separates column values
        :param lineterminator: character that signifies the end of line
        """
        self.readfilename = filename
        self._Pair = collections.namedtuple(typename="Pair", field_names=["key", "value"])
        self.delimiter = delimiter
        self.lineterminator = lineterminator
        self.commentPrepend = commentPrepend

    def _merge(self, partitions):
        """
        This method uses heap queue to merge sorted partitions of Pair tuples. For each Pair tuple in merged partitions,
        it writes out the value.

        :param partitions: collection of Pair tuples
        """
        for pair in heapq.merge(*partitions):
            yield pair.value

    def _yieldPartitions(self, iterable, func, fieldnameIndexes, length):
        """
        This method parses a set of lines for a partition, applies an anonymous function that converts each line of the
        partition to a key-value pair where key is of type tuple, sorts the key-value pairs on the keys and then yields
        the partition. Through this method, we obtain several sorted chunks.

        :param iterable: lines of text
        :param func: function that converts each row of the input file to an unique key
        :param fieldnameIndexes: dictionary of fieldnames and corresponding indexes
        :param length: determines the number of lines in the buffer
        """
        isKeyTuple = False

        # Take the first "length" number of items and return them as list.
        lines = list(itertools.islice(iterable, length))
        data = collections.OrderedDict()

        while len(lines) > 0:
            pairs = [None]*len(lines)

            # Create a list of (key, value) pairs
            # Each key consists of a tuple, value is the corresponding text
            # Note: CSV dictionary reader is NOT used because a chunk of text is parsed at a time, rather than a
            # line of text
            for i in xrange(len(lines)):
                # Note: CSV dictionary reader is NOT used because a chunk of text is parsed at a time, rather than a
                # line of text
                line = lines[i]
                tokens = MutUtils.getTokens(line, self.delimiter, self.lineterminator)

                for fieldname, index in fieldnameIndexes.items():
                    data[fieldname] = tokens[index]

                key = func(data)

                if not isKeyTuple:
                    isKeyTuple = isinstance(key, tuple)
                    if not isKeyTuple:
                        raise CallbackException("The value returned by the callback must be a tuple. Instead, a value "
                                                "of %s was returned." % (type(key)))
                pairs[i] = self._Pair(key, line)

            partition = sorted(pairs, key=operator.attrgetter("key"))

            lines = list(itertools.islice(iterable, length))

            yield partition

    def sortFile(self, filename, func, length=50000):
        """
        This method sorts the input file and writes out the sorted file to filename.

        :param filename: sorted filename
        :param func: function that converts each row of the input file to an unique, sortable key
        :param length: maximum number of lines in a partition
        """
        reader = GenericTsvReader(filename=self.readfilename, commentPrepend=self.commentPrepend,
                                  delimiter=self.delimiter)
        comments = reader.getComments()

        fieldnames = reader.getFieldNames()
        if fieldnames is None:
            fieldnames = []

        fieldnameIndexes = collections.OrderedDict()
        if fieldnames is not None:
            fieldnameIndexes = collections.OrderedDict([(x, i) for (i, x) in enumerate(fieldnames)])

        iterable = iter(reader.getInputContentFP())
        partitions = self._yieldPartitions(iterable, func, fieldnameIndexes, length)

        with open(name=filename, mode='wb', buffering=64 * 1024) as writer:
            writer.write(comments)
            writer.write(string.join(fieldnames, self.delimiter) + "\n")
            writer.writelines(self._merge(partitions))  # generators are allowed as inputs to writelines function