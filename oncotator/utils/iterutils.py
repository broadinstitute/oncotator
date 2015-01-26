__author__ = 'louisb'

from itertools import chain, imap

def flatmap(f, items):
    return list(chain.from_iterable(imap(f, items)))