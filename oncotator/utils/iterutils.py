import itertools

__author__ = 'louisb'

from itertools import chain, imap

def flatmap(f, items):
    return list(chain.from_iterable(imap(f, items)))

def chop(iterable, length=2):
    return itertools.izip(*(iter(iterable) , ) *length)