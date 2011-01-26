#!/usr/bin/env python
# $URL$
# $Rev$
#
# v2index.py
#
# David Jones, Clear Climate Code, 2010-02-08

"""
Create an index of the input/v2.mean file.  This allows other programs
to have faster random access.
"""

import itertools
import sys

def build(inp, out):
    """Build an index of the v2.mean file `inp` and write it to the file
    `out`."""

    out.writelines(v2index(inp))

def v2index(inp):
    """Read the v2.mean file from `inp` and yield a series of index
    lines.
    """

    def id12(t):
        """The 12-digit record identifier from a (location, line)
        pair."""
        return t[1][:12]

    for id12, group in itertools.groupby(tell_stream(inp), id12):
        whence,first = list(group)[0]
        yield "%s %s %d\n" % (id12, first[12:16], whence)

def tell_stream(f):
    """For a input file `f` yield a stream of (location, line) pairs
    where each location is the value of the file pointer (as returned by
    f.tell()) at the beginning of the line."""

    while True:
        whence = f.tell()
        line = f.readline()
        if line == '':
            break
        yield whence, line

class File:
    """A v2.mean file accessible randomly using an index."""

    def __init__(self, file):
        """`file` can be an open file object, as long as it has the
        .name attribute, or a filename."""

        try:
            self.name = file.name
        except:
            self.name = file
        self.indexname = self.name + '.index'
        try:
            self.index = index(open(self.indexname))
        except IOError:
            self.build()
        try:
            file.name
            self.file = file
        except:
            self.file = open(self.name)

    def build(self):
        """(Re-) build the index file.  Updates self.index."""

        sys.stderr.write("Building index...\n")
        inp = open(self.name)
        out = open(self.indexname, 'w')
        build(inp, out)
        out.close()
        inp.close()
        sys.stderr.write("Done building index...\n")
        self.index = index(open(self.indexname))

    def get_id12(self, id12):
        """For a given 12-digit record identifier, return an iterator
        that yields each datum.
        """

        assert 12 == len(id12)

        # Get the first line according to the index, and check it
        rebuiltindex = False
        while True:
            i = self.index[id12]
            self.file.seek(i.whence)
            first = self.file.readline()
            if first[:16] != i.match:
                if rebuiltindex:
                    raise Error(
                      "Index still wrong after rebuilding.  id12=%s" % id12)
                sys.stderr.write("Index is wrong, rebuilding it...\n")
                rebuiltindex = True
                self.build()
                continue # goto, really
            break

        def iter():
            """Iterator for the lines."""

            line = first
            while line[:12] == id12:
                yield line
                line = self.file.readline()
        return iter()

    def get_id11(self, id11):
        """For a given 11-digit station identifier, return a sequence of
        (id12, series) pairs where each corresponds to a 12-digit record
        identifier, and each series yields the data for that record.
        """

        assert 11 == len(id11)

        for id12 in self.index.get(id11, []):
            yield (id12, self.get_id12(id12))

    def get(self, id):
        """As get_id11, but either 11-digit or 12-digit identifier is
        allowed.  If a 12-digit identifier is used then the sequence
        returned has only one pair.
        """

        if 11 == len(id):
            for pair in self.get_id11(id):
                yield pair
        else:
            assert 12 == len(id)
            yield (id, self.get_id12(id))

def index(f):
    """Read the index file `f` and return a dictionary that maps from
    id12 to an Index object, and also maps id11 to a list of id12s.
    """

    d = {}
    for line in f:
        i = Index(line)
        d[i.id12] = i
        # Append the id12 to a list associated with each id11, creating
        # a list if necessary.
        id11 = i.id12[:11]
        v = d.setdefault(id11, [])
        v.append(i.id12)
    return d

class Index:
    """A single entry from the index file."""

    def __init__(self, line):
        self.id12, self.year, self.whence = line.split()
        # Concatentate these two to form the 16 character string that
        # must match the v2.mean file at the specified place.
        self.match = self.id12 + self.year
        self.whence = int(self.whence)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    return build(open('input/v2.mean'), open('work/v2.mean.index', 'w'))

if __name__ == '__main__':
    main()
