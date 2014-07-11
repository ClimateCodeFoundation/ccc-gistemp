#!/usr/bin/env python
# $URL$
# $Rev$
#
# v2index.py
#
# David Jones, Clear Climate Code, 2010-02-08

"""
Create an index of a file in GHCN-M format (typically either
input/ghcnm.tavg.qca.dat (GHCN-M v3) or input/v2.mean (GHCN-M
v2)).  This allows other programs to have faster random access.
"""

import itertools
import sys

class Error(Exception):
    pass

def build(inp, out):
    """Build an index of the GHCN-M file `inp` and write it to the file
    `out`."""

    index = None
    if 'v2' in inp.name:
        index = v2_index
    if inp.name.endswith('.dat'):
        index = v3_index

    if index is None:
        raise Error("Can't tell if input is GHCN-M v2 or v3")

    out.writelines(index(inp))

def v2_index(inp):
    return index_grouping(inp, lambda line: line[:12])

def v3_index(inp):
    return index_grouping(inp, lambda line: line[:11])

def index_grouping(inp, id):
    """Read the GHCN-M v2 file from `inp` and yield a series of index
    lines.
    """

    def id_of_pair(t):
        """The record identifier from a (location, line) pair."""
        return id(t[1])

    for ghcn_id, group in itertools.groupby(tell_stream(inp),
      id_of_pair):
        whence,first = list(group)[0]
        year = first[len(ghcn_id):][:4]
        yield "%s %s %d\n" % (ghcn_id, year, whence)

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
    """A GHCN-M file accessible randomly using an index."""

    def __init__(self, file):
        """
        `file` can be an open file object, as long as it has the
        .name attribute, or a filename.
        """

        try:
            self.name = file.name
            self.file = file
        except AttributeError:
            self.name = file
            self.file = open(self.name)
        self.index_name = self.name + '.index'
        try:
            self.index = index(open(self.index_name))
        except IOError:
            self.build()

    def build(self):
        """
        (Re-) build the index file.  Updates self.index.
        """

        sys.stderr.write("Building index...\n")
        inp = open(self.name)
        out = open(self.index_name, 'w')
        build(inp, out)
        out.close()
        inp.close()
        sys.stderr.write("Done building index...\n")
        self.index = index(open(self.index_name))

    def get_id12(self, id12):
        """For a given 12-digit record identifier, return an iterator
        that yields each datum.
        """

        assert 12 == len(id12)

        return get_single_id(self, id12)

    def get_single_id(self, id):
        """
        For a given record identifier (11-digit in GHCN-M v3,
        12-digit in GHCN-M v2), return an iterator that yields
        each datum.
        """

        # Get the first line according to the index, and check it.
        rebuilt_index = False
        while True:
            i = self.index[id]
            self.file.seek(i.whence)
            first_line = self.file.readline()
            if first_line.startswith(i.match):
                break
            if rebuilt_index:
                raise Error(
                  "Index still wrong after rebuilding.  id=%s" % id)
            sys.stderr.write("Index is wrong, rebuilding it...\n")
            rebuilt_index = True
            self.build()

        def iter():
            """Iterator for the lines."""

            # This is a pretty naughty capture of `first_line`
            # from the enclosing function.
            line = first_line
            while line.startswith(id):
                yield line
                line = self.file.readline()
        return iter()

    def get_many_id(self, id11):
        """
        For an identifier that maps to many records (this must
        be an 11-digit station identifier for a GHCN-M v2 file),
        return a sequence of (id12, series) pairs where each id12
        corresponds to a 12-digit record identifier, and each series
        yields the data for that record.
        """

        assert 11 == len(id11)

        for id12 in self.index.get(id11, []):
            yield (id12, self.get_id12(id12))

    def get(self, id):
        """
        As get_id11, but either 11-digit or 12-digit identifier is
        allowed.
        """

        item = self.index.get(id)
        if not item:
            return

        if isinstance(item, Index):
            yield (id, self.get_single_id(id))
            return

        for pair in self.get_many_id(id):
            yield pair

def index(f):
    """Read the index file `f` and return a dictionary that maps from
    id to an Index object, and also for GHCN-M v2 files that
    have 11-digit station identifiers and 12-digit record
    identifiers, maps id11 to a list of id12s.
    """

    d = {}
    for line in f:
        i = Index(line)
        d[i.id] = i
        if len(i.id) > 11:
            # Append the id12 to a list associated with each id11,
            # creating a list if necessary.
            id11 = i.id[:11]
            v = d.setdefault(id11, [])
            v.append(i.id)
    return d

class Index:
    """A single entry from the index file."""

    def __init__(self, line):
        self.id, self.year, self.whence = line.split()
        # Concatentate these two to form the 15 or 16 character
        # string that must match the GHCN-M file at the specified place.
        self.match = self.id + self.year
        self.whence = int(self.whence)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    name = "input/ghcnm.tavg.qca.dat"
    return build(open(name), open(name + '.index', 'w'))

if __name__ == '__main__':
    main()
