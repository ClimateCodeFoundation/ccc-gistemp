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

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    return build(open('input/v2.mean'), open('work/v2.mean.index', 'w'))

if __name__ == '__main__':
    main()
