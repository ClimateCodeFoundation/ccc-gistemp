#!/usr/bin/env python
#
# $Id: //info.ravenbrook.com/project/ccc/master/tool/fortrantotext.py#1 $
#
# Command line tool for converting Fortran binary data into a generic
# text format.
#
# David Jones.  Ravenbrook Limited.  2008-08-13
#
# Written according to the Python 2.3 documentation:
# http://www.python.org/doc/2.3.5/
# It should therefore work with Python version 2.3 and any future
# compatible versions.

# Ravenbrook
import fort
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.3.5/lib/module-sys.html
import sys

def totext(file, output=sys.stdout):
    """The file argument should be a binary file opened for reading.  It
    is treated as a Fortran binary file and converted to a text format,
    emitted on the file object output.  Each (binary) record is treated
    as a sequence of words (words being "standard-sized ints" in the native
    byte-ordering), followed by a possible remainder sequence of bytes
    (where the record length is not a multiple of a word).  The output
    format is one line per record, with each word being output as a
    fixed width hexadecimal number, and each trailing byte being output
    as a 2-digit hexadecimal number.  Spaces separate.  A "standard-size
    int" is interpreted the same way that struct.unpack('=I', x)
    interprets it.
    """

    # Compute the width of a standard word according to Python's struct
    # module...
    w = len(struct.pack('=I', 0))
    # and a suitable string format.
    # http://www.python.org/doc/2.3.5/lib/typesseq-strings.html
    # The string format is of the form '%08x' but the value of 8 may be
    # replaced.
    fmt = '%%0%dx' % (2*w)

    f = fort.File(file)
    # Iterate over all the records
    for r in f:
        # We unpack as much as the record as we can as a sequence of
        # binary words (typically 32-bits each); then the rest as a
        # sequence of bytes.
        # Number of words
        n = len(r) // w
        sep = ''
        for i in struct.unpack('%dI' % n, r[:n*w]):
            output.write(sep + (fmt % i))
            sep = ' '
        # Remainder of record, as bytes
        for c in r[n*w:]:
            output.write(sep + ('%02x' % ord(c)))
            sep = ' '
        output.write('\n')

def main():
    if len(sys.argv[1:]) == 0:
        totext(sys.stdin)
    else:
        for n in sys.argv[1:]:
            totext(open(n, 'rb'))

if __name__ == '__main__':
    main()
