#!/usr/bin/env python
# $URL$
# $Rev$
#
# ncartotext.py
#
# David Jones.  Ravenbrook Limited.  2008-08-13

"""
Command line tool for converting Fortran binary NCAR file
into text format.  The NCAR format is the input to STEP3 of GISTEMP
(hence the ouput of STEP2).  It is at least partly documented in
the header comment of the file to.SBBXgrid.f .

The name NCAR is taken from the comments at the top of to.SBBXgrid.f;
I have no idea if this is an actual documented format, but NCAR
probably stands for National Center for Atmospheric Research
(in the US).

The outpu is the station data in GHCN v2.mean format, unless the -m
option used, in which case the metadata is ouput and the data is not.

Metadata are output as one line per station:

ID 101603550002 SKIKDA                         *UC +36.90+007.00+0007 M1043

Considering the line as a series of words separated by spaces, the
words in order are:
ID: the literal word "ID" to remind you that the station ID follows.
record-id: 12-digit record identifier.  The first 3 digits identify the
  country; the first 11 digits (including the 3 for country) are the
  station identifier.  A station can have more than one ("duplicate")
  record, the 12th digit is used for different records for the same
  station.
name: The name of the station.  Any spaces internal to the station
  name are replaced with underscores.  Trailing spaces are preserved,
  making it neater, but also preserving both the fixed-width format
  and the "words separated by spaces" format.
misc: A 3-character word giving miscellaneous metadata.  The first
  character gives the USHCN-brightness index (0/1/2), '*' when this
  metadata is absent (typicaly for non-USHCN stations); second
  character gives the population index (U/S/R for Urban, Suburban,
  Rural, respectively); third character gives the GHCN-brightness
  index (A/B/C).
location: The location of the station in ISO 6709 notation.  Latitude
  and longitude (in decimal degrees), height.  All values signed.
  When height is not available it appears as -0999
first-month: the character "M" followed by a 4-digit number.  The
  number denotes the first month of the station data, where 1 denotes
  January of the first year (see YRBEG from the file metadata);
  subsequent months have subsequent numbers.  Thus 1043 denotes
  November of YRBEG+86.

Written according to the Python 2.3 documentation:
http://www.python.org/doc/2.3.5/
It should therefore work with Python version 2.3 and any future
compatible versions.
"""

# http://www.python.org/doc/2.3.5/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os
# http://www.python.org/doc/2.3.5/lib/module-re.html
import re
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.3.5/lib/module-sys.html
import sys

# Allows imports from the code directory.
sys.path.append(os.path.join(os.getcwd(),'code'))
# Clear Climate Code
import fort

# Move somewhere else?
def iso6709(lat, lon, height=None):
    """Convert geographic co-ordinates into ISO 6709 format."""

    assert -90 <= lat <= 90
    assert -180 <= lon <= 180
    s = '%+06.2f%+07.2f' % (lat, lon)
    if height is not None:
        s += '%+05d' % height
    return s

def totext(file, output=sys.stdout, error=sys.stderr, metaflag=False):
    """The file argument should be a binary file opened for reading.  It
    is treated as a Fortran binary file and converted to a text format,
    emitted on the file object output.
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

    r = f.readline()
    # First record contains various INFO items.
    # When referring to comments in to.SBBXgrid.f recall that Fortran
    # arrays are typically indexed from 1 onwards; Python from 0
    # onwards.  Therefore INFO(2) corresponds to a[0]
    a = struct.unpack('9i', r[:9*w])
    mfirst = a[0]
    mlast = a[8]
    kq = a[1]
    mavg = a[2]
    monm = a[3]
    recsize = a[4]
    yrbeg = a[5]
    bad = a[6]
    trace = a[7]

    if metaflag:
        output.write('KQ=%d MAVG=%d MONM=%d YRBEG=%d BAD=%d TRACE=%d\n' %
            (kq, mavg, monm, yrbeg, bad, trace))

    # Length of record trail, the non-variable part, in bytes.
    ltrail = 15*w
    # Iterate over remaining records
    for r in f:
        trail = r[-ltrail:]
        lat,lon,id,height = struct.unpack('4i', trail[:4*w])
        lat *= 0.1
        lon *= 0.1
        if len(r[:-ltrail]) != w*(mlast-mfirst+1):
            error.write(('Station ID %09d has suspect record length.' +
                'mfirst=%s mlast=%d record-length=%d\n') %
                (id, mfirst, mlast, len(r)))
        name = trail[4*w:-2*w]
        # Some metadata is stored in the name field, *sigh*
        meta = name[-6:]
        name = name[:-6]
        # 3 digit country code
        cc = meta[3:6]
        meta = meta[0:3]
        # Prepend country code to station ID.  Note: id becomes a string.
        id = '%s%09d' % (cc, id)
        # Replace any spaces in meta with '*', mostly to preserve the
        # "words separated by spaces" format.  Note that all these
        # replacements are where (non US) stations do not have
        # USHCN-brightness indexes.
        meta = meta.replace(' ', '*')
        # It just so happens that underscore does not appear in any of
        # the "name" fields, so we use that instead of space.  That
        # preserves the "words separated by spaces" format.
        name = name.replace(' ', '_')
        # It's tidier if we return the trailing underscores back into
        # spaces.
        m = re.search('_*$', name)
        name = name[:m.start()] + ' '*(m.end() - m.start())
        if metaflag:
            output.write('ID %s %s %s %s M%04d\n' %
                (id, name, meta, iso6709(lat, lon, height), mfirst))
        else:
            n = len(r[:-ltrail])//w
            data = struct.unpack('%di' % n, r[:-ltrail])
            # Convert to 0 == January indexing:
            mfirst -= 1
            assert mfirst >= 0
            year = mfirst // 12
            year += yrbeg
            # Pad data at beginning...
            m = mfirst % 12
            data = (bad,)*m + data
            # ... and at end.
            m = (-len(data)) % 12
            data = data + (bad,)*m
            assert 0 == len(data) % 12
            def changebad(datum):
                """Convert to GHCN v2.mean BAD format."""
                if datum == bad:
                    return -9999
                return datum
            data = map(changebad, data)
            for y,yeardata in enumerate(data[i:i+12]
                                        for i in range(0,len(data),12)):
                output.write('%s%d' % (id, year + y))
                output.write(('%5d'*12 + '\n') % tuple(yeardata))
        mfirst, mlast = struct.unpack('2I', trail[-2*w:])

def main(argv=None):
    if argv is None:
        argv = sys.argv
    meta = False
    opt,arg = getopt.getopt(argv[1:], 'm')
    for o,v in opt:
        if o == '-m':
            meta = True
    if len(arg) == 0:
        totext(sys.stdin, meta=meta)
    else:
        for n in arg:
            totext(open(n, 'rb'), metaflag=meta)

if __name__ == '__main__':
    main()
