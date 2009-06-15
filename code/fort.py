#!/usr/bin/env python
# $Id: //info.ravenbrook.com/project/ccc/master/code/fort.py#3 $

"""Handle (binary) Fortran files in Python.  A binary Fortran file can
be opened using the open method of this module; a File object is
returned that supports a writeline method (for writing records), a readline
method (for reading records) and the iterator protocol (for reading)."""

import struct

class Error(Exception):
    """An Exception."""
    pass

class FormatError(Error):
    """A problem with the expected format of a file."""
    pass

class File :
  """A Fortran file object."""
  
  #: Underlying file descriptor.  Usually returned by builtin open, but
  #: any object with a suitable read-like interface should do.
  fd = None

  #: The number of bytes in the natural word size used in the file.
  w = 4

  #: Byte Order and Size.  Specifies the byte order and sizes used in the
  #: underlying file.  This is used in a Python struct format so we use
  #: the same convention.  See
  #: http://docs.python.org/lib/module-struct.html
  #: '<' little-endian, '>' big-endian.
  bos = '<'

  def readi(self) :
    """Read a single signed integer from the file and return it.  If at
    end-of-file then None is returned.  Insufficient data (for example,
    only 3 bytes remaining but 4 required) will raise some sort of
    exception."""

    l = self.fd.read(self.w)
    if l == '' :
      return None
    l = struct.unpack(self.bos + 'i', l)[0]
    return l

  def __init__(self, fd, bos='@') :
    """A new Fortran File object that performs Fortran-style binary IO
    on the file object fd (normally this will be a file opened with the
    builtin 'open' in binary mode, but it doesn't have to be).  bos
    specifies the byte-order and size used in the record format of the
    Fortran binary-file.  It should be one of the characters used by the
    struct module to specify the byte-order and size (recap: '@' native,
    '<' little-endian, '>' big-endian).
    """

    self.fd = fd
    # We perform a test struct.pack.  This achieves two things: 1) sanity
    # checks the bos character; and, 2) tells us how many bytes we need
    # to read for the record length (which is stored in self.w).  We only
    # care about the result's length.
    try :
      self.w = len(struct.pack(bos + 'i', 7))
    except struct.error:
      raise 'suspect bos character: ' + bos
    # Byte Order and Size.  Specifies the byte order and sizes used in the
    # underlying file.  This is used in a Python struct format so we use
    # the same convention.  See
    # http://docs.python.org/lib/module-struct.html
    self.bos = bos

  def seek(self, offset, whence=0):
      """Proxy for the underlying ``fd`` seek.
      
      See ``file.seek`` in the standard python library for details.

      """
      return self.fd.seek(offset, whence)

  def close(self):
      """Close the underlying file.
      
      """
      return self.fd.close()

  def readline(self) :
    """For Fortran binary files this will return a binary record, as an
    arbitrary length binary string.  This is extremely unlikely to be
    terminated with a \\n.  On encountering end-of-file None is
    returned; this goes against the usual convention of returning '',
    but that could be returned as a valid 0-length record."""

    # This routine assumes each Fortran record is preceded and followed
    # by its length stored as an integer number of bytes.  This seems to
    # be fairly common.

    # Sadly, it's very difficult to make completely generic with respect
    # to self.w.  So, for now, we don't.
    assert self.w == 4

    # Location of start of record.
    at = 'unknown'
    try:
        at = self.fd.tell()
    except:
        pass
    # Get record length in bytes ...
    l = self.readi()
    if l is None :
      return None
    # then the l bytes of the record itself ...
    r = self.fd.read(l)
    # then check the terminating word is also equal to l
    check = self.readi()

    if check != l:
        at = str(at)
        raise FormatError(
          "Record prefix %d does not match suffix %r;"
          " record starting at %s." % (l, check, at))
    # :todo: Should really raise a (documented) exception here.
    # Instead...
    assert check == l
    return r

  def writeline(self, record) :
    """Writes the record onto the Fortran binary file."""

    # As per readline we assume a record is preceded and followed by its
    # length.

    assert self.w == 4
    l = struct.pack(self.bos + 'I', len(record))
    self.fd.write(l)
    self.fd.write(record)
    self.fd.write(l)

  def __iter__(self) :
    """Iterate over all records in file."""
    return self

  def next(self) :
    """Implements iterator protocol.  Avoid public use."""
    r = self.readline()
    if r is None :
      raise StopIteration
    return r

def open(name, mode='rb') :
  """Open the binary Fortran file called name.  mode is a mode string as
  per the builtin open function; for this version of this module it must
  be 'rb'."""

  assert 'b' in mode

  return File(file(name, mode))


def unpackRecord(line, start, fmt):
    """A simple function to unpack Fortran formatted records.

    This is only intended to be used during the conversion from Fortran
    Python code.

    The intention is to make is translation of Fortran such as::

      read(line(2:64), '(i4,i5,a12,i4,a36)') lat, lon, sid, ht, name
      read(line, '(i4,3i4,i4)') a, arr, b

    Into equivalent Python. In the above case you can write::

      lat, lon, sid, ht, name = unpackRecord(line, 2, 'i4,i5,a12,i4,a36')
      a, arr, b = unpackRecord(line, 1, 'i4,3i4,i4')

    This only supports a small subset of the Fortran format strings. Basically
    you can use:

    aN, iN
        A string or integer in a field width of ``N``. For example 'i5, s20'.
    MiN
        An array of ``M`` intgers, each in an ``N`` character field.
        The returned value is a list of integers.

    Note: This is throw-away code, which does not error checking.

    :Param line:
        The line of text to be 'unpacked'.
    :Param start:
        The start character in `line` to read from. Fortran indexing is used;
        i.e. the first character is as position 1.
    :Param format:
        The format string. This follows the Fortran conventions, but without the
        parentheses.
    """
    parts = []
    a = start - 1
    for f in fmt.replace(" ", "").split(","):
        count = 0
        if f.find('i') > 0:
            count, b = [int(x) for x in f.split('i')]
            typeCode = 'i'
        else:
            typeCode, b = f[0], int(f[1:])
        if count:
            arr = []
            for i in range(count):
                v = line[a:a + b]
                arr.append(int(v))
                a = a + b
            v = arr
        else:
            v = line[a:a + b]
            if typeCode == "i":
                v = int(v)
            a = a + b
        parts.append(v)
    return parts


def formatFloat(v):
    """Return a string representation of a float.

    The returned value is intended to match the default formatting that fortran
    code uses. This can be usefule during debugging, where it is necessary to
    compare printed values of floating point numbers.

    Notes:

    1. This has only been developed to be good enough for previous debugging
       tasks. It may not work as well on your data.
    2. Given the differences in precision between Fortran's and Python's
       floating point value, differences in the last digit are common.

    :Return:
        A string, within a 14 character field width.
    :Param v:
        The value to format. Should be a float.
    """
    s = "%.1f" % v
    if s.startswith("-"):
        s = s[1:]
    a, b = s.split(".")
    w = 7 - len(a)
    s = "%.*f" % (w, v) + "    "
    return s.rjust(14)
