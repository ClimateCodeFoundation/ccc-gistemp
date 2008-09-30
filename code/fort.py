#!/usr/bin/env python
# $Id: //info.ravenbrook.com/project/ccc/master/code/fort.py#3 $

"""Handle (binary) Fortran files in Python.  A binary Fortran file can
be opened using the open method of this module; a File object is
returned that supports a writeline method (for writing records), a readline
method (for reading records) and the iterator protocol (for reading)."""

import struct

class File :
  """A Fortran file object."""
  
  # Underlying file descriptor.  Usually returned by builtin open, but
  # any object with a suitable read-like interface should do.
  fd = None

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

  def readline(self) :
    """For Fortran binary files this will return a binary record, as an
    arbitrary length binary string.  This is extremely unlikely to be
    terminated with a \\n.  On encountering end-of-file None is
    returned; this goes against the usual convention of returning '',
    but that could be returned as a valid 0-length record."""

    # This routine assumes each Fortran record is preceded and followed
    # by its length stored as an integer number of bytes.  This seems to
    # be fairly common.

    # Get record length in bytes ...
    l = self.readi()
    if l is None :
      return None
    # then the l bytes of the record itself ...
    r = self.fd.read(l)
    # then check the terminating word is also equal to l
    check = self.readi()

    # :todo: Should really raise a (documented) exception here.
    # Instead...
    assert check == l
    return r

  def writeline(self, record) :
    """Writes the record onto the Fortran binary file."""

    # As per readline we assume a record is preceded and followed by its
    # length.

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

  assert mode == 'rb'

  return File(file(name, mode))
