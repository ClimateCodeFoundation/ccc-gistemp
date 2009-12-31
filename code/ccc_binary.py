"""Some common code for the CCC binary (unformatted) files.

Many (all?) of the CCC binary files share common features, such as the initial
header. This module provides support for these common parts.
"""
__docformat__ = "restructuredtext"

import struct
import copy

import fort

def prop(method):
    """Property decorator that keeps documentation.

    This simply invokes the ``property`` built-in function, but it passes the
    method's docstring as well; making generated API documentation more useful.

    :Param method:
        The method to convert to a property.

    """
    return property(method, doc=method.__doc__)


def mprop(attrExpr):
    """Like `prop`, but supports simple modification.

    Simple modification is supported, i.e. assignment of the underlying
    attribute.

    :Param attrExpr:
        An partial expression that can be used to set the property,
        as in self.<attrExpr> = v.

    """
    x = 'def doSet(self, value): self.%s = value' % attrExpr
    d = {}
    exec x in locals(), d
    doSet = d["doSet"]

    def makeProp(method):
        return property(fget=method, fset=doSet, doc=method.__doc__)
    return makeProp


class CCBin(object):
    """Base class for binary style file records.
    
    This simply abstracts out the common code from `CCHeader` and `CCRecord`.
    All clases derived from this class need to provide a member/property called
    ``_format``, which is used to pack/unpack the binary records using ithe
    ``struct`` module. They also need to over-ride the `getElements`,
    `setFromElements` and `initEmpty` methods.

    """
    def __init__(self, data=None, f=None, bos="@"):
        """Constructor:

        :Param data:
            Binary data used to set the contents. See `setFromBinary` for more
            details.

        :Param f:
            Not currently used. It is intended to allow the contents to be set
            by reading from a `fort.File`.

        """
        self._bos = bos
        if data is not None:
            self.setFromBinary(data)
        elif f is not None: #pragma: unsupported
            # TODO: Should provide proper error set
            raise Exception("Construction from file not (yet) supported")
        else:
            self.initEmpty()

    def copy(self):
        """Copy this instance to make a new record.
        
        This performs a deep copy so you can modify the copy without affecting
        the original.

        :Return:
            A new copy of the record.

        """
        return copy.deepcopy(self)

    def setFromBinary(self, data):
        """Set the contents from a binary string.

        :Param data:
            The binary data to use to set the content. This *must* be a string
            and be at least as long as the string returned by the `binary`
            property.

        """
        elements = struct.unpack(self._bos + self._format, data[:self.size])
        self.setFromElements(elements)

    @property
    def binary(self):
        """The record as a binary string.

        This is as the record will be when written to a Fortran binary file, but
        excludes record header and trailer.

        """
        try:
            return struct.pack(self._bos + self._format, *self.getElements())
        except struct.error, exc: #pragma: debug
            raise ValueError("%s\n  Provided %d elements" % (
                exc, len(self.getElements())))


    @property
    def size(self):
        """The size of the record, when converted to a binary string.

        This is functionally equivalent to len(self.binary), but should be
        quicker.

        """
        return struct.calcsize(self._bos + self._format)

    #{Interface methods
    def initEmpty(self): #pragma: unsupported
        """Initialise the record as empty.

        This is used when a record is created without any source data.
        It should set all members to suitable default/empty values.

        """
        raise NotImplementedError("Needs implementing in %s"
                % self.__class__.__name__)

    def getElements(self): #pragma: unsupported
        """Return list of data elements.

        The returned list is used to write binary records. This list must be
        suitable for passing to the ``struct.pack`` function.

        """
        raise NotImplementedError("Needs implementing in %s"
                % self.__class__.__name__)

    def setFromElements(self, elements): #pragma: unsupported
        """Set the record's members from a list of elements.
  
        This is the inverse of `getElements`. The elements will typically have come
        from ``struct.unpack`` and should be used to set the record's members.

        """
        raise NotImplementedError("Needs implementing in %s"
                % self.__class__.__name__)
    #}


class CCHeader(CCBin):
    """Holds the data for a binary temperature file.

    Many of the Fortran programs read and write binary files, which have
    a common format:
    ::

        +-----------+
        | CCHeader  |
        +-----------+
        | CCRecord  |
        +- - - - - -+
        | CCRecord  |
        +- - - - - -+
        :           :
        :           :
        +-----------+
     
    This class make accessing the header in the Python world relatively clean
    and easy.

    The header itself has the following structure::

        +-----------+ - - - - - - - - +-------------------+
        | info[9]   |              0  |  m1               |
        +-----------+ -               +-------------------+
        | title[80] |   \          1  |  kq               |
        +-----------+    .            +-------------------+
                         |         2  |  mavg             |
                         .            +-------------------+
                         |         3  |  monm             |
                         .            +-------------------+
                         |         4  |  recsize          |
                         .            +-------------------+
                         |         5  |  yrbeg            |
                         .            +-------------------+
                         |         6  |  bad              |
                         .            +-------------------+
                         |         7  |  trace            |
                         .            +-------------------+
                          \        8  |  m2               |
                            - - - - - +-------------------+

    The `title` should always be left justified and padded with spaces.
    The `title` property enforces this.

    In the Fortran, the info structure indices start at 1.

    .. Note::
        The names of the elements of the ``info`` are based on names used
        within certain Fortran files. The Fortran code does not have 
        consitent naming conventions and it is *very* likeky better names will
        emerge as the CCC code evolves.

    """
    # The header format is 9 integers followed by an 80 character string.
    _format = "9i80s"
    def setFromElements(self, elements):
        self.info, self.title = list(elements[:9]), elements[9]

    def initEmpty(self):
        self.info = [0] * 9
        self.title = ""
        
    def getElements(self):
        return self.info + [self.title]

    def _set_title(self, s):
        self._title = s.ljust(80)

    def _get_title(self):
        return self._title.rstrip()

    title = property(_get_title, _set_title,
        doc="""The title of the data file.

        This is a plain string which is left justified and automaticlly padded
        wih spaces.""")

    @property
    def m1(self):
        """The month (or year) number of the first idata element in the
        following record."""
        return self.info[0]

    @property
    def m2(self):
        """The month (or year) number of the last idata element in the
        following record."""
        return self.info[8]

    @property
    def kq(self):
        """Don't know yet!"""
        return self.info[1]

    @property
    def mavg(self):
        """Don't know yet!"""
        return self.info[2]

    @property
    def monm(self):
        """Don't know yet!"""
        return self.info[3]

    @property
    def recsize(self):
        """Don't know yet!"""
        return self.info[4]

    @property
    def yrbeg(self):
        """Don't know yet!"""
        return self.info[5]

    @property
    def bad(self):
        """Don't know yet!"""
        return self.info[6]

    @property
    def trace(self):
        """Don't know yet!"""
        return self.info[7]


class SBBX_Header(CCBin):
    """The header for SBBX files.

    This is very similar to `CCHeader`., but the info structure only
    has 8 elements rather than 9.

    The info elements are (see code/STEP4_5/trimSBBX.f for original
    details)::

        +-------------------+
      0 |  m1               | [Note 1]
        +-------------------+
      1 |  kq               | quantity flag [Note 2]
        +-------------------+
      2 |  mavg             | time avg flags [Note 3]
        +-------------------+
      3 |  nm               | length of each time record
        +-------------------+
      4 |  recsize          | size of data record length
        +-------------------+
      5 |  yrbeg            | first year of each time record
        +-------------------+
      6 |  bad              | value used to indicate missing data
        +-------------------+
      7 |  trace            | flag for precipitation trace
        +-------------------+

    Note 1:
        This is set to ``1`` to indicate that the file is not
        trimmed/reorganized.

    Note 2:
        The quantity flag values are not defined in
        ``code/STEP4_5/trimSBBX.f``.

    Note 3:
        The time average flag values are:
*           1-4  => DJF - SON
            5    => ANN
            6    => MONTHLY
            7    => SEAS
            8-19 =>JAN-DEC

    """
    # The header format is 8 integers followed by an 80 character string.
    _format = "8i80s"
    def setFromElements(self, elements):
        self.info, self.title = list(elements[:8]), elements[8]

    def initEmpty(self):
        self.info = [0] * 8
        self.title = ""
        
    def getElements(self):
        return self.info + [self.title]

    def _set_title(self, s):
        self._title = s.ljust(80)

    def _get_title(self):
        return self._title.rstrip()

    title = property(_get_title, _set_title,
        doc="""The title of the data file.

        This is a plain string which is left justified and automaticlly padded
        wih spaces.""")

    @mprop("info[0]")
    def m1(self):
        """Provisional: The month (or year) number of the first idata element
        in the following record."""
        return self.info[0]

    @mprop("info[7]")
    def trace(self):
        """Flag for precipitation trace."""
        return self.info[7]

    @mprop("info[1]")
    def kq(self):
        """Quantity flag."""
        return self.info[1]

    @mprop("info[2]")
    def mavg(self):
        """Time average flags."""
        return self.info[2]

    @mprop("info[3]")
    def nm(self):
        """Length of each time record."""
        return self.info[3]

    @mprop("info[4]")
    def recsize(self):
        """Size of data record length."""
        return self.info[4]

    @mprop("info[5]")
    def yrbeg(self):
        """First year of each time record."""
        return self.info[5]

    @mprop("info[6]")
    def bad(self):
        """Value used to indicate missing data."""
        return self.info[6]


class CCRecord(CCBin):
    """Holds the data for a binary temperature record.

    Many of the Fortran programs read and write binary files, which have
    a common format for each record. This class is designed to make it easy to
    read/write those records in the Python world.

    """
    def __init__(self, count, **kwargs):
        self.count = count
        super(CCRecord, self).__init__(**kwargs)

    def initEmpty(self):
        self.idata = [0] * self.count
        self.Lat = 0
        self.Lon = 0
        self.ID = 0
        self.iht = 0
        self.name = ""
        self.m1 = 0
        self.m2 = 0
        
    @property
    def _format(self):
        return "%diiiii36sii" % self.count

    def getElements(self):
        return self.idata + [
            self.Lat, self.Lon, self.ID, self.iht, self.name,
            self.m1, self.m2]

    def setFromElements(self, elements):
        self.idata, trailer = (list(elements[:self.count]),
                elements[self.count:])
        self.Lat = trailer[0]
        self.Lon = trailer[1]
        self.ID = trailer[2]
        self.iht = trailer[3]
        self.name = trailer[4]
        self.m1 = trailer[5]
        self.m2 = trailer[6]


class BufferedOutputRecordFile(object):
    """Simple buffer for output records.

    Much of the CCC code generates file where each record contains information
    about the following record. This naturally means that, effectively, two
    records need to be loaded in memory at once:

    - The records that is almost ready to be written.
    - The next rcord, which is needed to fill in a few details in the other
      record.

    This makes it fiddly to write natural loops, as in:

        for record in f:
            make_output_record()
            outF.writeRecord()

    This class solves this by looking like a file, but allowing the last record
    to be retrieved and modified using the `lastRecord` method. In practice it
    does this by buffering records and never committing the most recent record
    to the underlying file until either:
 
    - The `flush` method is invoked with arguments low=0, high=0 (the defaults)
    - The `RecordFile` is destroyed.
 
    Although this uses the `fort.File` class, it expects to buffer and write
    `CCBin` instances, not strings.

    """
    def __init__(self, path):
        """Constructor:

        :Param path:
            The path of the output file. It will be opened with a mode of "wb".

        """
        self.f = fort.open(path, "wb")
        self._buffer = []
        self._pos = 0

    def writeRecord(self, record):
        """Write a single record.

        Important. This does not copy the record.
        Calling this several times with the same record instance will result in
        filling the buffer with identical instances.

        :Param record:
            A single record. This must be a `CCBin` instance or equivalent.
            Specifically ``record.size`` must provide the number of bytes
            in the record and ``rec.binary`` a string containing those bytes.

        """
        self._buffer.append(record)
        self._pos += self.f.w * 2 + record.size
        self.flush(low=1, high=10)

    def lastRecord(self):
        """Return the last record written, but not flushed.

        :Return:
            The last record written if any. If no record has been written since
            the last time all records were flushed then it returns ``None``.

        """
        if self._buffer:
            return self._buffer[-1]

    def tell(self):
        """Like standard ``file`` tell method.

        :Return:
            The effective file's current position. The same value would be
            obtains by doing::

                self.flush()
                return self.f.fd.tell()

            But, of course, no flushing occurs.
        """
        return self._pos

    def __getattr__(self, name):
        """Provide proxy access to the underlying `fort.File` instance."""
        return getattr(self.f, name)

    def flush(self, low=0, high=0):
        """Force buffered records to be written to file.
        
        :Param low, high:
            Low and high water marks. If ``high`` records or fewer are buffered
            then nothing is flushed. If records are flushed ``low`` records will
            be retained in the buffers.

            These are mainly intended for internal use.
        """
        if len(self._buffer) <= high:
            return
        while len(self._buffer) > low:
            rec = self._buffer.pop(0)
            self.f.writeline(rec.binary)

    def __del__(self):
        """Ensure records get flushed at deletion."""
        self.flush()
        self.f.close()
