"""Some common code for the CCC binary (unformatted) files.

Many (all?) of the CCC binary files share common features, such as the initial
header. This module provides support for these common parts.
"""
__docformat__ = "restructuredtext"

import struct
import copy

import fort

class CCBin(object):
    """Base class for binary stile file records.
    
    This simply abstracts out the common code from `CCHeader` and `CCRecord`.
    All clases derived from this class need to provide a member/property called
    ``_format``, which is used to pack/unpack the binary records using ithe
    ``struct`` module. They also need to over-ride the `getElements`,
    `setFromElements` and `initEmpty` methods.

    """
    def __init__(self, data=None, f=None):
        if data is not None:
            self.setFromBinary(data)
        elif f is not None:
            # TODO: Should provide proper error set
            raise Exception("Construction from file not (yet) supported")
        else:
            self.initEmpty()

    def copy(self):
        """Return a copy of this instance"""
        return copy.deepcopy(self)

    def setFromBinary(self, data):
        elements = struct.unpack(self._format, data[:self.size])
        self.setFromElements(elements)

    @property
    def binary(self):
        try:
            return struct.pack(self._format, *self.getElements())
        except struct.error, exc:
            raise ValueError("%s\n  Provided %d elements" % (
                exc, len(self.getElements())))


    @property
    def size(self):
        return struct.calcsize(self._format)

    #{Interface methods
    def initEmpty(self):
        """Initialise the record as empty.

        This is used when a record is created without any source data.
        It should set all members to suitable default/empty values.

        """
        raise NotImplementedError("Needs implementing in %s"
                % self.__class__.__name__)

    def getElements(self):
        """Return list of data elements.

        The returned list is used to write binary records. This list must be
        suitable for passing to the ``struct.pack`` function.

        """
        raise NotImplementedError("Needs implementing in %s"
                % self.__class__.__name__)

    def setFromElements(self, elements):
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

    """
    _format = "9i80s"
    def setFromElements(self, elements):
        self.info, self.title = list(elements[:9]), elements[9]

    def initEmpty(self):
        self.info = [0] * 9
        self.title = ""
        
    def getElements(self):
        return self.info + [self.title]

    @property
    def m1(self):
        return self.info[0]

    @property
    def m2(self):
        return self.info[8]

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

    Some of the CCC code generates file where each record contins information
    about the following record. This naturally means that, effectively,
    two records need to be loaded at once:

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

    Although this uses the `fort.File` class, it expects to buffer and write
    `CCBin` instances, not strings.

    - The `flush` method is invoked.
    - The `RecordFile` is destroyed.

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
        self._pos += self.w * 2 + record.size
        self.flush(low=1, high=10)

    def lastRecord(self):
        """Return the last record written, but not flushed.

        :Return:
            The las record written if any. If no record has been written since
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
