#!/usr/bin/env python
# $URL$
# $Rev$
# 
# subbox.py
# Copyright (C) 2009 Ravenbrook Limited.  See end of file for license.
#
# For reading trimmed big-endian GISTEMP sub-box data files.
#
# Nick Barnes, Ravenbrook Limited

import struct
# Ravenbrook
import fort

w = len(struct.pack('=I', 0))
wf = len(struct.pack('f', 0.0))
header_length = 8*w

class Line(tuple):
    """Each line in a trimmed sub-box file is converted into a Line
    instance.  It's a subclass of tuple; the elements of the tuple are
    the data items.  The per-line metadata (latitude, longitude, and
    so on) become attributes of the object.
    """

    def __new__(cls, record, **kwargs):
        items = (len(record)-header_length)//wf
        return tuple.__new__(cls,
            struct.unpack('%s%df' % (kwargs['bos'], items),
                          record[header_length:]))
        
    def __init__(self, record, bos='>'):
        keys = 'n lat_S lat_N lon_W lon_E stations station_months d'.split()
        self.bos = bos
        self.__dict__.update(zip(keys,
                                 struct.unpack('%s7if' % bos,
                                               record[:header_length])))
        for k in ('lat_S','lat_N','lon_W','lon_E'):
            self.__dict__[k] /= 100.0

    def __repr__(self):
        return ('<subbox.Line (%+06.2f,%+06.2f) (%+07.2f,%+07.2f): %d>' %
                (self.lat_S, self.lat_N, self.lon_W, self.lon_E, self.n))

class File:
    """Represents a GISTEMP subbox file. The file's metadata appears as
    instance properties.  The sequential file data is available by using
    the instance as an iterator.
    """

    def __init__(self, file, bos='>'):
        """Takes an open file object."""

        self.fortfile = fort.File(file, bos=bos)
        self.bos = bos
        headerline = self.fortfile.readline()
        keys = '''mo1 kq mavg monm monm4 iyrbeg
                  missing_flag precipitation_flag title'''.split()
        # See SBBXotoBX.f line 79, 77, 113
        self.__dict__.update(zip(keys, struct.unpack('>8i80s', headerline)))

    def __repr__(self):
        return '<subbox.File %r' % self.__dict__

    def __iter__(self):
        return self.next()

    def next(self):
        for l in self.fortfile:
            yield Line(l, bos=self.bos)

# This file is copyright (C) 2009 Ravenbrook Limited.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1.  Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
# 2.  Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
