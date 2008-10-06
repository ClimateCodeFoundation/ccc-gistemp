#!/usr/bin/env python
# $Id$
# David Jones.  Copyright (C) 2008 Ravenbrook Limited.

"""Compatible implementation of uncompress(1) from traditional Unix
systems.
"""

import sys

# The reference for this implementation is
# http://www.opensource.apple.com/darwinsource/10.5/file_cmds-184/compress/zopen.c

# It seems to be difficult to find a description of the file format used
# by compress/uncompress.  Here is one reconstructed from the zcode.c source
# code:
#
# A compress'd file is a sequence of octets.  It comprises a header,
# followed by a body.  The header is 3 octets: 0x1F 0x9D BB. BB is the
# only variable part of the header; it specifies the maximum width of
# codes used and whether block compression (enabling the CLEAR code) is
# used.
#
# BB & 0x1f gives the maximum code word width in bits (which must be
# between 9 and 16, see notes below).  BB & 0x80 is 1 if block compression
# is enabled which means that the CLEAR code is used.
#
# The body consists of a sequence of codes packed into octets,
# essentially by filling up each octet from the right (least
# significant) end first.  Codes are packed into a bit sequence, the
# least significant bit of the code becomes the next bit of the
# sequence, then the next to least signiificant bit, and so on.
# The bit sequence is split into octets: bit i of the
# bit sequence corresponds to bit (i%8) of byte (i//8), where bit j of a
# byte is the bit corresponding to 2**j.
#
# codes start at 9 bits wide; subsequently the width of the codes in the
# body is always exactly enough to represent the largest possible code
# currently, up to a maximum width of 16.  If block compression is enabled
# (BB&0x80) then the CLEAR code, 257, resets the tables and therefore resets
# the current code width to 9.
#
# Notes of maximum code width
#
# Compressors may choose to restrict the maximum code width, it is
# possible to do this with the -b option to compress(1) for example.
# Some Unix systems are only capable of decompressing files with a
# maximum code width of 14 or less, so this is the portable maximum; in
# practice almost all reasonable Unix systems used a default maximum
# code width of 16, so this is what is typically found.
#
# This uncompressor copes with a maximum code width of up to 16 bits.

# zopen.c line 79
BITS = 16
# zopen.c line 90
BIT_MASK = 0x1f
# zopen.c line 91
BLOCK_MASK = 0x80
# zopen.c line 97
INIT_BITS = 9

# zopen.c line 201
FIRST = 257
# zopen.c line 202
CLEAR = 256

# zopen.c line 87
magic_header = '\x1f\x9d'
sizeof_header = 3

# zcode.c line 368
rmask = map(lambda s: ~(~0<<s), range(9))

def MAXCODE(n):
    """Maximum code value for a code comprised of n bits."""
    # zopen.c line 99
    return (1 << n) - 1

class Zfile:
    def __init__(self, name=None, fd=None):
        if name is None and fd is None:
            assert 'either name or fd should be supplied'
        if fd is None:
            fd = open(name, 'rb')
        self.f = fd

        # Mostly from zopen, zopen.c line 686
        self.maxmaxcode = 1 << BITS
        # Is 0/1 in C, but only used as a boolean
        self.clear_flg = False

        self.header = self.f.read(sizeof_header)
        # :todo: make exception
        assert len(self.header) == sizeof_header
        assert self.header.startswith(magic_header)
        self.maxbits = ord(self.header[2])
        self.block_compress = self.maxbits & BLOCK_MASK
        self.maxbits &= BIT_MASK
        if self.maxbits > BITS:
            raise 'Compressed codes have too many bits'


    # zcode.c line 560
    def getcode(self):
        """Yield input as an iterable sequence of codes."""

        n_bits = INIT_BITS
        maxcode = MAXCODE(n_bits)
        roffset = 0
        size = 0

        while True:
            if (self.clear_flg or
                    roffset >= size or
                    len(self.prefixof) > maxcode):
                if len(self.prefixof) > maxcode:
                    n_bits += 1
                    if n_bits == self.maxbits:
                        maxcode = self.maxmaxcode
                    else:
                        maxcode = MAXCODE(n_bits)
                if self.clear_flg:
                    n_bits = INIT_BITS
                    maxcode = MAXCODE(n_bits)
                    self.clear_flg = False
                gbuf = map(ord, self.f.read(n_bits))
                if not gbuf:
                    return
                roffset = 0
                # zcode.c line 597
                # Cute optimisation.  size is used to detect when the
                # buffer needs refilling (see roffset >= size, above).
                # It should be no more than the number of bits
                # occupied by entire codes and more than the number of bits
                # occupied by all the codes except the last.
                size = (len(gbuf) << 3) - (n_bits - 1)
            r_off = roffset
            bits = n_bits

            # In this Python version bp is a string index, from 0.
            bp = r_off >> 3
            r_off &= 7

            # print gbuf, roffset, bits, bp
            gcode = gbuf[bp] >> r_off ; bp += 1
            bits -= 8 - r_off
            r_off = 8 - r_off       # Now, roffset into gcode word

            # Middle 8-bit byte, if any
            if bits >= 8:
                gcode |= gbuf[bp] << r_off ; bp += 1
                r_off += 8
                bits -= 8
            # High order bits
            if bits:
                gcode |= (gbuf[bp] & rmask[bits]) << r_off
            roffset += n_bits

            yield gcode
        
    def read1(self):
        """Yield the characters of the uncompressed stream one at a
        time.  Returns an iterator.
        """

        self.prefixof = [0] * 256
        self.suffixof = range(256)
        if self.block_compress:
            self.prefixof.append(None)
            self.suffixof.append(None)

        # zcode.c line 505
        codes = self.getcode()
        finchar = oldcode = codes.next()
        yield chr(finchar)

        # In the C code, zcode.c, the variable free_ent records the
        # index of the next free entry in the coding table.  In this
        # Python version we dispense with this variable.  The coding
        # table is stretchy, it is a list; the len(self.prefixof)
        # replaces use of free_ent.  prefixof and suffixof should be the
        # same size throughout; they are each extended when a new code
        # is used.

        # A list of characters (length 1 strings)
        de_stack = []
        # This loop is exited when codes, which is an iterator, is
        # exhausted and the call to its next method (see codes.next(),
        # below) raises the StopIteration exception.  Sneaky.
        while True:
            code = codes.next()
            if code == CLEAR and self.block_compress:
                self.prefixof = [0]*256
                self.suffixof = self.suffixof[:256]
                self.clear_flg = True
                code = codes.next()
            incode = code

            # Special case for KwKwK string
            if code >= len(self.prefixof):
                de_stack += chr(finchar)
                code = oldcode

            # Generate output characters in reverse order
            while code >= 256:
                de_stack += chr(self.suffixof[code])
                code = self.prefixof[code]
            finchar = self.suffixof[code]
            de_stack += chr(finchar)

            while de_stack:
                yield de_stack.pop()

            # Generate the new entry
            code = len(self.prefixof)
            if code < self.maxmaxcode:
                self.prefixof.append(oldcode)
                self.suffixof.append(finchar)
            oldcode = incode


z = Zfile(fd=sys.stdin)
for i,c in enumerate(z.read1()):
    sys.stdout.write(c)

# TODO
# None for EOF in getcode

