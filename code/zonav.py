#!/usr/bin/env python
# $URL$
# $Rev$
#
# zonav.py
#
# David Jones, Ravenbrook Limited, 2009-12-11
#

"""
Perform Zonal Averaging.

(this is typically part of STEP 5)

The input is a box file, usually work/BX.Ts.GHCN.CL.PA.1200 .  The data
in the boxes are combined to produce averages over various latitudinal
zones.  The output is a zone file,
usually work/ZON.Ts.ho2.GHCN.CL.PA.1200.step1 (the output can be
converted to text with code/zontotext.py).

14 belts are produced.  The first 8 are the basic belts that are used
for the equal area grid, the remaining 6 are combinations:

(The Fortran uses 1 to 14, but here in Python we use 0 to 13)

  0 64N - 90N
  1 44N - 64N (asin 0.9)
  2 24N - 44N (asin 0.7)
  3 Equ - 24N (asin 0.4)
  4 24S - Equ
  5 44S - 24S
  6 64S - 44S
  7 90S - 64S
  8 24N - 90N (0 + 1 + 2)
  9 24S - 24N (3 + 4)
 10 90S - 24S (5 + 6 + 7)
 11 northern hemisphere (0 + 1 + 2 + 3)
 12 southern hemisphere (3 + 4 + 5 + 6)
 13 global (all belts 0 to 7)
"""

def zonav(inp, out, log):
    """Take an open file *inp* of boxed data and produce zonal means on
    *out*.  *log* is used for logging (mostly textual descriptions).
    """

    # local module
    import fort
    # http://www.python.org/doc/2.4.4/lib/module-struct.html
    import struct

    inp = fort.File(inp, bos='>')
    out = fort.File(out, bos='>')

    # Number of boxes (regions).  zonav.f line 71
    NRM = 80
    # Number of basic belts.  zonav.f line 71
    JBM = 8
    # Maximum number of boxes in a basic belt.  zonav.f line 71
    IBMM = 16
    # Number of special zones.  zonav.f line 71
    NZS = 3
    # Number of required common data in order to combine.
    NOVRLP = 20

    r = inp.readline()
    info = struct.unpack('>8i80s', r)
    titlei = info[8]
    info = info[:8]
    kq = info[1]
    km = 1
    if info[2] == 6:
        km = 12
    if info[2] == 7:
        km = 4
    ml = info[3]
    nyrsin = info[3]/km
    # The rest of the code is structured so that MONM and IYRBEG could
    # be different from these values (for example, to make the output
    # file smaller, I think).
    iyrbeg = info[5]
    monm = ml
    # One more than the last year with data
    yearlimit = monm/km + iyrbeg
    # Index of first month to be output
    mfout = (iyrbeg - info[5])*km

    # Initial and last years of the reference period
    IYBASE = 1951
    LYBASE = 1980
    # *reference* is a pair of (base,limit) that specifies the reference
    # period (base period).  In the Fortran NFB and NFL are used.  Line 110
    reference = (IYBASE-info[5], LYBASE-info[5]+1)
    # Typically 9999
    XBAD = info[6]
    last = info[6]
    trace = info[7]

    # Write output file header
    outinfo = list(info)
    outinfo[3] = monm
    outinfo[5] = iyrbeg
    out.writeline(struct.pack('>8i', *outinfo) + titlei +
     ' zones:  90->64.2->44.4->23.6->0->-23.6-'
     '>-44.4->-64.2->-90                      ')

    # In the Fortran, zonav.f line 127, there is a call to GRIDEA (Grid
    # Equal Area - an anachronism), and the arrays are filled in via
    # global variables.
    # Here in Python, we just return them.
    ibm,kzone,titlez = grid()
    assert len(ibm) == JBM

    # zonav.f line 83
    # :init:lenz:
    lenz = [None]*JBM

    # :todo: move combine and tavg into proper module
    from step5 import combine, tavg

    # Length, in bytes, of a float as stored in the binary file.
    w = len(struct.pack('>f', 0.0))
    # Length, in bytes, of an integer as stored in the binary file.
    wi = len(struct.pack('>i', 0))
    # If these are different, then things are likely to break.
    assert w == wi
    # Equivalent to zonav.f line 107, sort of.
    wt = [None] * JBM
    avg = [None] * JBM
    # zonav.f line 131 to 177
    for jb in range(JBM):
        ar = [None] * ibm[jb]
        wtr = [None] * ibm[jb]
        lenr = [None] * ibm[jb]
        for ib in range(ibm[jb]):
            record = inp.readline()
            arpart = record[:w*ml]
            wtrpart = record[w*ml:w*ml*2]
            metapart = record[w*ml*2:]
            assert len(arpart) == len(wtrpart) == w*ml
            # Not entirely clear from zonav.f but the metadata that
            # trails each record is 5 words:
            # NG, LATS, LATN, LONW, LONE
            # (We only use NG, Number Good; the location information is
            # ignored.)
            assert len(metapart) == wi*5
            ar[ib] = struct.unpack('>%df' % ml, arpart)
            wtr[ib] = struct.unpack('>%df' % ml, wtrpart)
            lenr[ib], = struct.unpack('>I', metapart[:wi])
        # zonav.f lines 137 to 140
        lntot = sum(lenr)
        if lntot == 0:
            log.write('**** NO DATA FOR ZONE %d %s' % (jb, titlez[jb]))
            continue
        z = zip(lenr, range(ibm[jb]))
        from step3 import sort
        sort(z, lambda x,y: y[0]-x[0])
        lenr,IORD = zip(*z)
        nr = IORD[0]
        # Copy the longest region's record into *wt* and *avg*.
        # Using list both performs a copy and converts into a mutable
        # list.
        wt[jb] = list(wtr[nr])
        avg[jb] = list(ar[nr])
        # And combine the remaining series.
        bias = [0]*12
        wtm = [0]*12
        # zonav.f line 153
        for n in range(1,ibm[jb]):
            nr = IORD[n]
            if lenr[n] == 0:
                break
            # zonav.f does not declare the variable WTM, but it uses a
            # modified version of CMBINE that does not touch it.
            # Amazing.
            # In Python we pass genuine lists, but we do not use the
            # results.
            # :todo: could the id parameter (the last one) be something
            # more useful, like "Zone 4"? Is it only used for
            # diagnostics?
            combine(avg[jb], bias, wt[jb], ar[nr], 1,nyrsin,
              wtr[nr], wtm, km, jb)
        # zonav.f line 161
        bias = tavg(avg[jb], km, nyrsin, reference[0], reference[1], True,
          "Belt %d" % jb)
        lenz[jb] = 0
        m = 0
        # zonav.f line 164
        for iy in range(nyrsin):
            for k in range(km):
                m = iy*km + k
                if avg[jb][m] == XBAD:
                    continue
                avg[jb][m] -= bias[k]
                lenz[jb] += 1

        zoneout(out, log, avg[jb], wt[jb], titlez[jb])

    # *lenz* contains the lemgths of each zone 0 to 7 (the number of
    # valid months in each zone).
    # zonav.f line 181
    lenz, iord = sort_perm(lenz)
    # zonav.f line 182
    for jz in range(NZS+3):
        if lenz[0] == 0:
            # zonav.f carries on round the loop if there is no data,
            # seems pointless to me, since no data on the zone with
            # the longest record means no data anywhere.
            raise Error('**** NO DATA FOR ZONE %d %s' %
              (JBM+jz, titlez[JBM+jz]))
        # Find the longest basic belt that is in the special zone.
        # zonav.f line 188.
        for j1 in range(JBM):
            # kzone records which basic belts (0 to 7) are in each
            # of the "special zones".
            if iord[j1] in kzone[jz]:
                break
        else:
            # Slightly obscure Python feature; we get here when the
            # "for" loop exits normally.  Which in this case, is a
            # problem (we always expect to find at least one basic
            # belt for every special zone).
            raise Error('No basic belt in special zone %d.' % jz)
        jb = iord[j1]
        # zonav.f line 192.
        wtg = list(wt[jb])
        avgg = list(avg[jb])
        # zonav.f line 195
        # Add in the remaining latitude belts JB with JB in kzone[jz].
        for j in range(j1+1,JBM):
            jb = iord[j]
            if jb not in kzone[jz]:
                continue
            if lenz[j] == 0:
                # Not convinced this behavious is either correct or
                # worth preserving from the zonav.f code.
                break
            combine(avgg, bias, wtg, avg[jb], 1,nyrsin,
              wt[jb], wtm, km, "Zone %d" % (JBM+jz))
        else:
            # zonav.f line 202
            # Set BIAS=time average over the base period if IYBASE > 0
            if reference[0] > 0:
                bias = tavg(avgg, km, nyrsin, reference[0], reference[1],
                  True, "Zone %d" % (JBM+jz))
            else:
                # Not sure what the Fortran code zonav.f does when
                # NFB <= 0, so assert here.
                bias = None
                assert 0
            m = 0
            for iy in range(nyrsin):
                for k in range(km):
                    if avgg[m] != XBAD:
                        avgg[m] -= bias[k]
                    m += 1
        zoneout(out, log, avgg, wtg, titlez[JBM+jz])
    out.flush()

def zoneout(out, log, average, weight, title):
    """Output, onto the files *out* and *log*, the zone record in
    the arrays *average* and *weight*, labelled *title*."""

    # http://www.python.org/doc/2.4.4/lib/module-math.html
    import math
    # http://www.python.org/doc/2.4.4/lib/module-struct.html
    import struct

    def f(x):
        """Format a monthly anomaly as a string, for logging."""
        return "%4d" % math.floor((10*x) + 0.5)

    def filerepr(a):
        """Take an array (list,sequence,tuple) of floats and return
        a string representing the binary file representation.  Which
        is simply each element as a 4-byte float.
        """

        return struct.pack('>%df' % len(a), *a)

    def swaw(s):
        """Swap the words of a string, so that every 4 bytes are
        reversed.  This emulates the way that Fortran strings get
        written out using the SWRITE subroutine in zonav.f and
        friends.
        """

        l = len(s)//4
        return struct.pack('<%dI' % l, *struct.unpack('>%dI' % l, s))

    print >> log, 'zonal mean', title
    spaces = [' '*5]
    for l in [map(f, average[i:i+12])+spaces+map(f, average[i+12:i+24])
                for i in range(0, len(average), 24)]:
        print >> log, ''.join(l)
    print >> log, 'weights'
    for l in [map(f, weight[i:i+12])+spaces+map(f, weight[i+12:i+24])
                for i in range(0, len(weight), 24)]:
        print >> log, ''.join(l)

    # zonav.f line 175 and line 215
    out.writeline(filerepr(average) + filerepr(weight) + swaw(title))

def sort_perm(a):
    """The array *a* is sorted into descending order.  The fresh sorted
    array and the permutation array are returned as a pair (*sorted*,
    *iord*).  The original *a* is not mutated.

    The *iord* array is such that `a[iord[x]] == sorted[x]`.
    """
    z = zip(a, range(len(a)))
    from step3 import sort
    sort(z, lambda x,y: y[0]-x[0])
    sorted,iord = zip(*z)
    return sorted,iord

# Roughly equivalent to GRIDEA from zonav.f.  The name GRIDEA (GRID
# Equal Area) is itself an anchronism.  The GRIDEA in zonav.f does not
# compute an equal area grid, but was probably copied and modified from
# a version of the routine (for example in STEP3) that did.
def grid():
    """Return the latitude and other parameters that define the 14 belts
    (8 basic belts and 6 zones).  A triple of (*ibm*,*kzone*,*titlez*)
    is returned.  `ibm[b]` gives the number of boxes (regions) in belt
    *b* for `b in range(8)`.  *kzone* defines how the 6 combined zones
    are made from the basic belts.  `b in kzone[k]` is true when basic
    belt *b* is in special zone *k* (*b* is in range(8), *k* is in
    range(6)).  `titlez[z]` is an 80 character string that is the title
    for zone *z* for `z in range(14)`.
    """

    ibm = [1,2,3,4,4,3,2,1]
    # multiply by 4 to get the number of boxes in each band.
    ibm = map(lambda x: x*x, ibm)

    # This is not the same way the zones are defined in Fortran, but it
    # should produce the same result.
    # See doc/step5-notes
    N = set(range(4))
    G = set(range(8))
    S = G - N
    T = set([3,4])
    kzone = [N-T, T, S-T, N, S, G]

    # Accumulate the titles here.
    titlez = []
    # Boundaries (degrees latitude, +ve North) of the 8 basic belts.
    # This literal list is copied directly from zonav.f line 248, that's
    # why it looks ugly.
    band = [90.,64.2,44.4,23.6,0.,-23.6,-44.4,-64.2,-90.]
    for j in range(8):
        # southern and northern edges of the belt
        sedge = band[j+1]
        nedge = band[j]
        title = '  LATITUDE BELT FROM  XX.X ? TO  XX.X ?'
        # Convert title to a list of characters to edit it using slice
        # assignment.
        title = list(title)
        # Generally, when converting from Fortran string slice assignment to
        # Python's list slice assignment, we have to subtract 1 just from
        # the beginning of the slice.
        title[22:26] = '%4.1f' % abs(sedge)
        title[27] = 'SN'[sedge > 0]
        if sedge == 0:
            title[21:28] = 'EQUATOR'
        title[33:37] = '%4.1f' % abs(nedge)
        title[38] = 'SN'[nedge > 0]
        if nedge == 0:
           title[32:39] = 'EQUATOR'
        # Convert back to string
        title = ''.join(title)
        titlez.append(title)

    titlez += [
        '  NORTHERN LATITUDES: 23.6 N TO  90.0 N',
        '       LOW LATITUDES: 23.6 S TO  23.6 N',
        '  SOUTHERN LATITUDES: 90.0 S TO  23.6 S'
      ]
    # zonav.f line 74
    titlez += ['NORTHERN HEMISPHERE','SOUTHERN HEMISPHERE','GLOBAL MEANS']

    # Ensure all titles are 80 characters long.
    titlez = map(lambda s: ('%-80s' % s)[:80], titlez)
    
    return ibm, kzone, titlez


def main(argv=None):
    import sys

    if argv is None:
        argv = sys.argv
    # getopt here, if you need it.

    out = open('work/ZON.Ts.ho2.GHCN.CL.PA.1200.step1', 'wb')
    inp = open('result/BX.Ts.ho2.GHCN.CL.PA.1200', 'rb')
    log = open('log/zonav.Ts.ho2.GHCN.CL.PA.log', 'w')
    return zonav(inp, out, log)

if __name__ == '__main__':
    main()
