#!/usr/bin/env python
# $URL$
# $Rev$
#
# annzon.py
#
# David Jones, Ravenbrook Limited, 2009-12-24

"""
annzon.py computes annual anomalies from monthly anomalies (and also
recomputes global and hemispheric averages).

As input it takes the output of zonav.py (a ZON file).

This is part of GISTEMP STEP5.

This is a work in progress (output of ZonAnn.* file is correct).
"""

def annzon(inp, log, out, zono, annzono,
  alternate={'global':2, 'hemi':True}):
    """Compute annual anomalies and write them out to various files.
    *inp* is the input file, a ZON file typically produced by `zonav`.
    *log* is an output file for logging.  *out* is a list (length 4) of
    binary output files for the result files:
    (ZonAnn|GLB|NH|SH).Ts.ho2.GHCN.CL.PA.txt .
    *zono* and *annzono* are binary output files where the (recomputed)
    zonal means and annual zonal means are written.

    The *alternate* argument (replaces IAVGGH in the Fortran) controls
    whether alternate algorithms are used to compute the global and
    hemispheric means (see annzon.f line 123 to line 162).
    alternate['global'] is 1 or 2, to select 1 of 2 different alternate
    computations, or false to not compute an alternative;
    alternate['hemi'] is true to compute an alternative, false
    otherwise.
    """

    # Clear Climate Code
    import fort
    # http://www.python.org/doc/2.4.4/lib/module-struct.html
    import struct

    bos = '>'

    inp = fort.File(inp, bos)
    zono = fort.File(zono, bos)
    annzono = fort.File(annzono, bos)

    # annzon.f line 41
    jzm = 14
    jzp = 14
    monmin = 6

    # annzon.f line 70
    header = struct.unpack(bos+'8i80s80s', inp.readline())
    info = header[:8]
    title = header[8]
    titl2 = header[9]
    kq = info[1]
    # annzon.f line 72
    # km: the number of time frames per year.
    km = 1
    if info[2] == 6:
        km = 12
    if info[2] == 7:
        km = 4
    iyrbeg = info[5]
    monm = info[3]
    iyrs = monm // km
    # Allocate the 2- and 3- dimensional arrays.
    # The *data* and *wt* arrays are properly 3 dimensional
    # ([zone][year][month]), but the inner frames are only allocated
    # when we read the data, see :read:zonal below.
    # annzon.f line 79
    data = [ None for _ in range(jzm)]
    wt =   [ None for _ in range(jzm)]
    ann =  [ [None]*iyrs for _ in range(jzm)]
    annw = [ [None]*iyrs for _ in range(jzm)]
    titlez = [None]*jzm
    # Here we use the Python convention, *iyrend* is one past the highest
    # year used.
    iyrend = info[3] // km + iyrbeg
    XBAD = float(info[6])
    # annzon.f line 84
    # Create and write out the header record of the output files.
    print >> log, title
    print >> out[0], ' Annual Temperature Anomalies (.01 C) - ' + title[28:80]
    for f in out[1:]:
        print >> f, title
    print >> log, info
    infoo = list(info)
    infoo[2] = 5
    infoo[3] //= 12
    # annzon.f line 42 and line 95 and following
    titleo = (
      'ANNUALLY AVERAGED (%d OR MORE MONTHS) TEMPERATURE ANOMALIES (C)' %
      monmin)
    # Ensure exactly 80 characters.
    titlelo = '%-80s' % titleo
    # Ah, the Fortran emits a '1' at the beginning of the title record
    # on the log file; this is an ANSI printer escape (for bold, if I
    # remember my IBM days rightly).  We Don't Do That.
    print >> log, titleo
    print >> log, infoo

    # annzon.f line 98
    # :read:zonal:
    # Collect JZM zonal means.
    fmt = bos + '%df' % monm
    w = len(struct.pack(bos + 'f', 0.0))
    for jz in range(jzm):
        record = inp.readline()
        tdata = struct.unpack(fmt, record[:w*monm])
        twt = struct.unpack(fmt, record[w*monm:2*w*monm])
        assert 80 == len(record[2*w*monm:])
        titlez[jz] = record[2*w*monm:]
        # Regroup the *data* and *wt* series so that they come in blocks of
        # 12 (*km*, really).
        # Uses essentially the same trick as the `grouper()` recipe in
        # http://docs.python.org/library/itertools.html#recipes
        data[jz] = zip(*[iter(tdata)]*km)
        wt[jz] = zip(*[iter(twt)]*km)

    # annzon.f line 104
    # Find (compute) the annual means.
    for jz in range(jzm):
        for iy in range(iyrs):
            ann[jz][iy] = XBAD
            annw[jz][iy] = 0.
            anniy = 0.
            annwiy = 0.
            mon = 0
            # annzon.f line 114
            for m in range(km):
                if data[jz][iy][m] == XBAD:
                    continue
                mon += 1
                anniy += data[jz][iy][m]
                annwiy += wt[jz][iy][m]
            if mon >= monmin:
                ann[jz][iy] = float(anniy)/mon
            # :todo: Surely KM is intended here, not 12?
            annw[jz][iy] = annwiy/12.

    # annzon.f line 123
    # Alternate global mean.
    if alternate['global']:
        glb = alternate['global']
        assert glb in (1,2)
        # Pick which "four" zones to use.
        # annzon.f line 50
        # (subtracting 1 from each zone to convert to Python convention)
        if glb == 1:
            zone = [8, 9, 9, 10]
        else:
            zone = [8, 3, 4, 10]
        wtsp = [3.,2.,2.,3.]
        # annzon.f line 128
        for iy in range(iyrs):
            glob = 0.
            # Fortran uses JZM to index the last zone, we use -1.
            ann[-1][iy] = XBAD
            for z,w in zip(zone, wtsp):
                if ann[z][iy] == XBAD:
                    # Note: Rather ugly use of "for...else" to emulate GOTO.
                    break
                glob += ann[z][iy]*w
            else:
                ann[-1][iy] = .1 * glob
        # annzon.f line 137
        for iy in range(iyrs):
            data[-1][iy] = [XBAD]*12
            for m in range(12):
                glob = 0.
                for z,w in zip(zone, wtsp):
                    if data[z][iy][m] == XBAD:
                        break
                    glob += data[z][iy][m]*w
                else:
                    data[-1][iy][m] = .1 * glob
    # annzon.f line 147
    # Alternate hemispheric means.
    if alternate['hemi']:
        # Note: In Fortran IHEM is either 1 or 2, in Python it is 0 or
        # 1.  This changes some computations but not others (because the
        # indexing is from 0 in Python).
        # For the computations it will be useful to recall how the zones
        # are numbered.  There is a useful docstring at the beginning of
        # zonav.py.
        for ihem in range(2):
            for iy in range(iyrs):
                ann[ihem+11][iy] = XBAD
                if ann[ihem+3][iy] != XBAD and ann[2*ihem+8][iy] != XBAD:
                    ann[ihem+11][iy] = (0.4*ann[ihem+3][iy] +
                                        0.6*ann[2*ihem+8][iy])
            for iy in range(iyrs):
                data[ihem+11][iy] = [XBAD]*12
                for m in range(12):
                    if (data[ihem+3][iy][m] != XBAD and
                      data[2*ihem+8][iy][m] != XBAD):
                        data[ihem+11][iy][m] = (
                          0.4*data[ihem+3][iy][m] +
                          0.6*data[2*ihem+8][iy][m])

    # iord literal borrowed exactly from Fortran...
    # annzon line 48
    iord = [14,12,13, 9,10,11, 1,2,3,4,5,6,7,8]
    # ... and then adjusted for Python index convention.
    iord = map(lambda x: x-1, iord)
    # annzon.f line 38
    iy1tab = 1880
    # annzon.f line 163
    # Display the annual means.
    def annasstr(z):
        """Helper function that returns the annual anomaly for zone *z*
        as a string representation of an integer (the integer is the
        anomaly scaled by 100 to convert to centikelvin).

        The returned value is a string that is 5 characters long.  If
        the integer will not fit into a 5 character string, '*****' is
        returned (this emulates the Fortran convention of formatting
        999900 (which is the XBAD value in centikelvin) as a '*****'.
        
        The year, *iy*, is lexically captured which is a bit horrible.
        """
        # http://www.python.org/doc/2.4.4/lib/module-math.html
        import math

        x = int(math.floor(100*ann[z][iy] + 0.5))
        x = '%5d' % x
        if len(x) > 5:
            return '*****'
        return x

    # annzon.f line 166
    for jz in range(jzp):
        print >> log, jz, titlez[iord[jz]]
    iyrsp = iyrs
    # Check (and skip) incomplete year.
    if data[-1][-1][-1] > 8000:
        iyrsp -= 1
    banner = """
                           24N   24S   90S     64N   44N   24N   EQU   24S   44S   64S   90S
Year  Glob  NHem  SHem    -90N  -24N  -24S    -90N  -64N  -44N  -24N  -EQU  -24S  -44S  -64S Year
""".strip('\n')
    for iy in range(iy1tab - iyrbeg, iyrsp):
        if (iy+iyrbeg >= iy1tab+5 and ((iy+iyrbeg) % 20 == 1) or
          iy == iy1tab - iyrbeg):
            print >> out[0]
            print >> out[0], banner
        iyr = iyrbeg+iy
        print >> out[0], ('%4d' + ' %s'*3 + '  ' + ' %s'*3 +
                          '  ' + ' %s'*8 + '%5d') % tuple([iyr] +
          [annasstr(iord[zone]) for zone in range(jzp)] + [iyr])
    # The trailing banner is just like the repeated banner, except that
    # "Year  Glob  NHem  SHem" appears on on the first line, not the
    # second line (and the same for the "Year" that appears at the end
    # of the line).  *sigh*.
    banner = banner.split('\n')
    banner[0] = banner[1][:24] + banner[0][24:] + ' Year'
    banner[1] = ' '*24 + banner[1][24:-5]
    banner = '\n'.join(banner)

    print >> out[0], banner
    print >> out[0]


def main():
    # http://www.python.org/doc/2.4.4/lib/module-os.html
    import os

    inp = open(os.path.join('work', 'ZON.Ts.ho2.GHCN.CL.PA.1200.step1'), 'rb')
    log = open(os.path.join('log', 'annzon.Ts.ho2.GHCN.CL.PA.log'), 'w')
    out = ['ZonAnn', 'GLB', 'NH', 'SH']
    out = [open(os.path.join('result', bit+'.Ts.ho2.GHCN.CL.PA.txt'), 'w')
            for bit in out]
    zono = open(os.path.join('work', 'ZON.Ts.ho2.GHCN.CL.PA.1200'), 'wb')
    annzono = open(os.path.join('work', 'ANNZON.Ts.ho2.GHCN.CL.PA.1200'), 'wb')

    return annzon(inp, log, out, zono, annzono)

if __name__ == '__main__':
    main()
