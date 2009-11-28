#!/usr/bin/env python
# $URL$
# $Rev$
# 
# step5.py
# 
# Python code reproducing the STEP5 part of the GISTEMP algorithm.
#
# David Jones, Ravenbrook Limited, 2009-10-27
#
# Work in progress.  Currently it replaces SBBXtoBX (but not yet 
#
# The code is derived from the Fortran GISTEMP code.
#
# Python notes:
#
# Much of the binary IO uses the struct module.  A quick review will
# help.
#
# I feel obliged to tell you that the Python expression N * list gives a
# new list that conists of list repeated N times.  This is used every
# now and then.
#
# We also make use of list slice assignment: a[2:4] = range(2)

# Ravenbrook
import eqarea
# Ravenbrook
import fort
# Ravenbrook
import subbox

# http://www.python.org/doc/2.3.5/lib/itertools-functions.html
import itertools
# http://www.python.org/doc/2.3.5/lib/module-math.html
import math
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct

# See SBBXotoBX.f
def SBBXtoBX(land, ocean, box, log, rland, intrp, base=(1961,1991)):
    """Simultaneously combine the land series and the ocean series and
    combine subboxes into boxes.  *land* and *ocean* should be open
    binary files for the land and ocean series (subbox gridded); *box*
    should be an open binary output file (where the output from this
    routine will be written); *log* should be an open text file (where
    diagnostic information will be logged); *rland* and *intrp* are
    equivalent to the command line arguments to SBBXotoBX.f.
    """

    bos ='>'
    boxf = fort.File(box, bos=bos)

    # RCRIT in SBXotoBX.f line 70
    radius = float(1200)
    rland = min(rland, radius)
    if rland < 0:
        rland = -9999
    # We make rland (and radius) a float so that any calculations done
    # using it are in floating point.
    rland = float(rland)
    # clamp intrp to [0,1]
    intrp = max(min(intrp, 1), 0)

    # Danger! Suspiciously similar to step3.py (because it's pasted):

    # number of boxes
    nbox = 80
    # number of subboxes within each box
    nsubbox = 100

    # Much of the Fortran code assumes that various "parameters" have
    # particular fixed values (probably accidentally).  I don't trust the
    # Python code to avoid similar assumptions.  So assert the
    # "parameter" values here.  Just in case anyone tries changing them.
    # Note to people reading comment becuse the assert fired:  Please
    # don't assume that the code will "just work" when you change one of
    # the parameter values.  It's supposed to, but it might not.
    assert nbox == 80
    # NCM in the Fortran
    assert nsubbox == 100

    # area of subbox in squared kilometres
    # Recall area of sphere is 4*pi*(r**2)
    # TOMETR in SBBXotoBX.f line 99
    import earth
    km2persubbox = (4*math.pi*earth.radius**2) / (nbox * nsubbox)

    novrlp=20

    landsubbox = subbox.File(land)
    oceansubbox = subbox.File(ocean)
    km = 1
    if landsubbox.mavg == 6:
        km = 12
    if landsubbox.mavg == 7:
        km = 4
    NYRSIN = landsubbox.monm/km
    # IYRBGC in the Fortran code
    combined_year_beg = min(landsubbox.iyrbeg, oceansubbox.iyrbeg)
    # index into the combined array of the first year of the land data.
    I1TIN = 12*(landsubbox.iyrbeg-combined_year_beg)
    # as I1TIN but for ocean data
    I1TINO = 12*(oceansubbox.iyrbeg-combined_year_beg)
    # combined_n_months is MONMC in the Fortran.
    combined_n_months = max(landsubbox.monm + I1TIN,
                            landsubbox.monm + I1TINO)
    # 
    # Indices into the combined arrays for the base period (which is a
    # parameter).
    nfb = base[0] - combined_year_beg
    nlb = base[1] - combined_year_beg

    info = landsubbox.info()
    info[3] = landsubbox.monm
    info[4] = 2*landsubbox.monm+5
    info[5] = combined_year_beg
    boxf.writeline(struct.pack('%s8i' % bos, *info) + landsubbox.title)

    XBAD = landsubbox.missing_flag

    # :todo: do we really need the area array to be 8000 cells long?
    for nr,box in enumerate(eqarea.grid()):
        # Averages for the land and ocean (one series per subbox)...
        avg = [[XBAD]*combined_n_months for _ in range(2*nsubbox)]
        area = [km2persubbox] * nsubbox
        wgtc = [0] * (nsubbox*2)
        # Eat the records from land and ocean 100 (nsubbox) at a time.
        # In other words, all 100 subboxes for the box (region).
        landsub = list(itertools.islice(landsubbox, nsubbox))
        oceansub = list(itertools.islice(oceansubbox, nsubbox))
        for i,l,o in zip(range(nsubbox),landsub,oceansub):
            l = tuple(l)
            o = tuple(o)
            avg[i][I1TIN:I1TIN+len(l)] = l
            avg[i+nsubbox][I1TINO:I1TINO+len(o)] = o
            # Count the number of valid entries.
            # :todo: perhaps do something fancy with ".count()"
            for x in l:
                if x != XBAD:
                    wgtc[i] += 1
            for x in o:
                if x != XBAD:
                    wgtc[i+nsubbox] += 1
            # :ocean:weight:a: Assign a weight to the ocean cell.
            # A similar calculation appears elsewhere.
            wocn = max(0, (landsub[i].d-rland)/(radius-rland))
            # Line 174.  Isn't normally so.
            if rland == XBAD:
                wocn = 0
            if wgtc[i+nsubbox] < 12*novrlp:
                wocn = 0
            # Normally *intrp* is 0
            if wocn > intrp:
                wocn = 1
            # The following code appears to ponderously multiply wgtc[i]
            # by either wocn or (1-wocn) as appropriate.
            # Line 177
            wgtc[i] = 0
            wgtc[i+nsubbox] = 0
            for x in l:
                if x != XBAD:
                    wgtc[i] += (1-wocn)
            for x in o:
                if x != XBAD:
                    wgtc[i+nsubbox] += wocn

        # GISTEMP sort.
        # We want to end up with IORDR, the permutation array that
        # represents the sorter order.  IORDR[0] is the index (into the
        # *wgtc* array) of the longest record, IORDR[1] the index of the
        # next longest record, and so on.  We do that by decorating the
        # *wgtc* array with indexes 0 to 199, and then extracting the
        # (permuted) indexes into IORDR.
        # :todo: should probably import from a purpopse built module.
        from step3 import sort
        z = zip(wgtc, range(2*nsubbox))
        # We only want to compare the weights (because we want to emulate
        # the GISTEMP sort exactly), so we use only the first part of the
        # tuple (that's why we can't just use `cmp`).
        sort(z, lambda x,y: y[0]-x[0])
        wgtc,IORDR = zip(*z)

        # From here to the for loop over the cells (below) we are
        # initialising data for the loop.  Primarily the AVGR and WTR
        # arrays.
        nc = IORDR[0]
        ncr = nc
        if ncr >= nsubbox:
            ncr = nc-nsubbox
        # :ocean:weight:b: Assign weight to ocean cell, see similar
        # calculation, above at :ocean:weight:a.
        # Line 191
        wocn = max(0, (landsub[ncr].d - rland)/(radius-rland))
        if rland == XBAD:
            wocn = 0
        if nc == ncr and wgtc[nc + nsubbox] < 12*novrlp:
            wocn = 0
        if wocn > intrp:
            wocn = 1
        wnc = wocn
        if nc < nsubbox:
            wnc = 1 - wocn

        # line 197
        # area is assumed to be constant, so we don't use it here.
        # (which actually avoids a bug that is present in the Fortran, see
        # doc/step5-notes for details).
        wtm = [wnc]*km
        bias = [0]*km
        
        # Weights for the region's record.
        wtr = [0]*combined_n_months
        for m,a in enumerate(avg[nc]):
            if a < XBAD:
                wtr[m] = wnc
        # Create the region (box) record by copying the subbox record
        # into AVGR
        avgr = avg[nc][:]

        # Loop over the remaining cells.
        for n in range(1,2*nsubbox):
            nc = IORDR[n]
            w = wgtc[n]
            # print "nr %(nr)d, n %(n)d, nc %(nc)d, wgtc %(w)d" % locals()
            # Line 207
            # :todo: Can it be correct to use [n]?  It's what the
            # Fortran does.
            if wgtc[n] < 12*novrlp:
                continue
            ncr = nc
            if nc >= nsubbox:
                ncr = nc - nsubbox
            wocn = max(0, (landsub[ncr].d - rland)/(radius-rland))
            if rland == XBAD:
                wocn = 0
            if nc == ncr and wgtc[nc + nsubbox] < 12*novrlp:
                wocn = 0
            if wocn > intrp:
                wocn = 1
            wnc = wocn
            if nc < nsubbox:
                wnc = 1 - wocn
            wt1 = wnc
            nsm = combine(avgr, bias, wtr, avg[nc], 0, combined_n_months/km,
              wt1, wtm, km, nc)
        if nfb > 0:
            bias = tavg(avgr, km, NYRSIN, nfb, nlb, nr, 0)
        ngood = 0
        m = 0
        for iy in range(combined_n_months/km):
            for k in range(km):
                m = iy*km + k
                if avgr[m] == XBAD:
                    continue
                avgr[m] -= bias[k]
                ngood += 1

        def nint10x(x):
            return int((10*x) + 0.5)
        def nint(x):
            return int(x+0.5)
        for l in [(map(nint10x, avgr[i:i+12]),map(nint10x, avgr[i+12:i+24]))
                    for i in range(0, combined_n_months, 24)]:
            print >> log, l
        # GISTEMP divdes the weight by TOMETR, but we never multiplied it by
        # the area in the first place.
        for l in [(map(nint, wtr[i:i+12]),map(nint, wtr[i+12:i+24]))
                    for i in range(0, combined_n_months, 24)]:
            print >> log, l
        # MONM or MONMC?
        fmt = '%s%df' % (bos, combined_n_months)
        boxf.writeline(
          struct.pack(fmt, *avgr) +
          struct.pack(fmt, *wtr) +
          struct.pack('%si' % bos, ngood) +
          struct.pack('%s4i' % bos, *box)
          )


# :todo: This was nabbed from code/step3.py.  Put it in one place and
# make it common.
#
# Equivalent of the Fortran subroutine CMBINE
# Parameters as per Fortran except that nsm is returned directly.
def combine(avg, bias, wt, dnew, nf1, nl1, wt1, wtm, km, id,
  NOVRLP=20, XBAD=9999):
    """Run the GISTEMP combining algorithm.  This combines the data
    in the *dnew* array into the *avg* array (also updating the *bias*
    array).

    Each of the arguments *avg*, *wt*, *dnew* is a linear array that is
    divided into "years" by considering each contiguous segment of
    *km* elements a year.  Only data for years in range(nf1, nl1) are
    considered and combined.  Note that range(nf1, nl1) includes *nf1*
    but excludes *nl1* (and that this differs from the Fortran
    convention).
    
    Each month (or other subdivision, such as season, according to
    *km*) of the year is considered separately.  For the set of times
    where both *avg* and *dnew* have data the mean difference (a bias)
    is computed.  If there are fewer than *NOVRLP* years in common the
    data (for that month of the year) are not combined.  The bias is
    subtracted from the *dnew* record and it is point-wise combined
    into *avg* according to the weight *wt1* and the exist
    weight for *avg*.

    *id* is an identifier used only when diagnostics are issued
    (when combining stations it is expected to be the station ID; when
    combining subboxes it is expected to be the subbox number (0 to 99)).
    """

    # In the absence of type checks, check that the arrays have an
    # accessible element.
    avg[0]
    bias[km-1]
    wt[0]
    dnew[0]
    wtm[km-1]
    assert nf1 < nl1
    assert wt1 >= 0


    # See to.SBBXgrid.f lines 519 and following

    # The Fortran code handles the arrays by just assuming them to
    # be 2D arrays of shape (*,KM).  Sadly Python array handling
    # just isn't that convenient, so look out for repeated uses of
    # "[k+km*n]" instead.

    nsm = 0
    missed = km
    missing = [True]*km
    for k in range(km):
        sumn = 0    # Sum of data in dnew
        sum = 0     # Sum of data in avg
        ncom = 0    # Number of years where both dnew and avg are valid
        for n in range(nf1, nl1):
            kn = k+km*n     # CSE for array index
            # Could specify that arguments are array.array and use
            # array.count(BAD) and sum, instead of this loop.
            if avg[kn] >= XBAD or dnew[kn] >= XBAD:
                continue
            ncom += 1
            sum += avg[kn]
            sumn += dnew[kn]
        if ncom < NOVRLP:
            continue
        biask = float(sum-sumn)/ncom
        # Find mean bias
        wtmnew = wtm[k]+wt1
        # :todo: remove after debugging.
        if wtmnew == 0:
            print wtm, wt1, k, wtmnew
            l = dict(locals())
            del l['dnew']
            del l['wt']
            del l['avg']
            print l
        bias[k] = float(wtm[k]*bias[k]+wt1*biask)/wtmnew
        wtm[k]=wtmnew
        # Update period of valid data, averages and weights
        for n in range(nf1, nl1):
            kn = k+km*n     # CSE for array index
            if dnew[kn] >= XBAD:
                continue
            wtnew = wt[kn] + wt1
            avg[kn] = float(wt[kn]*avg[kn] + wt1*(dnew[kn]+biask))/wtnew
            wt[kn] = wtnew
            nsm += 1
        missed -= 1
        missing[k] = False
    if False and missed > 0:
        print "Unused data - ID/SUBBOX,WT", id, wt1, missing
    return nsm

# :todo: make common with step3.py
# Equivalent to Fortran subroutine TAVG.  Except the bias array
# (typically 12 elements) is returned.
# See to.SBBXgrid.f lines 563 and following.
def tavg(data, km, nyrs, base, limit, nr, nc, deflt=0.0, XBAD=9999):
    """:meth:`tavg` computes the time averages (separately for each calendar
    month if *km*=12) over the base period (year *base* to *limit*) and
    saves them in *bias*. In case of no data, the average is set to
    *deflt* if nr=0 or computed over the whole period if nr>0.

    Similarly to :meth:`combine` *data* is treated as a linear array divided
    into years by considering contiguous chunks of *km* elements.

    Note: the Python convention for *base* and *limit* is used, the base
    period consists of the years starting at *base* and running up to,
    but including, the year *limit*.
    """

    bias = [0.0]*km
    missed = km
    len = km*[0]    # Warning: shadows builtin "len"
    for k in range(km):
        bias[k] = deflt
        sum = 0.0
        m = 0
        for n in range(base, limit):
            kn = k+km*n     # CSE for array index
            if data[kn] >= XBAD:
                continue
            m += 1
            sum += data[kn]
        len[k] = m
        if m == 0:
            continue
        bias[k] = float(sum)/float(m)
        missed -= 1
    if nr*missed == 0:
        return bias
    # Base period is data free (for at least one month); use bias
    # with respect to whole series.
    for k in range(km):
        if len[k] > 0:
            continue
        print "No data in base period - MONTH,NR,NC", k, nr, nc
        sum = 0.0
        m = 0
        for n in range(nyrs):
            kn = k+km*n     # CSE for array index
            if data[kn] >= XBAD:
                continue
            m += 1
            sum += data[kn]
        if m == 0:
            continue
        bias[k] = float(sum)/float(m)
    return bias


def step5():
    land = open('work/SBBX1880.Ts.GHCN.CL.PA.1200', 'rb')
    ocean = open('input/SBBX.HadR2', 'rb')
    box = open('result/BX.Ts.ho2.GHCN.CL.PA.1200', 'wb')
    log = open('log/SBBXotoBX.log', 'w')
    SBBXtoBX(land, ocean, box, log, rland=100, intrp=0)
    # 
    # zonav
    # annzon

def main(argv=None):
    import sys

    if argv is None:
        argv = sys.argv
    return step5()

if __name__ == '__main__':
    main()
