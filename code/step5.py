#!/usr/bin/env python
# $URL$
# $Rev$
# 
# step5.py
#
# David Jones, Ravenbrook Limited, 2009-10-27

"""
STEP5 of the GISTEMP algorithm.

In STEP5: 8000 subboxes are combined into 80 boxes, and ocean data is
combined with land data; boxes are combined into latitudinal zones
(including hemispheric and global zones); annual and seasonal anomalies
are computed from monthly anomalies.
"""

# Clear Climate Code
import eqarea
import fort
import earth

# http://www.python.org/doc/2.3.5/lib/itertools-functions.html
import itertools
# http://www.python.org/doc/2.3.5/lib/module-math.html
import math
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os
import sys

def SBBXtoBX(data, box, rland, intrp, base=(1961,1991), ignore_land=False):
    """Simultaneously combine the land series and the ocean series and
    combine subboxes into boxes.  *data* should be an iterator of
    (land, ocean) subbox series pairs; *box* should be an open binary
    output file (where the output from this routine will be written);
    *rland* and *intrp* are equivalent to the command line arguments
    to SBBXotoBX.f.  *ignore_land* is a flag that when True means cells
    with both ocean and land data will ignore their land data.
    """
    
    land_meta,ocean_meta = data.next()

    bos ='>'
    boxf = fort.File(box, bos=bos)

    # RCRIT in SBBXotoBX.f line 70
    # We make rland (and radius) a float so that any calculations done
    # using it are in floating point.
    radius = float(1200)
    if not ignore_land:
        rland = min(rland, radius)
    else:
        # Leave rland unassigned.  Any attempt to use it will fail (with
        # an UnboundLocalError), so that provides a check that we don't
        # use it.
        pass
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
    km2persubbox = (4*math.pi*earth.radius**2) / (nbox * nsubbox)

    novrlp=20

    # TODO: Formalise use of only monthlies, see step 3.
    assert land_meta.mavg == 6
    km = 12
    NYRSIN = land_meta.monm/km
    # IYRBGC in the Fortran code
    combined_year_beg = min(land_meta.yrbeg, ocean_meta.yrbeg)
    # Index into the combined array of the first year of the land data.
    land_offset = 12*(land_meta.yrbeg-combined_year_beg)
    # As I1TIN but for ocean data.
    ocean_offset = 12*(ocean_meta.yrbeg-combined_year_beg)
    # combined_n_months is MONMC in the Fortran.
    combined_n_months = max(land_meta.monm + land_offset,
                            land_meta.monm + ocean_offset)

    # Indices into the combined arrays for the base period (which is a
    # parameter).
    nfb = base[0] - combined_year_beg
    nlb = base[1] - combined_year_beg

    info = [land_meta.mo1, land_meta.kq, land_meta.mavg, land_meta.monm,
            land_meta.monm4, combined_year_beg, land_meta.missing_flag,
            land_meta.precipitation_flag]

    info[4] = 2 * land_meta.monm + 5
    boxf.writeline(struct.pack('%s8i' % bos, *info) + land_meta.title)

    # TODO: Use giss_data
    XBAD = land_meta.missing_flag

    # :todo: do we really need the area array to be 8000 cells long?
    for nr,box in enumerate(eqarea.grid()):
        # Averages for the land and ocean (one series per subbox)...
        avg = [[XBAD]*combined_n_months for _ in range(2*nsubbox)]
        area = [km2persubbox] * nsubbox
        wgtc = [0] * (nsubbox*2)
        # Eat the records from land and ocean 100 (nsubbox) at a time.
        # In other words, all 100 subboxes for the box (region).
        landsub,oceansub = zip(*itertools.islice(data, nsubbox))
        # :todo: combine below zip with above zip?
        for i,l,o in zip(range(nsubbox),landsub,oceansub):
            avg[i][land_offset:land_offset+len(l.series)] = l.series
            avg[i+nsubbox][ocean_offset:ocean_offset+len(o.series)] = o.series
            # Count the number of valid entries.
            wgtc[i] = l.good_count
            wgtc[i+nsubbox] = o.good_count
            # :ocean:weight:a: Assign a weight to the ocean cell.
            # A similar calculation appears elsewhere.
            if ignore_land:
                wocn = 0
            else:
                wocn = max(0, (landsub[i].d-rland)/(radius-rland))
            if wgtc[i+nsubbox] < 12*novrlp:
                wocn = 0
            # Normally *intrp* is 0
            if wocn > intrp:
                wocn = 1
            wgtc[i] *= (1 - wocn)
            wgtc[i+nsubbox] *= wocn

        # GISTEMP sort.
        # We want to end up with IORDR, the permutation array that
        # represents the sorter order.  IORDR[0] is the index (into the
        # *wgtc* array) of the longest record, IORDR[1] the index of the
        # next longest record, and so on.  We do that by decorating the
        # *wgtc* array with indexes 0 to 199, and then extracting the
        # (permuted) indexes into IORDR.
        # :todo: should probably import from a purpose built module.
        from step3 import sort
        z = zip(wgtc, range(2*nsubbox))
        # We only want to compare the weights (because we want to emulate
        # the GISTEMP sort exactly), so we use only the first part of the
        # tuple (that's why we can't just use `cmp`).
        sort(z, lambda x,y: y[0]-x[0])
        wgtc,IORDR = zip(*z)

        # From here to the "for" loop over the cells (below) we are
        # initialising data for the loop.  Primarily the AVGR and WTR
        # arrays.
        nc = IORDR[0]
        ncr = nc
        if ncr >= nsubbox:
            ncr = nc-nsubbox
        # :ocean:weight:b: Assign weight to ocean cell, see similar
        # calculation, above at :ocean:weight:a.
        # Line 191
        if ignore_land:
            wocn = 0
        else:
            wocn = max(0, (landsub[ncr].d - rland)/(radius-rland))
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
            if ignore_land:
                wocn = 0
            else:
                wocn = max(0, (landsub[ncr].d - rland)/(radius-rland))
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
            bias = tavg(avgr, km, NYRSIN, nfb, nlb, True, "Region %d" % nr)
        ngood = 0
        m = 0
        for iy in range(combined_n_months/km):
            for k in range(km):
                m = iy*km + k
                if avgr[m] == XBAD:
                    continue
                avgr[m] -= bias[k]
                ngood += 1

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

    # This is somewhat experimental.  *wt1* the weight of the incoming
    # data, *dnew*, can either be a scalar (applies to the entire *dnew*
    # series) or a list (in which case each element is the weight of the
    # corresponding element of *dnew*).  (zonav uses the list form).
    # In the body of this function we treat *wt1* as an indexable
    # object.  Here we convert scalars to an object that always returns
    # a constant.
    try:
        wt1[0]
        def update_bias():
            pass
    except TypeError:
        wt1_constant = wt1
        def update_bias():
            """Find mean bias."""
            wtmnew = wtm[k]+wt1_constant
            bias[k] = float(wtm[k]*bias[k]+wt1_constant*biask)/wtmnew
            wtm[k]=wtmnew
            return
        class constant_list:
            def __getitem__(self, i):
                return wt1_constant
        wt1 = constant_list()

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
        update_bias()

        # Update period of valid data, averages and weights
        for n in range(nf1, nl1):
            kn = k+km*n     # CSE for array index
            if dnew[kn] >= XBAD:
                continue
            wtnew = wt[kn] + wt1[kn]
            avg[kn] = float(wt[kn]*avg[kn] + wt1[kn]*(dnew[kn]+biask))/wtnew
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
def tavg(data, km, nyrs, base, limit, nr, id, deflt=0.0, XBAD=9999):
    """:meth:`tavg` computes the time averages (separately for each calendar
    month if *km*=12) over the base period (year *base* to *limit*) and
    saves them in *bias* (a fresh array that is returned).
    
    In case of no data, the average is set to
    *deflt* if nr=0 or computed over the whole period if nr>0.

    Similarly to :meth:`combine` *data* is treated as a linear array divided
    into years by considering contiguous chunks of *km* elements.

    *id* is an arbitrary printable value used for identification in
    diagnostic output (for example, the cell number).

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
        print "No data in base period - MONTH,NR,ID", k, nr, id
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


def zonav(inp):
    """ 
    Perform Zonal Averaging.

    The input *inp* is a box file, usually work/BX.Ts.GHCN.CL.PA.1200 .
    The data in the boxes are combined to produce averages over
    various latitudinal zones.  Returns an iterator of
    (averages, weights, title) tuples, one per zone.

    14 belts are produced.  The first 8 are the basic belts that are used
    for the equal area grid, the remaining 6 are combinations:

    (The Fortran uses 1 to 14, but here in Python we use 0 to 13)

      0 64N - 90N               \
      1 44N - 64N (asin 0.9)    |-  8 24N - 90 N  (0 + 1 + 2)
      2 24N - 44N (asin 0.7)    /
      3 Equ - 24N (asin 0.4)    \_  9 24S - 24 N  (3 + 4)
      4 24S - Equ               /
      5 44S - 24S               \
      6 64S - 44S               |- 10 90S - 24 S  (5 + 6 + 7)
      7 90S - 64S               /

     11 northern hemisphere (0 + 1 + 2 + 3)
     12 southern hemisphere (4 + 5 + 6 + 7)
     13 global (all belts 0 to 7)
    """

    inp = fort.File(inp, bos='>')

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
    yield (outinfo, titlei, ' zones:  90->64.2->44.4->23.6->0->-23.6->-44.4->-64.2->-90                      ')

    ibm,kzone,titlez = zones()

    JBM = len(ibm)

    # :init:lenz:
    lenz = [None]*JBM

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
            continue
        lenr,IORD = sort_perm(lenr)
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
            combine(avg[jb], bias, wt[jb], ar[nr], 0,nyrsin,
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

        yield (avg[jb], wt[jb], titlez[jb])

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
            combine(avgg, bias, wtg, avg[jb], 0,nyrsin,
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
        yield(avgg, wtg, titlez[JBM+jz])

def sort_perm(a):
    """The array *a* is sorted into descending order.  The fresh sorted
    array and the permutation array are returned as a pair (*sorted*,
    *iord*).  The original *a* is not mutated.

    The *iord* array is such that `a[iord[x]] == sorted[x]`.
    """
    from step3 import sort
    z = zip(a, range(len(a)))
    sort(z, lambda x,y: y[0]-x[0])
    sorted,iord = zip(*z)
    return sorted,iord

def zones():
    """Return the latitude and other parameters that define the 14 belts
    (8 basic belts and 6 zones).  A triple of (*ibm*,*kzone*,*titlez*)
    is returned.  `ibm[b]` gives the number of boxes (regions) in belt
    *b* for `b in range(8)`.  *kzone* defines how the 6 combined zones
    are made from the basic belts.  `b in kzone[k]` is true when basic
    belt *b* is in special zone *k* (*b* is in range(8), *k* is in
    range(6)).  `titlez[z]` is an 80 character string that is the title
    for zone *z* for `z in range(14)`.
    """

    # Number of boxes (regions) in each band.
    ibm = [4,8,12,16,16,12,8,4]

    N = set(range(4)) # Northern hemisphere
    G = set(range(8)) # Global
    S = G - N         # Southern hemisphere
    T = set([3,4])    # Tropics
    kzone = [N-T, T, S-T, N, S, G]

    # Boundaries (degrees latitude, +ve North) of the 8 basic belts.
    band = ['90.0 N',
            '64.2 N',
            '44.4 N',
            '23.6 N',
            'EQUATOR',
            '23.6 S',
            '44.4 S',
            '64.2 S',
            '90.0 S']
    # Accumulate the titles here.
    titlez = ['  LATITUDE BELT FROM %7s TO %7s' % (band[j+1], band[j]) for j in range(8)]

    titlez += [
        '  NORTHERN LATITUDES: 23.6 N TO  90.0 N',
        '       LOW LATITUDES: 23.6 S TO  23.6 N',
        '  SOUTHERN LATITUDES: 90.0 S TO  23.6 S',
        'NORTHERN HEMISPHERE',
        'SOUTHERN HEMISPHERE',
        'GLOBAL MEANS']

    # Ensure all titles are 80 characters long.
    titlez = map(lambda s: ('%-80s' % s)[:80], titlez)
    
    return ibm, kzone, titlez


def annzon(zoned_averages, out, zono, annzono, alternate={'global':2, 'hemi':True}):
    """Compute annual anomalies and write them out to various files.
    *zoned_averages* is an iterator of zoned averages produced by `zonav`.
    *out* is a list (length 4) of binary output files for the result
    files:
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

    bos = '>'

    zono = fort.File(zono, bos)
    annzono = fort.File(annzono, bos)

    # annzon.f line 41
    jzm = 14
    jzp = 14
    monmin = 6

    # annzon.f line 70
    (info, title, titl2) = zoned_averages.next()
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
    print >> out[0], ' Annual Temperature Anomalies (.01 C) - ' + title[28:80]
    for f in out[1:]:
        print >> f, title
    infoo = list(info)
    infoo[2] = 5
    infoo[3] //= 12
    # annzon.f line 42 and line 95 and following
    titleo = (
      'ANNUALLY AVERAGED (%d OR MORE MONTHS) TEMPERATURE ANOMALIES (C)' %
      monmin)
    # Ensure exactly 80 characters.
    titlelo = '%-80s' % titleo

    # annzon.f line 98
    # :read:zonal:
    # Collect JZM zonal means.
    w = len(struct.pack(bos + 'f', 0.0))
    for jz in range(jzm):
        (tdata, twt, titlez[jz]) = zoned_averages.next()
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
        x = int(math.floor(100*ann[z][iy] + 0.5))
        x = '%5d' % x
        if len(x) > 5:
            return '*****'
        return x

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

    # annzon.f line 44
    tit = ['    GLOBAL','N.HEMISPH.','S.HEMISPH.']
    # Shift the remaining 3 output files so that the indexing works out.
    out = out[1:]
    banner = 'Year   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug' + \
      '  Sep  Oct  Nov  Dec    J-D D-N    DJF  MAM  JJA  SON  Year'
    # All the "WRITE(96+J" stuff in the Fortran is replaced with this
    # enumeration into the *out* array (an array of file descriptors).
    # annzon.f line 191
    for j,outf in enumerate(out):
        print >> outf, (tit[j] + ' Temperature Anomalies' + 
          ' in .01 C     base period: 1951-1980')
        for iy in range(iy1tab-iyrbeg, iyrs):
            iout = [100*XBAD]*18
            if (iy+iyrbeg >= iy1tab+5 and ((iy+iyrbeg) % 20 == 1) or
              iy == iy1tab - iyrbeg):
                print >> outf
                print >> outf, banner
            # *data* for this zone, avoids some duplication of code.
            zdata = data[iord[j]]
            # :todo: Would probably be better to have a little 4-long
            # seasonal array to do the computation in.
            awin = 9999
            if iy > 0:
                awin = zdata[iy-1][11] + zdata[iy][0] + zdata[iy][1]
            aspr = sum(zdata[iy][2:5])
            asmr = sum(zdata[iy][5:8])
            afl  = sum(zdata[iy][8:11])
            if awin < 8000:
                iout[14] = nint(100.0*awin/3)
            if aspr < 8000:
                iout[15] = nint(100.0*aspr/3)
            if asmr < 8000:
                iout[16] = nint(100.0*asmr/3)
            if afl < 8000:
                iout[17] = nint(100.0*afl/3)
            # annzon.f line 216
            ann2=awin+aspr+asmr+afl
            if ann2 < 8000:
                iout[13] = nint(100.0*ann2/12)
            ann1=ann[iord[j]][iy]
            if iy == iyrs-1 and zdata[iy][-1] > 8000:
                ann1 = 9999
            if ann1 < 8000:
                iout[12] = nint(100.0*ann[iord[j]][iy])
            for m in range(12):
                iout[m] = nint(100.0*zdata[iy][m])
            iyr = iyrbeg+iy
            # Convert each of *iout* to a string, storing the results in
            # *sout*.
            sout = [None]*len(iout)
            for i,x in enumerate(iout):
                # All the elements of iout are formatted to width 5,
                # except for element 13 (14 in the Fortran code), which
                # is length 4.
                if i == 13:
                    x = '%4d' % x
                    if len(x) > 4:
                        x = '****'
                else:
                    x = '%5d' % x
                    if len(x) > 5:
                        x = '*****'
                sout[i] = x
            print >> outf, (
              '%4d ' + '%s'*12 + '  %s%s  ' + '%s'*4 + '%6d') % tuple(
              [iyr] + sout + [iyr])
        print >> outf, banner

    # annzon.f line 229
    # Save annual means on disk.
    annzono.writeline(struct.pack(bos+'8i', *infoo) +
                      titleo +
                      title[28:80] +
                      ' '*28)
    for jz in range(jzm):
        fmt = bos + '%df' % (monm//12)
        annzono.writeline(struct.pack(fmt, *ann[jz]) +
                          struct.pack(fmt, *annw[jz]) +
                          titlez[jz])
    # annzon.f line 236
    # Save monthly means on disk
    zono.writeline(struct.pack(bos + '8i', *info) +
                   title + titl2)

    for jz in range(jzm):
        fmt = bos + '%df' % monm
        zono.writeline(struct.pack(fmt, *itertools.chain(*data[jz])) +
                       struct.pack(fmt, *itertools.chain(*wt[jz])) +
                       titlez[jz])


def nint(x):
    """Nearest integer.  Reasonable approximation to Fortran's NINT
    routine."""

    # http://www.python.org/doc/2.4.4/lib/module-math.html

    return int(math.floor(x+0.5))


def step5(data):
    """Step 5 of GISTEMP.

    This step takes input provided by steps 3 and 4 (zipped together).

    :Param data:
        These are the land and ocean sub-box data, zipped together.
        They need to support the protocol defined by
        `code.giss_data.SubboxSetProtocol`.

    """
    box = open(os.path.join('result', 'BX.Ts.ho2.GHCN.CL.PA.1200'), 'wb')
    SBBXtoBX(data, box, rland=100, intrp=0)
    # Necessary, because box is an input to the next stage, so the file
    # needs to be fully written.
    box.close()

    inp = open(os.path.join('result', 'BX.Ts.ho2.GHCN.CL.PA.1200'), 'rb')
    zoned_averages = zonav(inp)

    out = ['ZonAnn', 'GLB', 'NH', 'SH']
    out = [open(os.path.join('result', bit+'.Ts.ho2.GHCN.CL.PA.txt'), 'w')
            for bit in out]
    zono = open(os.path.join('work', 'ZON.Ts.ho2.GHCN.CL.PA.1200'), 'wb')
    annzono = open(os.path.join('work', 'ANNZON.Ts.ho2.GHCN.CL.PA.1200'), 'wb')
    annzon(zoned_averages, out, zono, annzono)

    yield "Step 5 Completed"
