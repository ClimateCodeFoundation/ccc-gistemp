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
import giss_data
import parameters
from tool import giss_io

# http://www.python.org/doc/2.3.5/lib/itertools-functions.html
import itertools
# http://www.python.org/doc/2.3.5/lib/module-math.html
import math
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os
import sys

def SBBXtoBX(data):
    """Simultaneously combine the land series and the ocean series and
    combine subboxes into boxes.  *data* should be an iterator of
    (land, ocean) subbox series pairs. Returns an iterator of box data.
    """
    
    # First item from iterator is normally a pair of metadataobjects,
    # one for land, one for ocean.  If we are piping step3 straight into
    # step5 then it is not a pair.  In that case we synthesize missing
    # ocean data.
    meta = data.next()
    try:
        land_meta,ocean_meta = meta
    except:
        # Use the land meta object for both land and ocean data
        land_meta,ocean_meta = meta,meta
        print "No ocean data; using land data only"
        def blank_ocean_data(data):
            """Augment a land-only data series with blank ocean data."""
            for land_box in data:
                series = [giss_data.XMISSING] * len(land_box.series)
                ocean_box = giss_data.SubboxRecord(
                    lat_S=land_box.lat_S,
                    lat_N=land_box.lat_N,
                    lon_W=land_box.lon_W,
                    lon_E=land_box.lon_E,
                    stations=0, station_months=0,
                    d=giss_data.XMISSING, series=series)
                yield land_box, ocean_box
        data = blank_ocean_data(data)

    # number of subboxes within each box
    nsubbox = 100

    # TODO: Formalise use of only monthlies, see step 3.
    assert land_meta.mavg == 6
    NYRSIN = land_meta.monm/12
    combined_year_beg = min(land_meta.yrbeg, ocean_meta.yrbeg)
    # Index into the combined array of the first year of the land data.
    land_offset = 12*(land_meta.yrbeg-combined_year_beg)
    # As I1TIN but for ocean data.
    ocean_offset = 12*(ocean_meta.yrbeg-combined_year_beg)
    combined_n_months = max(land_meta.monm + land_offset,
                            land_meta.monm + ocean_offset)

    info = [land_meta.mo1, land_meta.kq, land_meta.mavg, land_meta.monm,
            land_meta.monm4, combined_year_beg, land_meta.missing_flag,
            land_meta.precipitation_flag]

    info[4] = 2 * land_meta.monm + 5
    yield(info, land_meta.title)

    # TODO: Use giss_data
    XBAD = land_meta.missing_flag

    for box_number,box in enumerate(eqarea.grid()):
        # Averages for the land and ocean (one series per subbox)...
        avg = [[XBAD]*combined_n_months for _ in range(2*nsubbox)]
        wgtc = [0] * (nsubbox*2)
        # Eat the records from land and ocean 100 (nsubbox) at a time.
        # In other words, all 100 subboxes for the box.
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
            if wgtc[i+nsubbox] < parameters.subbox_min_valid or landsub[i].d < parameters.subbox_land_range:
                wocn = 0
            else:
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
        if landsub[ncr].d < parameters.subbox_land_range or (nc == ncr and wgtc[nc+nsubbox] < parameters.subbox_min_valid):
            wocn = 0
        else:
            wocn = 1

        wnc = wocn
        if nc < nsubbox:
            wnc = 1 - wocn

        wtm = [wnc]*12
        bias = [0]*12
        
        # Weights for the box's record.
        wtr = [0]*combined_n_months
        for m,a in enumerate(avg[nc]):
            if a < XBAD:
                wtr[m] = wnc
        # Create the box record by copying the subbox record
        # into AVGR
        avgr = avg[nc][:]

        # Loop over the remaining cells.
        for n in range(1,2*nsubbox):
            nc = IORDR[n]
            w = wgtc[n]
            # :todo: Can it be correct to use [n]?  It's what the
            # Fortran does.
            if wgtc[n] < parameters.subbox_min_valid:
                continue
            ncr = nc
            if nc >= nsubbox:
                ncr = nc - nsubbox
            if landsub[ncr].d < parameters.subbox_land_range or (nc == ncr and wgtc[nc + nsubbox] < parameters.subbox_min_valid):
                wocn = 0
            else:
                wocn = 1
            wnc = wocn
            if nc < nsubbox:
                wnc = 1 - wocn
            wt1 = wnc
            nsm = combine(avgr, bias, wtr, avg[nc], 0, combined_n_months/12, wt1, wtm, nc)
        bias = tavg(avgr, NYRSIN,
                    parameters.subbox_reference_first_year - combined_year_beg,
                    parameters.subbox_reference_last_year - combined_year_beg + 1,
                    True, "Box %d" % box_number)
        ngood = 0
        m = 0
        for iy in range(combined_n_months/12):
            for k in range(12):
                m = iy*12 + k
                if avgr[m] == XBAD:
                    continue
                avgr[m] -= bias[k]
                ngood += 1

        yield (avgr, wtr, ngood, box)


# :todo: This was nabbed from code/step3.py.  Put it in one place and
# make it common.
#
def combine(avg, bias, wt, dnew, nf1, nl1, wt1, wtm, id, XBAD=9999):
    """Run the GISTEMP combining algorithm.  This combines the data
    in the *dnew* array into the *avg* array (also updating the *bias*
    array).

    Each of the arguments *avg*, *wt*, *dnew* is a linear array that
    is divided into "years" by considering each contiguous segment of
    12 elements a year.  Only data for years in range(nf1, nl1) are
    considered and combined.  Note that range(nf1, nl1) includes *nf1*
    but excludes *nl1* (and that this differs from the Fortran
    convention).
    
    Each month of the year is considered separately.  For the set of
    times where both *avg* and *dnew* have data the mean difference (a
    bias) is computed.  If there are fewer than
    *parameters.box_min_overlap* years in common the data (for that
    month of the year) are not combined.  The bias is subtracted from
    the *dnew* record and it is point-wise combined into *avg*
    according to the weight *wt1* and the exist weight for *avg*.

    *id* is an identifier used only when diagnostics are issued
    (when combining stations it is expected to be the station ID; when
    combining subboxes it is expected to be the subbox number (0 to 99)).
    """

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

    nsm = 0
    missed = 12
    missing = [True]*12
    for k in range(12):
        sumn = 0    # Sum of data in dnew
        sum = 0     # Sum of data in avg
        ncom = 0    # Number of years where both dnew and avg are valid
        for n in range(nf1, nl1):
            kn = k+12*n     # CSE for array index
            # Could specify that arguments are array.array and use
            # array.count(BAD) and sum, instead of this loop.
            if avg[kn] >= XBAD or dnew[kn] >= XBAD:
                continue
            ncom += 1
            sum += avg[kn]
            sumn += dnew[kn]
        if ncom < parameters.box_min_overlap:
            continue

        biask = float(sum-sumn)/ncom
        update_bias()

        # Update period of valid data, averages and weights
        for n in range(nf1, nl1):
            kn = k+12*n     # CSE for array index
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
def tavg(data, nyrs, base, limit, nr, id, deflt=0.0, XBAD=9999):
    """:meth:`tavg` computes the time averages for each calendar
    month over the base period (year *base* to *limit*) and
    saves them in *bias* (a fresh array that is returned).
    
    In case of no data, the average is set to
    *deflt* if nr=0 or computed over the whole period if nr>0.

    *id* is an arbitrary printable value used for identification in
    diagnostic output (for example, the cell number).
    """

    bias = [0.0]*12
    missed = 12
    len = 12*[0]    # Warning: shadows builtin "len"
    for k in range(12):
        bias[k] = deflt
        sum = 0.0
        m = 0
        for n in range(base, limit):
            kn = k+12*n     # CSE for array index
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
    for k in range(12):
        if len[k] > 0:
            continue
        print "No data in base period - MONTH,NR,ID", k, nr, id
        sum = 0.0
        m = 0
        for n in range(nyrs):
            kn = k+12*n     # CSE for array index
            if data[kn] >= XBAD:
                continue
            m += 1
            sum += data[kn]
        if m == 0:
            continue
        bias[k] = float(sum)/float(m)
    return bias


def zonav(boxed_data):
    """ 
    Perform Zonal Averaging.

    The input *boxed_data* is an iterator of boxed time series.
    The data in the boxes are combined to produce averages over
    various latitudinal zones.  Returns an iterator of
    (averages, weights, title) tuples, one per zone.

    14 belts are produced.  The first 8 are the basic belts that are used
    for the equal area grid, the remaining 6 are combinations:

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

    # Number of special zones.
    NZS = 3
    # Number of required common data in order to combine.

    (info, titlei) = boxed_data.next()
    kq = info[1]
    ml = info[3]
    nyrsin = info[3]/12
    # The rest of the code is structured so that MONM and IYRBEG could
    # be different from these values (for example, to make the output
    # file smaller, I think).
    iyrbeg = info[5]
    monm = ml
    # One more than the last year with data
    yearlimit = monm/12 + iyrbeg
    # Index of first month to be output
    mfout = (iyrbeg - info[5])*12

    # *reference* is a pair of (base,limit) that specifies the reference
    # period (base period).
    reference = (parameters.box_reference_first_year-iyrbeg,
                 parameters.box_reference_last_year-iyrbeg+1)
    # Typically 9999
    XBAD = info[6]
    last = info[6]
    trace = info[7]

    # Write output file header
    outinfo = list(info)
    outinfo[3] = monm
    outinfo[5] = iyrbeg
    yield (outinfo, titlei)

    ibm,kzone = zones()

    JBM = len(ibm)

    # :init:lenz:
    lenz = [None]*JBM

    wt = [None] * JBM
    avg = [None] * JBM
    for jb in range(JBM):
        ar = [None] * ibm[jb]
        wtr = [None] * ibm[jb]
        lenr = [None] * ibm[jb]
        for ib in range(ibm[jb]):
            (ar[ib], wtr[ib], lenr[ib], _) = boxed_data.next()
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
        for n in range(1,ibm[jb]):
            nr = IORD[n]
            if lenr[n] == 0:
                break
            # :todo: could the id parameter (the last one) be something
            # more useful, like "Zone 4"? Is it only used for
            # diagnostics?
            combine(avg[jb], bias, wt[jb], ar[nr], 0,nyrsin, wtr[nr], wtm, jb)
        bias = tavg(avg[jb], nyrsin, reference[0], reference[1], True, "Belt %d" % jb)
        lenz[jb] = 0
        m = 0
        for iy in range(nyrsin):
            for k in range(12):
                m = iy*12 + k
                if avg[jb][m] == XBAD:
                    continue
                avg[jb][m] -= bias[k]
                lenz[jb] += 1

        yield (avg[jb], wt[jb])

    # *lenz* contains the lemgths of each zone 0 to 7 (the number of
    # valid months in each zone).
    lenz, iord = sort_perm(lenz)
    for jz in range(NZS+3):
        if lenz[0] == 0:
            raise Error('**** NO DATA FOR ZONE %d' % JBM+jz)
        # Find the longest basic belt that is in the special zone.
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
        wtg = list(wt[jb])
        avgg = list(avg[jb])
        # Add in the remaining latitude belts JB with JB in kzone[jz].
        for j in range(j1+1,JBM):
            jb = iord[j]
            if jb not in kzone[jz]:
                continue
            if lenz[j] == 0:
                # Not convinced this behavious is either correct or
                # worth preserving.
                break
            combine(avgg, bias, wtg, avg[jb], 0,nyrsin, wt[jb], wtm, "Zone %d" % (JBM+jz))
        else:
            bias = tavg(avgg, nyrsin, reference[0], reference[1], True, "Zone %d" % (JBM+jz))
            m = 0
            for iy in range(nyrsin):
                for k in range(12):
                    if avgg[m] != XBAD:
                        avgg[m] -= bias[k]
                    m += 1
        yield(avgg, wtg)

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
    """Return the parameters of the 14 zones (8 basic belts and 6
    additional).  A pair of (*ibm*,*kzone*) is returned.  `ibm[b]`
    gives the number of boxes (regions) in belt *b* for `b in
    range(8)`.  *kzone* defines how the 6 combined zones are made from
    the basic belts.  `b in kzone[k]` is true when basic belt *b* is
    in special zone *k* (*b* is in range(8), *k* is in range(6)).
    """

    # Number of boxes (regions) in each band.
    ibm = [4,8,12,16,16,12,8,4]

    N = set(range(4)) # Northern hemisphere
    G = set(range(8)) # Global
    S = G - N         # Southern hemisphere
    T = set([3,4])    # Tropics
    kzone = [N-T, T, S-T, N, S, G]

    return ibm, kzone


def annzon(zoned_averages, alternate={'global':2, 'hemi':True}):
    """Compute annual zoned anomalies. *zoned_averages* is an iterator
    of zoned averages produced by `zonav`.

    The *alternate* argument controls whether alternate algorithms are
    used to compute the global and hemispheric means.
    alternate['global'] is 1 or 2, to select 1 of 2 different
    alternate computations, or false to not compute an alternative;
    alternate['hemi'] is true to compute an alternative, false
    otherwise.
    """

    jzm = 14
    monmin = 6

    (info, title) = zoned_averages.next()
    kq = info[1]
    iyrbeg = info[5]
    monm = info[3]
    iyrs = monm // 12
    # Allocate the 2- and 3- dimensional arrays.
    # The *data* and *wt* arrays are properly 3 dimensional
    # ([zone][year][month]), but the inner frames are only allocated
    # when we read the data, see :read:zonal below.
    data = [ None for _ in range(jzm)]
    wt =   [ None for _ in range(jzm)]
    ann =  [ [None]*iyrs for _ in range(jzm)]
    annw = [ [None]*iyrs for _ in range(jzm)]
    # Here we use the Python convention, *iyrend* is one past the highest
    # year used.
    iyrend = info[3] // 12 + iyrbeg
    XBAD = float(info[6])

    # Collect JZM zonal means.
    for jz in range(jzm):
        (tdata, twt) = zoned_averages.next()
        # Regroup the *data* and *wt* series so that they come in blocks of 12.
        # Uses essentially the same trick as the `grouper()` recipe in
        # http://docs.python.org/library/itertools.html#recipes
        data[jz] = zip(*[iter(tdata)]*12)
        wt[jz] = zip(*[iter(twt)]*12)

    # Find (compute) the annual means.
    for jz in range(jzm):
        for iy in range(iyrs):
            ann[jz][iy] = XBAD
            annw[jz][iy] = 0.
            anniy = 0.
            annwiy = 0.
            mon = 0
            for m in range(12):
                if data[jz][iy][m] == XBAD:
                    continue
                mon += 1
                anniy += data[jz][iy][m]
                annwiy += wt[jz][iy][m]
            if mon >= monmin:
                ann[jz][iy] = float(anniy)/mon
            annw[jz][iy] = annwiy/12.

    # Alternate global mean.
    if alternate['global']:
        glb = alternate['global']
        assert glb in (1,2)
        # Pick which "four" zones to use.
        # (subtracting 1 from each zone to convert to Python convention)
        if glb == 1:
            zone = [8, 9, 9, 10]
        else:
            zone = [8, 3, 4, 10]
        wtsp = [3.,2.,2.,3.]
        for iy in range(iyrs):
            glob = 0.
            ann[-1][iy] = XBAD
            for z,w in zip(zone, wtsp):
                if ann[z][iy] == XBAD:
                    # Note: Rather ugly use of "for...else" to emulate GOTO.
                    break
                glob += ann[z][iy]*w
            else:
                ann[-1][iy] = .1 * glob
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
    # Alternate hemispheric means.
    if alternate['hemi']:
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
    return (info, data, wt, ann, annw, monmin, title)

def step5(data):
    """Step 5 of GISTEMP.

    This step takes input provided by steps 3 and 4 (zipped together).

    :Param data:
        These are the land and ocean sub-box data, zipped together.
        They need to support the protocol defined by
        `code.giss_data.SubboxSetProtocol`.

    """
    boxed = SBBXtoBX(data)
    boxed = giss_io.step5_bx_output(boxed)
    zoned_averages = zonav(boxed)
    return annzon(zoned_averages)
