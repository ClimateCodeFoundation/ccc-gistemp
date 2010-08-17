#!/usr/bin/env python
# $URL$
# $Rev$
#
# step5.py
#
# David Jones, Ravenbrook Limited, 2009-10-27

"""
Step 5 of the GISTEMP algorithm.

In Step 5: 8000 subboxes are combined into 80 boxes, and ocean data is
combined with land data; boxes are combined into latitudinal zones
(including hemispheric and global zones); annual and seasonal anomalies
are computed from monthly anomalies.
"""

# Clear Climate Code
import eqarea
import giss_data
import parameters
import series
from tool import giss_io
from giss_data import valid, invalid, MISSING

# http://www.python.org/doc/2.3.5/lib/itertools-functions.html
import itertools
# http://www.python.org/doc/2.3.5/lib/module-math.html
import math
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os

def blank_ocean_data(data):
    """Augment a land-only data series with blank ocean data."""
    for land_box in data:
        ocean_series = [MISSING] * len(land_box.series)
        ocean_box = giss_data.Series(series=ocean_series,
            box=land_box.box,
            stations=0, station_months=0,
            d=MISSING)
        yield land_box, ocean_box


def select(data):
    """Taking a stream of (land,ocean) record pairs, decide how to
    weight the two series in subsequent combination and yield a stream
    of (weight,land,ocean) triples.  *weight* will be 1 when the land
    record is to be used, and 0 if the ocean record is to be used.
    """

    # First item from iterator is normally a pair of metadata objects,
    # one for land, one for ocean.  If we are piping step3 straight into
    # step5 then it is not a pair.  In that case we synthesize missing
    # ocean data.
    meta = data.next()
    try:
        land_meta, ocean_meta = meta
    except (TypeError, ValueError):
        # Use the land meta object for both land and ocean data
        land_meta,ocean_meta = meta, meta
        print "No ocean data; using land data only"
        data = blank_ocean_data(data)

    yield land_meta, ocean_meta

    for land,ocean in data:
        if (ocean.good_count < parameters.subbox_min_valid
            or land.d < parameters.subbox_land_range):
            landmask = 1.0
        else:
            landmask = 0.0

        yield landmask, land, ocean


def SBBXtoBX(data):
    """Simultaneously combine the land series and the ocean series and
    combine subboxes into boxes.  *data* should be an iterator of
    (weight, land, ocean) triples (a land and ocean series for each
    subbox).  When *weight* is 1 the land series is selected; when
    *weight* is 0 the ocean series is selected.  Currently no
    intermediate weights are supported.

    Returns an iterator of box data.
    """

    meta = data.next()
    land_meta, ocean_meta = meta

    # TODO: Formalise use of only monthlies, see step 3.
    assert land_meta.mavg == 6
    combined_year_beg = min(land_meta.yrbeg, ocean_meta.yrbeg)
    # Index into the combined array of the first year of the land data.
    land_offset = 12*(land_meta.yrbeg-combined_year_beg)
    # As land_offset but for ocean data.
    ocean_offset = 12*(ocean_meta.yrbeg-combined_year_beg)
    combined_n_months = max(land_meta.monm + land_offset,
                            land_meta.monm + ocean_offset)

    info = [land_meta.mo1, land_meta.kq, land_meta.mavg, land_meta.monm,
            land_meta.monm4, combined_year_beg, land_meta.missing_flag,
            land_meta.precipitation_flag]

    info[4] = 2 * land_meta.monm + 5
    yield(info, land_meta.title)

    # Number of subboxes within each box.
    nsubbox = 100

    for box_number,box in enumerate(eqarea.grid()):
        # Averages for the land and ocean (one series per subbox)...
        avg = []
        wgtc = []
        # Eat the records from land and ocean 100 (nsubbox) at a time.
        # In other words, all 100 subboxes for the box.
        landweight,landsub,oceansub = zip(*itertools.islice(data, nsubbox))
        # :todo: combine below zip with above zip?
        for i, t, l, o in zip(range(nsubbox), landweight, landsub, oceansub):
            # Simple version only selects either land or ocean
            assert t in (0,1)
            a = [MISSING]*combined_n_months
            if t:
                # use land series for this subbox
                a[land_offset:land_offset+len(l.series)] = l.series
                wgtc.append(l.good_count)
            else:
                # use ocean series for this subbox
                a[ocean_offset:ocean_offset+len(o.series)] = o.series
                wgtc.append(o.good_count)
            avg.append(a)

        # GISTEMP sort.
        # We want to end up with IORDR, the permutation array that
        # represents the sorted order.  IORDR[0] is the index (into the
        # *wgtc* array) of the longest record, IORDR[1] the index of the
        # next longest record, and so on.
        # :todo: should probably import from a purpose built module.
        from step3 import sort
        IORDR = range(nsubbox)
        sort(IORDR, lambda x,y: wgtc[y] - wgtc[x])

        # From here to the "for" loop over the cells (below) we are
        # initialising data for the loop.  Primarily the AVGR and WTR
        # arrays.
        nc = IORDR[0]

        # Weights for the box's record.
        wtr = [a != MISSING for a in avg[nc]]
        # Box record
        avgr = avg[nc][:]

        # Loop over the remaining cells.
        for nc in IORDR[1:]:
            if wgtc[nc] >= parameters.subbox_min_valid:
                series.combine(avgr, wtr, avg[nc], 1, 0,
                           combined_n_months/12, parameters.box_min_overlap)

        series.anomalize(avgr, parameters.subbox_reference_period,
                         combined_year_beg)
        ngood = sum(valid(a) for a in avgr)
        yield (avgr, wtr, ngood, box)
    # We've now consumed all 8000 input boxes and yielded 80 boxes.  We
    # need to tickle the input to check that it is exhausted and to
    # cause it to run the final tail of its generator.
    # We expect the call to .next() to raise StopIteration, which is
    # just what we want.
    data.next()
    # Ordinarily we never reach here.
    assert 0, "Too many input records"


def zonav(boxed_data):
    """Zonal Averaging.

    The input *boxed_data* is an iterator of boxed time series.
    The data in the boxes are combined to produce averages over
    various latitudinal zones.  Returns an iterator of
    (averages, weights, title) tuples, one per zone.

    14 zones are produced.  The first 8 are the basic belts that are used
    for the equal area grid, the remaining 6 are combinations:

      0 64N - 90N               \
      1 44N - 64N (asin 0.9)     -  8 24N - 90 N  (0 + 1 + 2)
      2 24N - 44N (asin 0.7)    /
      3 Equ - 24N (asin 0.4)    \_  9 24S - 24 N  (3 + 4)
      4 24S - Equ               /
      5 44S - 24S               \
      6 64S - 44S                - 10 90S - 24 S  (5 + 6 + 7)
      7 90S - 64S               /

     11 northern hemisphere (0 + 1 + 2 + 3)
     12 southern hemisphere (4 + 5 + 6 + 7)
     13 global (all belts 0 to 7)
    """

    (info, titlei) = boxed_data.next()
    iyrbeg = info[5]
    monm = info[3]
    nyrsin = monm/12
    # One more than the last year with data
    yearlimit = nyrsin + iyrbeg

    yield (info, titlei)

    boxes_in_band,band_in_zone = zones()

    bands = len(boxes_in_band)

    lenz = [None] * bands
    wt = [None] * bands
    avg = [None] * bands
    # For each band, combine all the boxes in that band to create a band
    # record.
    for band in range(bands):
        # The temperature (anomaly) series for each of the boxes in this
        # band.
        box_series = [None] * boxes_in_band[band]
        # The weight series for each of the boxes in this band.
        box_weights = [None] * boxes_in_band[band]
        # "length" is the number of months (with valid data) in the box
        # series.  For each box in this band.
        box_length = [None] * boxes_in_band[band]
        for box in range(boxes_in_band[band]):
            # The last element in the tuple is the boundaries of the
            # box.  We ignore it.
            box_series[box], box_weights[box], box_length[box], _ = (
              boxed_data.next())
        # total number of valid data in band's boxes
        total_length = sum(box_length)
        if total_length == 0:
            wt[band] = [0.0]*monm
            avg[band] = [MISSING]*monm
        else:
            box_length,IORD = sort_perm(box_length)
            nr = IORD[0]
            # Copy the longest box record into *wt* and *avg*.
            # Using list both performs a copy and converts into a mutable
            # list.
            wt[band] = list(box_weights[nr])
            avg[band] = list(box_series[nr])
            # And combine the remaining series.
            for n in range(1,boxes_in_band[band]):
                nr = IORD[n]
                if box_length[n] == 0:
                    # Nothing in this box, and since we sorted by length,
                    # all the remaining boxes will also be empty.  We can
                    # stop combining boxes.
                    break
                series.combine(avg[band], wt[band],
                  box_series[nr], box_weights[nr], 0, nyrsin,
                  parameters.box_min_overlap)
        series.anomalize(avg[band], parameters.box_reference_period, iyrbeg)
        lenz[band] = sum(valid(a) for a in avg[band])
        yield (avg[band], wt[band])

    # We expect to have consumed all the boxes (the first 8 bands form a
    # partition of the boxes).  We check that the boxed_data stream is
    # exhausted and contains no more boxes.
    try:
        boxed_data.next()
        assert 0, "Too many boxes found"
    except StopIteration:
        # We fully expect to get here.
        pass

    # *lenz* contains the lengths of each zone 0 to 7 (the number of
    # valid months in each zone).
    lenz, iord = sort_perm(lenz)
    for zone in range(len(band_in_zone)):
        if lenz[0] == 0:
            raise Error('**** NO DATA FOR ZONE %d' % bands+zone)
        # Find the longest band that is in the compound zone.
        for j1 in range(bands):
            if iord[j1] in band_in_zone[zone]:
                break
        else:
            # Should be an assertion really.
            raise Error('No band in compound zone %d.' % zone)
        band = iord[j1]
        wtg = list(wt[band])
        avgg = list(avg[band])
        # Add in the remaining bands, in length order.
        for j in range(j1+1,bands):
            band = iord[j]
            if band not in band_in_zone[zone]:
                continue
            series.combine(avgg, wtg, avg[band], wt[band], 0,nyrsin,
                           parameters.box_min_overlap)
        series.anomalize(avgg, parameters.box_reference_period, iyrbeg)
        yield(avgg, wtg)

def sort_perm(a):
    """The array *a* is sorted into descending order.  The fresh sorted
    array and the permutation array are returned as a pair (*sorted*,
    *indexes*).  The original *a* is not mutated.

    The *indexes* array is such that `a[indexes[x]] == sorted[x]`.
    """
    from step3 import sort
    z = zip(a, range(len(a)))
    sort(z, lambda x,y: y[0]-x[0])
    sorted, indexes = zip(*z)
    return sorted, indexes

def zones():
    """Return the parameters of the 14 zones (8 basic bands and 6
    additional compound zones).
    
    A pair of (*boxes_in_band*,*band_in_zone*) is returned.
    `boxes_in_band[b]` gives the number of boxes in band
    *b* for `b in range(8)`.  *band_in_zone* defines how the 6
    combined zones are made from the basic bands.  `b in
    band_in_zone[k]` is true when basic band *b* is in compound zone
    *z* (*b* is in range(8), *z* is in range(6)).

    Implicit (in this function and its callers) is that 8 basic bands
    form a decomposition of the 80 boxes.  All you need to know is the
    number of boxes in each band; simply take the next N boxes to make
    the next band.
    """

    # Number of boxes (regions) in each band.
    boxes_in_band = [4,8,12,16,16,12,8,4]

    N = set(range(4)) # Northern hemisphere
    G = set(range(8)) # Global
    S = G - N         # Southern hemisphere
    T = set([3,4])    # Tropics
    band_in_zone = [N-T, T, S-T, N, S, G]

    return boxes_in_band, band_in_zone


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

    zones = 14

    (info, title) = zoned_averages.next()
    iyrbeg = info[5]
    monm = info[3]
    iyrs = monm // 12
    iyrend = iyrs + iyrbeg

    # Allocate the 2- and 3- dimensional arrays.
    # The *data* and *wt* arrays are properly 3 dimensional
    # ([zone][year][month]), but the inner frames are only allocated
    # when we read the data, see :read:zonal below.
    data = [ None for _ in range(zones)]
    wt =   [ None for _ in range(zones)]
    ann =  [ [MISSING]*iyrs for _ in range(zones)]

    # Collect zonal means.
    for zone in range(zones):
        (tdata, twt) = zoned_averages.next()
        # Regroup the *data* and *wt* series so that they come in blocks of 12.
        # Uses essentially the same trick as the `grouper()` recipe in
        # http://docs.python.org/library/itertools.html#recipes
        data[zone] = zip(*[iter(tdata)]*12)
        wt[zone] = zip(*[iter(twt)]*12)

    # Find (compute) the annual means.
    for zone in range(zones):
        for iy in range(iyrs):
            anniy = 0.
            mon = 0
            for m in range(12):
                if data[zone][iy][m] == MISSING:
                    continue
                mon += 1
                anniy += data[zone][iy][m]
            if mon >= parameters.zone_annual_min_months:
                ann[zone][iy] = float(anniy)/mon

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
            ann[-1][iy] = MISSING
            for z,w in zip(zone, wtsp):
                if ann[z][iy] == MISSING:
                    # Note: Rather ugly use of "for...else" to emulate GOTO.
                    break
                glob += ann[z][iy]*w
            else:
                ann[-1][iy] = .1 * glob
        for iy in range(iyrs):
            data[-1][iy] = [MISSING]*12
            for m in range(12):
                glob = 0.
                for z,w in zip(zone, wtsp):
                    if data[z][iy][m] == MISSING:
                        break
                    glob += data[z][iy][m]*w
                else:
                    data[-1][iy][m] = .1 * glob

    # Alternate hemispheric means.
    if alternate['hemi']:
        # For the computations it will be useful to recall how the zones
        # are numbered.  There is a useful docstring at the beginning of
        # zonav().
        for ihem in range(2):
            for iy in range(iyrs):
                ann[ihem+11][iy] = MISSING
                if (ann[ihem+3][iy] != MISSING and
                  ann[2*ihem+8][iy] != MISSING):
                    ann[ihem+11][iy] = (0.4*ann[ihem+3][iy] +
                                        0.6*ann[2*ihem+8][iy])
            for iy in range(iyrs):
                data[ihem+11][iy] = [MISSING]*12
                for m in range(12):
                    if (data[ihem+3][iy][m] != MISSING and
                      data[2*ihem+8][iy][m] != MISSING):
                        data[ihem+11][iy][m] = (
                          0.4*data[ihem+3][iy][m] +
                          0.6*data[2*ihem+8][iy][m])
    return (info, data, wt, ann, parameters.zone_annual_min_months, title)

def step5(data):
    """Step 5 of GISTEMP.

    This step takes input provided by steps 3 and 4 (zipped together).

    :Param data:
        These are the land and ocean sub-box data, zipped together.
        They need to support the protocol defined by
        `code.giss_data.SubboxSetProtocol`.

    """
    boxed = SBBXtoBX(select(data))
    boxed = giss_io.step5_bx_output(boxed)
    zoned_averages = zonav(boxed)
    return annzon(zoned_averages)
