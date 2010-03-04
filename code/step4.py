#!/usr/bin/env python
# $URL$
# $Rev$

"""Incorporates recent sea-surface temperature records into the
sea-surface temperature anomaly boxed dataset.
"""
__docformat__ = "restructuredtext"


import sys
from tool import giss_io
import itertools

import giss_data


# Constants (from original fortran parameters)
IM = 360         # degrees of longitude
JM = 180         # degrees of latitude
IOFF = IM // 2   # half longitude
IYRBEG = 1880    # first year

def merge_ocean(ocean, sst, dates):
    """Adds the array *sst* of new monthly sea-surface temperature readings,
    which has data for the dates *dates*, to the boxed iterator *ocean*.
    Returns a new boxed iterator.
    """

    clim = giss_io.step4_load_clim()

    first_new_year = dates[0][0]
    last_new_year = dates[-1][0]
    last_new_month = dates[-1][1]

    reader = iter(ocean)
    meta = reader.next()
    meta.monm = 12 * (last_new_year - IYRBEG + 1)
    meta.monm4 = meta.monm + 8
    meta.title = (meta.title[:40] + " Had: 1880-11/1981, oi2: 12/1981-%2d/%04d" % (last_new_month, last_new_year))
    yield meta

    # Interpolate to Sergej's subbox grid
    for box in reader:
        box.pad_with_missing(meta.monm)

        js = (18001 + (box.lat_S + 9000) * JM)/18000
        jn = (17999 + (box.lat_N + 9000) * JM)/18000
        iw = (36001 + (box.lon_W + 18000) * IM)/36000 + IOFF
        ie = (35999 + (box.lon_E + 18000) * IM)/36000 + IOFF
        if ie > IM:
            iw = iw - IM
            ie = ie - IM
        assert 1 <= iw <= IM and 1 <= ie <= IM 

        for y, m in dates:
            mm = (y - first_new_year) * 12 + m
            month = (m - 1) % 12
            count = 0
            sum = 0.0
            for j in range(js -1, jn):
                for i in range(iw - 1, ie):
                    if sst[i][j][mm-1] < -1.77 or clim[i][j][month] == giss_data.XMISSING:
                        continue
                    count += 1
                    sum += sst[i][j][mm-1] - clim[i][j][month]

            index = (y - IYRBEG) * 12 + mm - 1
            box.set_value(index, giss_data.XMISSING)
            if count > 0:
                box.set_value(index, sum / count)

        box.trim()
        yield box


def step4(data):
    """Step 4 of GISTEMP processing.  This is a little unusual
    compared to the other steps.  The input data is a 3-tuple *(land,
    ocean, monthlies)*, *land* and *ocean* being iterables of boxed
    data sets and *monthlies* being data read from any sea-surface
    monthly datasets which are present (and more recent than the
    *ocean* data).
    """

    land, ocean, monthlies = data
    if monthlies is not None:
        sst, dates = monthlies
        ocean = merge_ocean(ocean, sst, dates)

    return itertools.izip(land, ocean)
