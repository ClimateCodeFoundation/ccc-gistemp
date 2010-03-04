#!/usr/bin/env python
# $URL$
# $Rev$

"""Incorporates recent sea-surface temperature records into the
sea-surface temperature anomaly boxed dataset.
"""
__docformat__ = "restructuredtext"


import sys
import itertools

import giss_data
import parameters
from tool import giss_io


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

    # Average into Sergej's subbox grid
    for box in reader:
        box.pad_with_missing(meta.monm)

        # identify all the degree boxes which are included in this subbox
        js = int(box.lat_S + 90.01)
        jn = int(box.lat_N + 89.99)
        iw = int(box.lon_W + 360.01)
        ie = int(box.lon_E + 359.99)
        if ie >= 360:
            iw = iw - 360
            ie = ie - 360

        for y, m in dates:
            mm = (y - first_new_year) * 12 + m
            month = (m - 1) % 12
            count = 0
            sum = 0.0
            for j in range(js, jn+1):
                for i in range(iw, ie+1):
                    if sst[i][j][mm-1] < parameters.sea_surface_cutoff_temp or clim[i][j][month] == giss_data.XMISSING:
                        continue
                    count += 1
                    sum += sst[i][j][mm-1] - clim[i][j][month]

            index = (y - IYRBEG) * 12 + m - 1
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
