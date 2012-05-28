#! /usr/bin/env python
# $URL$
# $Rev$
#
# step1.py
#
# Nick Barnes and David Jones.
# Copyright (C) Ravenbrook Limited, 2008-2010.
# Copyright (C) Climate Code Fonudation, 2010-2012.
#
# BSD license, see license.txt

"""
Python code reproducing the Step 1 part of the GISTEMP algorithm.  In
this step, certain records (or parts) are adjusted or dropped, under
the control of a configuration file.

Requires this file from GISTEMP STEP0/input_files/:

Ts.strange.v3.list.IN_full
"""

# Clear Climate Code
import read_config
from giss_data import MISSING, BASE_YEAR

def drop_strange(data):
    """Drops station records, or parts of records, under control of
    the file 'input/Ts.strange.v3.list.IN_full' file.
    """

    changes_dict = read_config.get_changes_dict()
    for record in data:
        changes = changes_dict.get(record.uid, [])
        series = record.series
        begin = record.first_year
        # :todo: Use record.last_year
        end = begin + (len(series)//12) - 1
        for (kind, year, x) in changes:
            if kind == 'years':
                # omit all the data from year1 to year2, inclusive
                year1 = year
                year2 = x
                if year1 <= begin and year2 >= end:
                    # Drop this whole record.  Note: avoids "else:"
                    # clause at end of "for" loop.
                    break

                # Clamp range of deleted years to the range of the
                # series.
                year1 = max(year1, begin)
                year2 = min(year2, end)
                if year2 < year1:
                    # Happens when deleted range is entirely outside the
                    # range of the series.  In which case, pass record
                    # unchanged.
                    continue
                # Invalidate the data.
                nmonths = (year2 + 1 - year1) * 12
                series[(year1-begin)*12:(year2+1-begin)*12] = [
                        MISSING] * nmonths

            else: # remove a single month
                # It can happen that the datum-to-be-removed is outside
                # the date range for this record (if we are using new
                # config files, and old data).  So we check.
                if begin <= year <= end:
                    series[(year-begin)*12 + x-1] = MISSING
        else:
            record.set_series(begin * 12 + 1, series)
            yield record

def step1(records):
    """An iterator for step 1.  Produces a stream of
    `giss_data.Series` instances.

    :Param records:
        An iterable source of `giss_data.Series` instances (which it
        will assume are station records).
    """
    without_strange = drop_strange(records)
    for record in without_strange:
        assert record.first_year == BASE_YEAR
        yield record
