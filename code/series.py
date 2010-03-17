#!/usr/bin/env python
# $URL$
# $Rev$
#
# series.py
#
# Nick Barnes, Ravenbrook Limited, 2010-03-08

import itertools
from giss_data import valid, invalid

"""
Shared series-processing code in the GISTEMP algorithm.
"""

def container(item):
    """Takes either a container or a non-container *item*,
    returns either the same *item* or a container which
    contains that item at every index."""
    try:
        item[0]
    except TypeError:
        constant = item
        class constant_list:
            def __getitem__(self, i):
                return constant
        item = constant_list()
    return item


def combine(average, weight, new, new_weight,
            first_year, last_year, min_overlap):
    """Run the GISTEMP combining algorithm.  This combines the data
    in the *new* array into the *average* array.  *new* has weight
    *new_weight*, *average* has weights in the *weight* array.

    Only data for years in *range(first_year, last_year)* are
    considered and combined.

    *new_weight* can be either a constant or an array of weights for
     each datum in *new*.

    The number of month records combined is returned.

    Each month of the year is considered separately.  For the set of
    times where both *average* and *new* have data the mean difference
    (a bias) is computed.  If there are fewer than *min_overlap* years
    in common the data (for that month of the year) are not combined.
    The bias is subtracted from the *new* record and it is point-wise
    combined into *average* according to the weight *new_weight* and
    the existing weights for *average*.
    """

    new_weight = container(new_weight)

    months_combined = 0
    for m in range(12):
        sum_new = 0.0    # Sum of data in new
        sum = 0.0        # Sum of data in average
        count = 0    # Number of years where both new and average are valid
        for a,n in itertools.izip(average[first_year*12+m: last_year*12: 12],
                                  new[first_year*12+m: last_year*12: 12]):
            if invalid(a) or invalid(n):
                continue
            count += 1
            sum += a
            sum_new += n
        if count < min_overlap:
            continue
        bias = (sum-sum_new)/count

        # Update period of valid data, averages and weights
        for i in range(first_year*12+m, last_year*12, 12):
            if invalid(new[i]):
                continue
            new_month_weight = weight[i] + new_weight[i]
            average[i] = (weight[i]*average[i]
                          + new_weight[i]*(new[i]+bias))/new_month_weight
            weight[i] = new_month_weight
            months_combined += 1
    return months_combined

def anomalize(data, reference_period=None, base_year=-9999):
    """Turn the series *data* into anomalies, based on monthly
    averages over the *reference_period*, for example (1951, 1980).
    *base_year* is the first year of the series.  If *reference_period*
    is None then the averages are computed over the whole series.
    Similarly, If any month has no data in the reference period,
    the average for that month is computed over the whole series.

    The *data* sequence is mutated.
    """
    if reference_period:
        base = reference_period[0] - base_year
        limit = reference_period[1] - base_year + 1
    else:
        # Setting base, limit to (0,0) is a bit of a hack, but it
        # does work.
        base = 0
        limit = 0
    for m in range(12):
        sum = 0.0
        count = 0
        for datum in data[m+12*base:m+12*limit:12]:
            if invalid(datum):
                continue
            count += 1
            sum += datum
        if count == 0:
            for datum in data[m::12]:
                if invalid(datum):
                    continue
                count += 1
                sum += datum
        if count == 0:
            average = 0.0
        else:
            average = sum/count
        index = m
        while index < len(data):
            if valid(data[index]):
                data[index] -= average
            index += 12
