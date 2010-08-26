#!/usr/bin/env python
# $URL$
# $Rev$
#
# series.py
#
# Nick Barnes, Ravenbrook Limited, 2010-03-08

import itertools
from giss_data import valid, invalid, MISSING

"""
Shared series-processing code in the GISTEMP algorithm.
"""


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

    new_weight = ensure_array(weight, new_weight)

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

def ensure_array(exemplar, item):
    """Coerces *item* to be an array (linear sequence); if *item* is
    already an array it is returned unchanged.  Otherwise, an array of
    the same length as exemplar is created which contains *item* at
    every index.  The fresh array is returned.
    """

    try:
        item[0]
    except TypeError:
        return (item,)*len(exemplar)

def anomalize(data, reference_period=None, base_year=-9999):
    """Turn the series *data* into anomalies, based on monthly
    averages over the *reference_period*, for example (1951, 1980).
    *base_year* is the first year of the series.  If *reference_period*
    is None then the averages are computed over the whole series.
    Similarly, If any month has no data in the reference period,
    the average for that month is computed over the whole series.

    The *data* sequence is mutated.
    """

    means, anoms = monthly_anomalies(data, reference_period, base_year)
    # Each of the elements in *anoms* are the anomalies for one of the
    # months of the year (for example, January).  We need to splice each
    # month back into a single linear series.
    for m in range(12):
        data[m::12] = anoms[m]

def valid_mean(seq, min=1):
    """Takes a sequence, *seq*, and computes the mean of the valid
    items (using the valid() function).  If there are fewer than *min*
    valid items, the mean is MISSING."""

    count = 0
    sum = 0.0
    for x in seq:
        if valid(x):
            sum += x
            count += 1
    if count >= min:
        return sum/float(count)
    else:
        return MISSING

def monthly_anomalies(data, reference_period=None, base_year=-9999):
    """Calculate monthly anomalies, by subtracting from every datum
    the mean for its month.  A pair of (monthly_mean, monthly_anom) is
    returned.  *monthly_mean* is a 12-long sequence giving the mean for
    each of the 12 months; *monthly_anom* is a 12-long sequence giving
    the anomalized series for each of the 12 months.

    If *reference_period* is supplied then it should be a pair (*first*,
    *last) and the mean for a month is taken over the period (an example
    would be reference_period=(1951,1980)).  *base_year* specifies the
    first year of the data.
    
    The input data is a flat sequence, one datum per month.
    Effectively the data changes shape as it passes through this
    function.
    """

    years = len(data) // 12
    if reference_period:
        base = reference_period[0] - base_year
        limit = reference_period[1] - base_year + 1
    else:
        # Setting base, limit to (0,0) is a bit of a hack, but it
        # does work.
        base = 0
        limit = 0
    monthly_mean = []
    monthly_anom = []
    for m in range(12):
        row = data[m::12]
        mean = valid_mean(row[base:limit])
        if invalid(mean):
            # Fall back to using entire period
            mean = valid_mean(row)
        monthly_mean.append(mean)
        if valid(mean):
            def asanom(datum):
                """Convert a single datum to anomaly."""
                if valid(datum):
                    return datum - mean
                return MISSING
            monthly_anom.append(map(asanom, row))
        else:
            monthly_anom.append([MISSING]*years)
    return monthly_mean, monthly_anom
