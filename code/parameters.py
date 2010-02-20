#! /usr/bin/env python
# $URL$
# $Rev$
# 
# parameters.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

"""Parameters controlling the GISTEMP algorithm.

Various numeric parameters controlling each phase of the algorithm are
collected and documented here.
"""
__docformat__ = "restructuredtext"

# TODO: Go through the whole algorithm collecting parameters and
# moving them here.

USHCN_offset_start_year = 1980
"""The first year of the period considered when calculating the
offsets between GHCN and USHCN records, to apply to a USHCN
temperature series.
"""

USHCN_offset_max_months = 10
"""The maximum number of records used when computing a monthly offset
between GHCN and USHCN, to apply to a USHCN temperature series.  The
algorithm starts in the current year and works back to
*USHCN_offset_start_year*, computing an offset for each month using up
to *USHCN_offset_max_months* records."""

station_combine_min_overlap = 4
"""The minimum number of years of overlap, between a combined record
and a candidate record, to allow the candidate to be added into the
combined record."""

station_combine_bucket_radius = 10
"""Used when deciding whether to add a non-overlapping station record
to a combined record.  This number of years is considered, either side
of the centre of the potential new combined record.  If there are
enough valid years in each record (see *station_combine_min_mid_years*),
and the difference between the average anomaly of the combined record
and the average anomaly of the candidate record, over those years, is
less than the standard deviation of the combined record, then the
record is combined."""

station_combine_min_mid_years = 5
"""The minimum number of acceptable year records in the central
"bucket", for both the combined record and the candidate record, when
combining non-overlapping station records."""

station_drop_minimum_months = 20
"""A combined station record must have at least one month of the year
with at least this many valid data values, otherwise it is dropped
immediately prior to the peri-urban adjustment step."""

urban_adjustment_min_years = 20
"""When trying to calculate an urban station adjustment, at least this
many years have to have sufficient rural stations (if there are
not enough qualifying years, we may try again at a larger radius)."""

urban_adjustment_proportion_good = 2.0 / 3.0
"""When trying to calculate an urban station adjustment, at least this
proportion of the years to which the fit applies have to have
sufficient rural stations (if there are insufficient stations, we may
try again at a larger radius)."""

urban_adjustment_min_rural_stations = 3
"""When trying to calculate an urban station adjustment, a year
without at least this number of valid readings from rural stations is
not used to calculate the fit."""

urban_adjustment_min_leg = 5
"""When finding a two-part adjustment, only consider knee years which
have at least this many data points (note: not years) on each side."""

urban_adjustment_short_leg = 7
"""When a two-part adjustment has been identified, if either leg is
shorter than this number of years, a one-part adjustment is applied
instead."""

urban_adjustment_steep_leg = 0.1
"""When a two-part adjustment has been identified, if the gradient of
either leg is steeper than this (in absolute degrees Celsius per
year), or if the difference between the leg gradients is greater than
this, a one-part adjustment is applied instead."""

urban_adjustment_leg_difference = 0.05
"""When a two-part adjustment has been identified, if the difference
in gradient between the two legs is greater than this (in absolute
degrees Celsius per year), it is counted separately for statistical
purposes."""

urban_adjustment_reverse_gradient = 0.02
"""When a two-part adjustment has been identified, if the two
gradients have opposite sign, and both gradients are steeper than this
(in absolute degrees Celsius per year), a one-part adjustment is
applied instead."""

urban_adjustment_full_radius = 1000.0
"""Range in kilometres within which a rural station will be considered
for adjusting an urban station record.  Half of this radius will be
attempted first."""

rural_station_min_overlap = 20
"""When combining rural station annual anomaly records to calculate
urban adjustment parameters, do not combine a candidate rural record
if it has fewer than this number years of overlap."""

