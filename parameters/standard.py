#! /usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/parameters.py $
# $Rev: 618 $
#
# parameters/standard.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

"""Parameters controlling the standard GISTEMP algorithm.

Various parameters controlling each phase of the algorithm are
collected and documented here.  They appear here in approximately the
order in which they are used in the algorithm.

Parameters controlling ccc-gistemp extensions to the standard GISTEMP
algorithm, or obsolete features of GISTEMP, are in other parameter
files.
"""
__docformat__ = "restructuredtext"

station_drop_minimum_months = 20
"""A station record must have at least one month of the year with at
least this many valid data values, otherwise it is dropped immediately
prior to the peri-urban adjustment step."""

rural_designator = "global_light <= 10"
"""Describes the test used to determine whether a station is rural or
not, in terms of the station metadata fields.  Relevant fields are:
'global_light' (global satellite nighttime radiance value); 'popcls'
(GHCN population class flag; the value 'R' stands for rural);
'us_light' (class derived from satellite nighttime radiance covering
the US and some neighbouring stations), 'berkeley' (a field of unknown
provenance which seems to be related to the Berkeley Earth Surface
Temperature project).

The value of this parameter may be a comma separated sequence.  Each
member in that sequence can either be a metadata field name, or a
numeric comparison on a metadata field name (e.g. "global_light <= 10",
the default).  If a field name appears on its own, the meaning is
field-dependent.

The fields are consulted in the order specified until one is found
that is not blank, and that obeys the condition (the only field which
is likely to be blank is 'us_light': this sequential feature is
required to emulate a previous version of GISTEMP).

Previous versions of GISTEMP can be "emulated" as follows:
"popcls" GISTEMP 1999 to 2001
"us_light, popcls" GISTEMP 2001 to 2010
"global_light <= 10" GISTEMP 2010 onwards
"global_light <= 0" GISTEMP 2011 passing 2 as second arg to do_comb_step2.sh 
"berkeley <= 0" GISTEMP 2011 passing 3 as second arg to do_comb_step2.sh 
"""

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

gridding_min_overlap = 20
"""When combining station records to give a grid record, do not
combine a candidate station record if it has fewer than this number of
years of overlap with the combined grid record."""

gridding_radius = 1200.0
"""The radius in kilometres used to find and weight station records to
give a grid record."""

gridding_reference_period = (1951, 1980)
"""When gridding, temperature series are turned into anomaly series by
subtracting monthly means computed over a reference period.  This is
the first and last years of that reference period."""

sea_surface_cutoff_temp = -1.77
"""When incorporating monthly sea-surface datasets, treat any
temperature colder than this as missing data."""

subbox_min_valid = 240
"""When combining the sub-boxes into boxes, do not use any sub-box
record, either land or ocean, which has fewer than this number of
valid data."""

subbox_land_range = 100
"""If a subbox has both land data and ocean data, but the distance
from the subbox centre to the nearest station used in its record is
less than this, the land data is used in preference to the ocean data
when calculating the box series. Note: the distance used is actually a
great-circle chord length."""

subbox_reference_period = (1961, 1990)
"""When combining subbox records into box records, temperature series
are turned into anomaly series by subtracting monthly means computed
over a reference period.  This is the first and last years of that
reference period."""

box_min_overlap = 20
"""When combining subbox records to make box records, do not combine a
calendar month from a candidate subbox record if it has fewer than
this number of years of overlap with the same calendar month in the
combined box record.  Also used when combining boxes into zones."""

box_reference_period = (1951, 1980)
"""When combining box records into zone records, temperature series
are turned into anomaly series by subtracting monthly means computed
over a reference period.  This is the first and last years of that
reference period."""

zone_annual_min_months = 6
"""When computing zone annual means, require at least this many valid
month data."""
