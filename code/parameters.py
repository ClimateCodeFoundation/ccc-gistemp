#! /usr/bin/env python
# $URL$
# $Rev$
#
# parameters.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

"""Parameters controlling the GISTEMP algorithm.

Various numeric parameters controlling each phase of the algorithm are
collected and documented here.  They appear here in approximately the
order in which they are used in the algorithm.
"""
__docformat__ = "restructuredtext"

data_sources = "ghcn ushcn hohenpeissenberg scar"
"""Data sources that are used for the analysis (space separated string).
'ghcn' is the Global Historical Climate Network (NOAA); 'ushcn' is
the US Historical Climate Network (NOAA); 'hohenpeissenberg' is
Wiljens data for Hohenpeissenberg; 'scar' is the READER data from
Scientific Committee on Antarctic Research.
"""

work_file_format = "v2"
"""The format of the intermediate files written to the 'work' directory:
'v2' for GHCN v2, 'v3' for GHCN v3.
"""

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

retain_contiguous_US = True
"""Whether to retain GHCN records that are in the contiguous US.
Ordinarily, USHCN replaces GHCN records, but: some stations have
duplicates; some GHCN stations have no USHCN counterpart.  This
parameter affects those records.
"""

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

rural_designator = "global_light"
"""Specifies which station metadata fields are used to determine whether
a station is rural or not.  Only certain fields may be specified:
'global_light' (global satellite nighttime radiance value); 'popcls'
(GHCN population class); 'us_light' (class derived from satellite nighttime
radiance covering the US and some neighbouring stations).

The value of this field may be a comma separated sequence, in which case
the fields are consulted in the order specified until one is found that
is not blank (the only field for which this is useful is the 'us_light'
field, and the only intended use of this feature is to emulate a
previous version of GISTEMP).

Previous versions of GISTEMP can be "emulated" as follows:
"popcls" GISTEMP 1999 to 2001
"us_light,popcls" GISTEMP 2001 to 2010
"global_light" GISTEMP 2010 onwards
"""

max_rural_brightness = 10
"""The maximum brightness value for a station to be considered rural,
when 'global_light' is the field used in determing whether a station
is rural or not (see *rural_designator*)."""

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
