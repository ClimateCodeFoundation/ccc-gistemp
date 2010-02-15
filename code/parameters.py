#! /usr/bin/env python
# $URL$
# $Rev$
# 
# parameters.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

#!/usr/bin/env python
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
