#! /usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/parameters.py $
# $Rev: 618 $
#
# parameters/obsolete.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

"""Parameters controlling obsolete features of the GISTEMP algorithm.

Parameters controlling the standard GISTEMP algorithm, or ccc-gistemp
extensions to the GISTEMP algorithm, are in other parameter files.
"""
__docformat__ = "restructuredtext"

USHCN_convert_id = True
"""Whether to convert USHCN station identifiers to GHCN station
identifiers (using tables supplied by GISS), when using USHCN as a
data source.  In the usual analysis prior to using GHCN v3, this had
to be True.  Irrelevant in the usual analysis using GHCN v3, as USHCN
is not used."""

USHCN_meta = ''
"""Specifies what file to use for USHCN metadata (latitude,
longitude), when using USHCN as a data source.  If empty then the GHCN
metadata file is used (and this will only work when USHCN_convert_id
is True).  Specify "ushcn-v2-stations.txt" and place that file in the
input/ directory to use the USHCN metadata.  Irrelevant in the usual
analysis using GHCN v3, as USHCN is not used."""

USHCN_offset_start_year = 1980
"""The first year of the period considered when calculating the
offsets between GHCN and USHCN records, to apply to a USHCN
temperature series.  Irrelevant in the usual analysis using GHCN v3,
as USHCN is not used.
"""

USHCN_offset_max_months = 10
"""The maximum number of records used when computing a monthly offset
between GHCN and USHCN, to apply to a USHCN temperature series.  The
algorithm starts in the current year and works back to
*USHCN_offset_start_year*, computing an offset for each month using up
to *USHCN_offset_max_months* records.  Irrelevant in the usual
analysis using GHCN v3, as USHCN is not used."""

retain_contiguous_US = True
"""Whether to retain GHCN records that are in the contiguous US.  When
using USHCN ordinarily, USHCN replaces GHCN records, but: some
stations have duplicates; some GHCN stations have no USHCN
counterpart.  This parameter affects those records.  When not using
USHCN, this should be True.
"""

# STEP 1
#
# The following parameters control behaviour which used to be a part
# of GISTEMP step 1.

combine_records=False
"""Whether to attempt to combine station records.  The first part of
STEP 1 used to attempt to combine any non-overlapping records from the
same station into a combined station record, and then to combine
overlapping records.  This was dropped by GISTEMP v3 in 2011-12, when
GHCN-M 3 was adopted as the source data set: this sort of combination
is done in the preparation of GHCN-M 3.
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

