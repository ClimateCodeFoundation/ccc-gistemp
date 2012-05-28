#! /usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/parameters.py $
# $Rev: 618 $
#
# parameters/extensions.py
#
# Nick Barnes, Ravenbrook Limited, 2010-02-15

"""Parameters controlling ccc-gistemp extensions to the standard
GISTEMP algorithm.

Parameters controlling the standard GISTEMP algorithm, or obsolete
features of GISTEMP, are in other parameter files.
"""
__docformat__ = "restructuredtext"

data_sources = "ghcn hohenpeissenberg scar"
"""Data sources that are used for the analysis (space separated string).
'ghcn' is the Global Historical Climate Network (NOAA/NCDC), version 3;
'hohenpeissenberg' is Wiljens data for Hohenpeissenberg;
'scar' is the READER data from Scientific Committee on Antarctic
Research.

These sources are no longer used by GISTEMP but should still work with
ccc-gistemp:

'ghcn.v2' is GHCN version 2;
'ushcn' is the United Stated Historical Climate Network (NOAA), version 2;
"""

augment_metadata = ''
"""(In the usual analysis this parameter is empty) This parameter enables
additional metadata fields to be read from a file.  The format is
"name=colA,colB,colC" (with an arbitrary number of comma separated
columuns).  The file called *name* is opened; each row is a comma
separated sequence of field values, with the fields being 'colA',
'colB', and so on.  There must be exactly one column called 'uid'.
If a station with the same uid appears in the ordinary metadata
(usually sourced from the v2.inv file) then the extra fields are
associated with the station.
"""

work_file_format = "v2"
"""The format of the intermediate files written to the 'work' directory:
'v2' for GHCN v2, 'v3' for GHCN v3.
"""
