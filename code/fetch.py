#!/usr/bin/env python
# $Id$
# David Jones.
# Copyright 2008 Ravenbrook Limited.

"""Script to fetch the inputs required for the GISTEMP program.  The
inputs on are documented in the gistemp.txt file:
http://clearclimatecode.org/master/test/GISTEMP/gistemp.txt
"""

# http://www.python.org/doc/2.4.4/lib/module-urllib.html
import urllib
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys

# See
# http://groups.google.com/group/ccc-gistemp-discuss/web/compiling-gistemp-source?hl=en
# (But station_inventory.Z corrected to station.inventory.Z)

noaa = """
ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.mean.Z
ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.temperature.inv
ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/hcn_doe_mean_data.Z
ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/station.inventory.Z
""".split()

all = noaa

output = sys.stdout

for url in all:
    def hook(n, bs, ts):
        output.write("\r%s %d" % (url, n*bs))
        output.flush()
    urllib.urlretrieve(url, url.split('/')[-1], hook)
    output.write('\n')
