#!/usr/bin/env python
# $URL$
# $Rev$
#
# rectmask.py
#
# David Jones, Climate Code Foundation, 2010-12-16

"""
rectmask.py [--latitude south,north] [--longitude west,east]

Create a step5mask file (one row for each of 8000 cells) in the "shape"
of a rectangle.  The output file has a "1.000" when the cell centre is
inside the specified rectangle, and "0.000" otherwise.

If only --latitude is specified, then a zonal strip is created; if only
--longitude is specified then meridional strip is created.  If neither
option is used the mask has a "1.000" everywhere.

The output file can then be copied to the input/ directory and this can
be used to mask the land series in Step 5 of the ccc-gistemp analysis.
"""

import extend_path

def rectangle(out, latbound=(-90.0,+90.0), lonbound=(-180.0,+180.0)):
    """*latbound* should be a pair of (southernbound, northernbound),
    *lonbound* should be a pair of (westernbound, easternbound).
    On *out* will be written a step5mask file where every cell whose
    centre is within the specified rectangle will be marked as
    "1.000" and every other cell marked as "0.000".
    """

    from code import eqarea
    from code import giss_data

    s,n = latbound
    w,e = lonbound

    for cell in eqarea.grid8k():
        lat,lon = eqarea.centre(cell)
        if s <= lat < n and w <= lon < e:
            m = 1.0
        else:
            m = 0.0
        out.write("%sMASK%.3f\n" % (giss_data.boxuid(cell), m))


def main(argv=None):
    import getopt
    import sys
    if argv is None:
        argv = sys.argv
    try:
        key = {}
        opts,arg = getopt.getopt(argv[1:], '', ['latitude=', 'longitude='])
        for o,v in opts:
            if o == '--latitude':
                key['latbound'] = tuple(float(x) for x in v.split(','))
            if o == '--longitude':
                key['lonbound'] = tuple(float(x) for x in v.split(','))
    except getopt.GetoptError:
        print __doc__
        return 2

    return rectangle(out=sys.stdout, **key)

if __name__ == '__main__':
    main()
