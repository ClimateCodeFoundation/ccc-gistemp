#!/usr/bin/env python
# $URL$
# $Rev$
#
# subboxv2.py
#
# David Jones, Ravenbrook Limited, 2010-03-12

import extend_path

# Clear Climate Code
import giss_io

def convert(inp, outpath):
    """Convert a file inp from subbox to V2 mean format."""

    # Clear Climate Code
    from code import eqarea

    v2 = giss_io.V2MeanWriter(path=outpath, scale=0.01)

    subbox = iter(giss_io.SubboxReader(inp))
    # First record is metadata, which we ignore.
    subbox.next()
    for record in subbox:
        lat,lon = eqarea.centre(record.box)
        record.uid = '%+05.1f%+06.1fC' % (lat,lon)
        v2.write(record)
    v2.close()

def main(argv=None):
    if argv is None:
        import sys
        argv = sys.argv

    arg = argv[1:]
    
    return convert(open(arg[0]), arg[1])

if __name__ == '__main__':
    main()
