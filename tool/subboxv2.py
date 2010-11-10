#!/usr/bin/env python
# $URL$
# $Rev$
#
# subboxv2.py
#
# David Jones, Ravenbrook Limited, 2010-03-12

import extend_path

# Clear Climate Code
import io

def convert(inp, out):
    """Convert a file inp from subbox to V2 mean format."""

    # Clear Climate Code
    from code import eqarea

    v2 = io.V2MeanWriter(file=out, scale=0.01)

    subbox = iter(io.SubboxReader(inp))
    # First record is metadata, which we ignore.
    subbox.next()
    for record in subbox:
        lat,lon = eqarea.centre(record.box)
        record.uid = '%+05.1f%+06.1fC' % (lat,lon)
        v2.write(record)
    v2.close()

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if len(arg) <= 1:
        out = sys.stdout
    else:
        out = open(arg[1])
    
    return convert(open(arg[0]), out)

if __name__ == '__main__':
    main()
