#!/usr/bin/env python

"""
Show the first year that each box has data.
"""

import glob
import itertools
import os
import sys

# ccc-gistemp

from compare_results import box_series

def box_starts(inp):
    def get_box(p):
        return p[0][0]
    for box, data in itertools.groupby(box_series(inp), get_box):
        print box, list(data)[0][0][1]

def main(argv=None):
    if argv is None:
        argv = sys.argv

    box_file = glob.glob(os.path.join('result', 'landBX*'))[0]
    sys.stderr.write("Using %s\n" % box_file)
    with open(box_file, 'rb') as inp:
        box_starts(inp)


if __name__ == '__main__':
    main()
