#!/usr/bin/env python

"""
Show the first year that each box has data.
"""

# REFERENCES

# [HANSENLEBEDEFF1987]
#   "Global Trends of Measured Surface Air Temperature"
#   Journal of Geophysical Research, Vol. 92, No. D11,
#   Pages 13,345-13,372;
#   1987-11-20


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
        yield box, list(data)[0][0][1]

def as_text(data):
    for box, year in data:
        # To be commensurate with [HANSENLEBEDEFF1987] figure 2,
        # we start numbering boxes from 1.
        print box+1, year

def main(argv=None):
    if argv is None:
        argv = sys.argv

    box_file = glob.glob(os.path.join('result', 'landBX*'))[0]
    sys.stderr.write("Using %s\n" % box_file)
    with open(box_file, 'rb') as inp:
        as_text(box_starts(inp))


if __name__ == '__main__':
    main()
