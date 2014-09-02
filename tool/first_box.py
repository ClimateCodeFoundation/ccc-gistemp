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


import getopt
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

def as_html(data):
    # `data` is a sequence, with potentially missing data.
    box_first = dict(data)
    n_boxes_in_band = [4, 8, 12, 16, 16, 12, 8, 4]

    print '<table border="1" style="border-collapse:collapse; border-width:1px; border-color:#cccccc; border-style:solid; font-size:small">'

    box = 0
    for band, n_boxes in enumerate(n_boxes_in_band):
        box_span = 48 / n_boxes
        print "<tr>"
        for i in range(n_boxes):
            box_content = box_first.get(box, '')
            print '<td align="center" colspan="%d">' % (box_span, )
            print '%s</td>' % box_content
            box += 1
        print "</tr>"
    print "</table>"

def main(argv=None):
    if argv is None:
        argv = sys.argv

    show = as_text

    opt,arg = getopt.getopt(argv[1:], '', ['html'])
    for k,v in opt:
        if k == '--html':
            show = as_html

    box_file = glob.glob(os.path.join('result', 'landBX*'))[0]
    sys.stderr.write("Using %s\n" % box_file)
    with open(box_file, 'rb') as inp:
        show(box_starts(inp))


if __name__ == '__main__':
    main()
