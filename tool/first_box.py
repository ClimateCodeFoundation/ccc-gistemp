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
    """
    Return a sequence of (box, first_year) pairs for each of 80
    boxes (from 0 to 79). If a box has no data (so no first
    year) then `first_year` is None for that box.
    """

    def get_box(p):
        return p[0][0]
    # First year for each box
    starts = [None] * 80
    for box, data in itertools.groupby(box_series(inp), get_box):
        starts[box] = list(data)[0][0][1]
    return enumerate(starts)

def multi_box_files(inps):
    result = []
    for inp in inps:
        result.append(box_starts(inp))
    return result

def as_text(datas):
    for pairs in zip(*datas):
        boxes = set(p[0] for p in pairs)
        # All pairs should have the same box number.
        assert len(boxes) == 1
        (box,) = boxes
        years = [p[1] for p in pairs]
        # To be commensurate with [HANSENLEBEDEFF1987] figure 2,
        # we start numbering boxes from 1.
        out = sys.stdout
        out.write("%d" % (box+1))
        for year in years:
            out.write(" %d" % year)
        out.write("\n")

def as_html(datas):
    # A sequence of pairs, one pair for each input file.
    zippers = itertools.izip(*datas)
    n_boxes_in_band = [4, 8, 12, 16, 16, 12, 8, 4]

    print '<table border="1" style="border-collapse:collapse; border-width:1px; border-color:#cccccc; border-style:solid; font-size:small">'

    for band, n_boxes in enumerate(n_boxes_in_band):
        box_span = 48 / n_boxes
        print "<tr>"
        for i in range(n_boxes):
            pairs = next(zippers)
            box_index = set(p[0] for p in pairs)
            assert len(box_index) == 1
            box_index, = box_index
            starts = [str(p[1]) for p in pairs]
            box_content = '<br />'.join(starts)
            print '<td align="center" colspan="%d">' % (box_span, )
            print box_content
            print '<span class="label" style="text-align: right; display: block; font-size: smaller">%d</span>' % (box_index + 1)
            print '</td>'
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

    if not arg:
        box_files = [glob.glob(os.path.join('result', 'landBX*'))[0]]
    else:
        box_files = arg
    sys.stderr.write("Using %s\n" % box_files)
    inps = [open(f, 'rb') for f in box_files]
    show(multi_box_files(inps))


if __name__ == '__main__':
    main()
