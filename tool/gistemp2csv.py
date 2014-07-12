#!/usr/bin/env python
#
# gistemp2csv.py
#
# Filipe Fernandes, 2011-08-01
# David Jones, 2014-07-11

"""
Convert ccc-gistemp result file in text format to
comma separated value format.
"""

# http://docs.python.org/release/2.4.4/lib/module-os.path.html
import os
# http://docs.python.org/release/2.4.4/lib/module-re.html
import re
# http://docs.python.org/release/2.4.4/lib/module-csv.html
import csv


def chunks(l, n):
    """
    Return n-sized chunks from the string (or list) *l*; the
    final chunk may have fewer than `n` elements.
    """
    return [l[i:i + n] for i in range(0, len(l), n)]


def non_zonal(line):
    """Break Non-Zonal file format into specified chunks."""

    # The first 12 entries (monthly anomaly) are each 5 characters wide.
    month = chunks(line[:65], 5)
    # After 3 spaces there are 2 entries (yearly anomaly
    # for Jan-Dec and Dec-Nov) each 4 characters wide.
    year = chunks(line[68:76], 4)
    # After 2 spaces there are 4 entries (seasonal anomaly)
    # each 5 characters wide.
    # 4 seasons, each 5 characters wide.
    season = chunks(line[78:98], 5)
    # Bundle into single row.
    row = month + year + season
    row = [x.strip() for x in row]
    return row


def gistemp2csv(fin, one_header=False):
    """
    Write a comma separated value file from ccc-gistemp text
    output; such files are typically found in the result/
    directory.

    if `one_header` is true then just one header line is written
    (the last one, appearing immediately before data).
    """

    if os.path.getsize(fin) == 0:
        # Skip empty text file.
        return

    # Pick an output filename by changing the extension to '.csv'
    path, fname = os.path.split(fin)
    basename, ext = os.path.splitext(fname)
    foutname = os.path.join(path, basename + '.csv')

    # Data for each row.  We also use this to detect the header
    # rows that come before the data.
    data = None

    with open(fin) as lines, open(foutname, 'wb') as fout:
        csv_out = csv.writer(fout, delimiter=',')
        for line in lines:
            # Data lines start with a 4 digit number (the year).
            if re.match(r"^\d{4}", line):
                # Non-Zonal avg files are treated differently.
                if not 'Zon' in basename:
                    data = non_zonal(line)
                else:
                    data = line.split()[:-1]
                if one_header and header:
                    csv_out.writerow(header)
                    header = None
                csv_out.writerow(data)
            # We assume all the header lines appear before any
            # data, so data will not be assigned yet.
            if data is None:
                # The header and blurb lines are somewhat
                # special cased.
                if 'Year' in line:
                    # Ignore final "Year" column (same as first column).
                    cells = line.split()[:-1]
                elif 'EQU' in line:
                    cells = ['', '', '', ''] + line.split()
                else:
                    cells = [line.strip()]
                header = cells
                if not one_header:
                    csv_out.writerow(cells)

    return foutname

def main(argv=None):
    # https://docs.python.org/release/2.4.4/lib/module-getopt.html
    import getopt
    import sys

    if argv is None:
        argv = sys.argv

    opt,arg = getopt.getopt(argv[1:], '', ['one-header'])
    keys = dict()
    for o,v in opt:
        if o == '--one-header':
            keys['one_header'] = True

    if not arg:
        import glob
        arg = glob.glob("*.txt")

    for f in arg:
        gistemp2csv(f, **keys)

if __name__ == '__main__':
    main()
