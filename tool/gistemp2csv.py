#!/usr/bin/env python
# $URL$
# $Rev$
#
# gistemp2csv.py
#
# Filipe Fernandes, 2011-08-01

"""Convert ccc-gistemp text files output to comma separated value format."""

# http://docs.python.org/release/2.4.4/lib/module-os.path.html
import os
# http://docs.python.org/release/2.4.4/lib/module-re.html
import re
# http://docs.python.org/release/2.4.4/lib/module-csv.html
import csv


def chunks(l, n):
    """
    return n-sized chunks from  the list l.
    """
    return [l[i:i + n] for i in range(0, len(l), n)]


def non_zonal(line):
    """Break Non-Zonal file format into specified chunks."""

    # The first 12 entries are 5 characters width columns.
    month = chunks(line[:65], 5)
    # After 3 white spaces it becomes 4 characters width columns.
    year = chunks(line[68:76], 4)
    # After 2 white spaces it is back to 5 characters width columns.
    # the -5 removes the last redundant year entry.
    season = chunks(line[78:-5], 5)
    # Bundle into single row.
    row = month + year + season
    row = [x.strip() for x in row]
    return row


def gistemp2csv(fin):
    """Write a comma separated value file from ccc-gistemp text
    output; such files are typically found in the result/
    directory."""

    path, fname = os.path.split(fin)
    basename, ext = os.path.splitext(fname)
    foutname = os.path.join(path, basename + '.csv')

    fout = open(foutname, 'wb')

    csv_out = csv.writer(fout, delimiter=',')

    header = []
    field1 = ''
    field2 = ''

    # Data for each row.  We also use this to detect the header
    # rows that come before the data.
    data = None

    with open(fin) as lines:
        for line in lines:
            # Data lines start with a 4 digit number.
            if re.match(r"^\d{4}", line):
                # Non-Zonal avg files are treated differently.
                if not 'Zon' in basename:
                    data = non_zonal(line)
                else:
                    data = line.split()[:-1]
                csv_out.writerow(data)
            # We assume all the header lines appear before any
            # data, so data will not be assigned yet.
            if data is None:
                cells = line.split()[:-1]
                if cells and 'Year' == cells[0] or 'EQU' in cells:
                    # The header rows; need to be arranged in cells.
                    csv_out.writerow(cells)
                else:
                    # Just a "blurb" line.  Pass as-is.
                    csv_out.writerow([line.strip()])

    return foutname

if __name__ == '__main__':
    import glob
    files = glob.glob("*.txt")
    for f in files:
        gistemp2csv(f)
