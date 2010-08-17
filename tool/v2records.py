#!/usr/bin/env python
# $URL$
# $Rev$
#
# records.py

"""
Process a v2.mean file to find:
- longest record (latest year minus oldest year)
- shortest record
- largest range (highest minus lowest temperature)
- least range
- fattest (largest ratio of record length to range)
- skinniest (smallest ratio of record length to range)
"""

import itertools
import os
import struct
import sys

sys.path.append(os.path.join(os.getcwd(),'code'))

def collect(v2):
    """v2 should be a v2.mean file"""

    def id11(x):
        return x[:11]
    def year(x):
        """Given a line of v2.mean, return the year."""
        return int(x[12:16])

    longest = (0, '')
    shortest = (9999, '')
    largest = (0, '')
    flattest = (9999, '')
    fattest = (0, '')
    skinniest = (9999, '')
   
    for id11,lines in itertools.groupby(v2, id11):
        sys.stdout.write('\r' + id11 + ' ')
        sys.stdout.write(' '.join([longest[1], shortest[1], largest[1],
        flattest[1]]))
        sys.stdout.flush()
        lines = list(lines)
        yearmin = min(map(year, lines))
        yearmax = max(map(year, lines))
        length = yearmax - yearmin + 1
        tmin = 9999
        tmax = -9999
        for row in lines:
            data = struct.unpack('5s'*12, row[16:-1])
            data = map(int, data)
            data = filter(lambda x: x != -9999, data)
            tmin = min([tmin]+data)
            tmax = max([tmax]+data)
        range = tmax - tmin
        aspect = float(length)/float(range)
        if length > longest[0]:
            longest = (length, id11)
        if length < shortest[0]:
            shortest = (length, id11)
        if range > largest[0]:
            largest = (range, id11)
        if range < flattest[0]:
            flattest = (range, id11)
        if aspect > fattest[0]:
            fattest = (aspect, id11)
        if aspect < skinniest[0]:
            skinniest = (aspect, id11)

    sys.stdout.write('\n')
    return longest, shortest, largest, flattest, fattest, skinniest

if __name__ == '__main__':
    from step0 import open_or_uncompress
    print collect(open_or_uncompress('input/v2.mean'))
