#! /usr/bin/env python
# $URL$
# $Rev$
#
# read_config.py
#
# Nick Barnes, Ravenbrook Limited, 2010-01-16

"""
Python code to read the various config and station files used by
GISTEMP:
"""

def v2_get_sources():
    """Reads the three tables mcdw.tbl, ushcn3.tbl, sumofday.tbl and
    return a dictionary that maps from 12-digit (string) station ID to
    the source (which is one of the strings 'MCDW', 'USHCN2',
    'SUMOFDAY').
    """

    sources = {}
    for source in ['MCDW', 'USHCN3', 'SUMOFDAY']:
        for line in open('input/%s.tbl' % source.lower()):
            _, id11, duplicate = line.split()
            sources[id11 + duplicate] = source
    return sources


def step1_adjust():
    """Reads the file config/step1_adjust into a dict,
    mapping a record identifier to a tuple (year, month, summand).
    By convention the month is 1 for January."""

    adjust = {}
    for line in open('config/step1_adjust', 'r'):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        id, _, year, month, summand = line.split()
        adjust[id] = (int(year), int(month), float(summand))
    return adjust


def get_changes_dict():
    """Reads the file input/Ts.strange.v3.list.IN_full and returns a
    dict result.  Each line in that file begins with a 12-digit
    station ID - actually the tuple (country-code, WMO station,
    modifier, duplicate) - and ends with either yyyy/mm, specifying a
    month datum to omit or with xxxx-yyyy, specifying years to omit.
    xxxx can be 0, meaning from the beginning. yyyy can be 9999,
    meaning to the end.  The dict is a map from ID to
    ('month',yyyy,mm) or ('years',xxxx,yyyy).
    """

    dict = {}
    for line in open('input/Ts.strange.v3.list.IN_full', 'r'):
        split_line = line.split()
        id = split_line[0]
        try:
            year1, year2 = map(int, split_line[-1].split("-"))
            val = ("years", year1, year2)
        except ValueError:
            year, month = map(int, split_line[-1].split("/"))
            val = ("month", year, month)
        dict[id] = dict.get(id,[])
        dict[id].append(val)
    return dict

