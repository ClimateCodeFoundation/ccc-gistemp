#! /usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/step1.py $
# $Rev: 193 $
# 
# read_config.py
#
# Nick Barnes, Ravenbrook Limited, 2010-01-16

"""
Python code to read the various config and station files used by
GISTEMP:
"""

import struct

def v2_get_sources():
    """Reads the three tables mcdw.tbl, ushcn2.tbl, sumofday.tbl and
    return a dictionary that maps from 12-digit (string) station ID to
    the source (which is one of the strings 'MCDW', 'USHCN2',
    'SUMOFDAY').
    """

    sources = {}
    for source in ['MCDW', 'USHCN2', 'SUMOFDAY']:
        for line in open('input/%s.tbl' % source.lower()):
            _, id, rec_no = line.split()
            sources[id + rec_no] = source
    return sources


def v2_get_info():
    """Open the input/v2.inv file and convert it into a dictionary, *d*,
    that is returned.

    d[id] contains the metadata of the station with identification *id*
      where *id* is an 11-digit string.

    For a list of metadata fields, see *keys*.

    The input file is in the same format as the GHCN V2 file
    v2.temperature.inv (in fact, it's the same file, but with records
    added for the Antarctic stations that GHCN doesn't have).  The best
    description of that file's format is the Fortran program:
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.read.inv.f

    Here is a typical line, with a record diagram

    40371148001 ALMASIPPI,MA                    49.55  -98.20  274  287R   -9FLxxno-9x-9COOL FIELD/WOODSA1   0

    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    id---------xname--------------------------xlat---xlon----x1---2----34----5-6-7-8-910grveg-----------GU

    id: 40371148001
    name: ALMASIPPI,MA
    lat: 49.55
    lon: -98.20
     1: elevs: 274
     2: elevg: 287R
     3: pop: R
     4: ipop: -9
     5: topo: FL
     6: stveg: xx
     7: stloc: no
     8: iloc: -9
     9: airstn: x
    10: itowndis: -9
    grveg: COOL FIELD/WOODS
     G: GHCN-brightness: A
     U: US-brightness:1
    """

    keys = 'id name lat lon elevs elevg pop ipop topo stveg stloc iloc airstn itowndis grveg GHCN-brightness US-brightness'.split()
    info = {}
    for row in open('input/v2.inv', 'r'):
        values = struct.unpack('11sx30sx6sx7sx4s5sc5s2s2s2s2sc2s16scc', row[:102])
        # Would use struct.unpack_from if we knew we had Python 2.5
        hash = dict(zip(keys, values))
        hash['ipop'] = int(hash['ipop'])
        hash['lat'] = float(hash['lat'])
        hash['lon'] = float(hash['lon'])
        info[hash['id']] = hash
    return info


def get_helena_dict():
    """Reads the file config/combine_pieces_helena.in into a dict,
    mapping a station id to a tuple (ID with duplicate marker, year,
    month, summand)."""
    
    helena_ds = {}
    for line in open('config/combine_pieces_helena.in', 'r'):
        id, _, year, month, summand = line.split()
        helena_ds[id] = (int(year), int(month), float(summand))
    return helena_ds


def get_changes_dict():
    """Reads the file config/Ts.strange.RSU.list.IN and returns a dict
    result.  Each line in that file begins with a 12-digit station ID
    - actually the tuple (country-code, WMO station, modifier,
    duplicate) - and ends with either yyyy/mm, specifying a month
    datum to omit or with xxxx-yyyy, specifying years to omit.  xxxx
    can be 0, meaning from the beginning. yyyy can be 9999, meaning to
    the end.  The dict is a map from ID to ('month',yyyy,mm) or
    ('years',xxxx,yyyy).
    """

    dict = {}
    for line in open('config/Ts.strange.RSU.list.IN', 'r'):
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


def get_alter_dict():
    """Reads the file config/Ts.discont.RS.alter.IN into a dict.  Each
    line has a 12 digit station ID, a month, a year, and a
    floating-point temperature delta.  The dict maps the ID to (month,
    year, delta).
    """

    dict = {}
    for line in open('config/Ts.discont.RS.alter.IN'):
        id, month, year, num = line.split()
        dict[id] = [int(month), int(year), float(num)]
    return dict
