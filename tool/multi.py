#!/usr/bin/env python
# $URL$
# $Rev$

"""
Multi-tool.  A collection of miscellaneous commands and procedures.
Mostly intended for debugging and introspection.

Code is collected here to avoid cluttering up tool with many one-off
scripts.

Run "python tool/multi.py commands" to see command list.
"""

# http://docs.python.org/release/2.6.6/library/json.html
import json
import re

class Fatal(Exception):
    """An error occurred."""

def cellstations(arg):
    """[command] Compares the lists of stations for each cell.  arg[1]
    and arg[2] should be the names of two different log/step3.log files.
    """

    fs = map(open, arg[1:3])
    stations = map(celldict, fs)
    if set(stations[0]) != set(stations[1]):
        print "Sets of cells differ"
    common = set(stations[0]) & set(stations[1])
    print "%d common cells" % len(common)
    for cell in common:
        # Get the two collections of stations...
        # as lists...
        stationsa = stations[0][cell]
        stationsb = stations[1][cell]
        # as dicts...
        dicta = dict(stationsa)
        dictb = dict(stationsb)
        # and sets.
        seta = set(dicta)
        setb = set(dictb)
        if stationsa != stationsb:
            ina = seta - setb
            inb = setb - seta
            if ina:
                print "Cell %s %s only " % (cell, fs[0].name,
                  ' '.join(sorted(ina)))
            if inb:
                print "Cell %s %s only " % (cell, fs[1].name,
                  ' '.join(sorted(inb)))
            commonstations = seta & setb
            for station in commonstations:
                if dicta[station] != dictb[station]:
                    print "Cell %s station %s has different weights (%f and %f)." %(
                      cell, station, dicta[station], dictb[station])
            # Split into stations with 0 weight and those with non-zero
            # weight.
            zeroa = [s for s,w in stationsa if not w and s in commonstations]
            zerob = [s for s,w in stationsb if not w and s in commonstations]
            posa = [s for s,w in stationsa if w and s in commonstations]
            posb = [s for s,w in stationsb if w and s in commonstations]
            for i,(a,b) in enumerate(zip(posa, posb)):
                if a != b:
                    print "Cell %s order differs: %s and %s at index %d" % (
                       cell, a, b, i)
                break
            # Not normally reported.
            if False and zeroa != zerob:
                print "Cell %s (the unused stations are in different order; result not affected)" % (
                  cell)


def celldict(f):
    """From the (open) step3.log file *f* return a dict that maps
    from cell identifier (12 characters) to a list of (station,weight)
    pairs.
    """

    result = {}
    for row in f:
        row = row.split(' ', 2)
        if row[1] != 'stations':
            continue
        stations = json.loads(row[2])
        result[row[0]] = stations
    return result

def commands(arg):
    """[command] Show command list."""

    rx = re.compile(r'\[command] *')

    for name,f in members.items():
        doc = getattr(f, '__doc__')
        if doc and rx.match(doc):
            doc = rx.sub('', doc)
            doc = doc.replace('\n', '')
            doc = re.sub(r'[;.] +.*', '.', doc)
            print name, ':', doc

members = locals()

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if not arg:
        # With no args, run "commands" to show the command list.
        arg = ['commands']
    if arg[0] in members:
        members[arg[0]](arg)
    else:
        raise Fatal("Unknown command %s." % arg[0])

if __name__ == '__main__':
    main()
