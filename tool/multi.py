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
    and arg[2] should be the names of two different log/step3.log files."""

    fs = map(open, arg[1:3])
    stations = map(celldict, fs)
    if set(stations[0]) != set(stations[1]):
        print "Sets of cells differ"
    common = set(stations[0]) & set(stations[1])
    print "%d common cells" % len(common)
    for cell in common:
        if stations[0][cell] != stations[1][cell]:
            ina = stations[0] - stations[1]
            inb = stations[1] - stations[0]
            if ina:
                print "Cell %s %s only " % (cell, fs[0].name,
                  ' '.join(sorted(ina)))
            if inb:
                print "Cell %s %s only " % (cell, fs[1].name,
                  ' '.join(sorted(inb)))

def celldict(f):
    """From the (open) step3.log file *f* return a dict that maps
    from cell identifier (12 characters) to a set of stations."""

    result = {}
    for row in f:
        row = row.split(' ', 2)
        if row[1] != 'stations':
            continue
        stations = set(station for station,weight_ in json.loads(row[2]))
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
