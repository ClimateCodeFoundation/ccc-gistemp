#!/usr/bin/env python
# $URL$
# $Rev$

# SBBXtotext.py
#
# David Jones, Ravenbrook Limited
#

"""
Convert GISTEMP subbox file to text.

Similar job to subboxtotext.py, but this file uses the subbox module.
subboxtotext.py is more complete (it automagically handles
trimmed/untrimmed for example) but does not use the subbox module.
The plan is for this tool (SBBXtotext) to evolve into the preferred tool.
"""

# Ravenbrook
sys.path.append(os.path.join(os.getcwd(),'code'))
import subbox

def totext(inp, out, metaonly=False, bos='>', trimmed=None):
    """Display a GISTEMP subbox file in text format. *metaonly*, when
    True, means that only the metadata for each subbox record (basically
    its location on Earth and its length) are displayed.  *bos* specifies the
    Byte Order and Size: '>' for big-endian (GISS), '<' for
    little-endian. *trimmed* is not used yet.
    """

    boxes = subbox.File(inp, bos=bos)
    print >> out, repr(boxes)
    for cell in boxes:
        print >> out, cell
        if metaonly:
            continue
        print >> out, cell[:]

def main(argv=None):
    import getopt
    import sys

    if argv is None:
        argv = sys.argv

    metaonly = False
    bos = '>'
    trimmed = 'fromfile'
    opt, arg = getopt.getopt(argv[1:], 'b:mtu')
    for o,v in opt:
        if o == '-m':
            metaonly = True
        if o == '-b':
            bos = v
        if o == '-t':
            trimmed = True
        if o == '-u':
            trimmed = False
    if len(arg) == 0:
        totext(sys.stdin, metaonly=metaonly, bos=bos, trimmed=trimmed)
    else:
        for n in arg:
            totext(open(n, 'rb'), sys.stdout, metaonly=metaonly, bos=bos,
              trimmed=trimmed)

if __name__ == '__main__':
    main()
