#!/usr/bin/env python
# $URL$
# $Rev$
#
# topng.py

"""
Convert slides to PNG, using Inkscape.
"""

def topng():
    import glob
    import os

    for name in glob.glob('*.svg'):
        pngname = name.replace('.svg','.png')
        cmd = "inkscape --export-background=#ffffffff --export-png=%s %s" % (pngname, name)
        print cmd
        os.system(cmd)

def main():
    topng()

if __name__ == '__main__':
    main()
