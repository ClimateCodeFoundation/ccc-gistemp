#!/usr/bin/env python
# $URL$
# $Id$
#
# zonav.py
#
# David Jones, Ravenbrook Limited, 2009-12-11
#

"""
Perform Zonal Averaging.
"""

def main():
    import os
    os.system('GFORTRAN_CONVERT_UNIT="big_endian:10,11" bin/zonav.exe > log/zonav.Ts.ho2.GHCN.CL.PA.log')

if __name__ == '__main__':
    main()
