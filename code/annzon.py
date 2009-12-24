#!/usr/bin/env python
# $URL$
# $Rev$
#
# annzon.py
#
# David Jones, Ravenbrook Limited, 2009-12-24

"""
annzon.py computes annual anomalies from monthly anomalies (and also
recomputes global and hemispheric averages).

As input it takes the output of zonav.py (a ZON file).

This is part of GISTEMP STEP5.

This is a placeholder implementation that just calls out to a
(precompiled) Fortran program.
"""

def main():
    import os
    os.system('GFORTRAN_CONVERT_UNIT="big_endian:10,11,12" bin/annzon.exe  > log/annzon.Ts.ho2.GHCN.CL.PA.log')

if __name__ == '__main__':
    main()
