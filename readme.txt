CLEAR CLIMATE CODE GISTEMP README FOR RELEASE x.x.x

Nick Barnes, Ravenbrook Limited

$Date: 2008/09/19 $


CONTENTS

  1. Introduction
  2. Installation
  A. References
  B. Document history
  C. Copyright and license

1. INTRODUCTION

This is release x.x.x of the Clear Climate Code GISTEMP project.

Clear Climate Code is reimplementing GISTEMP (the GISS surface
temperature analysis system) in Python, to make it clearer.

CCC GISTEMP release x.x.x is a release of CCC GISTEMP version x.x.
The purpose of version x.x is: to publish our early progress towards
the project goals.

There are still some pieces of GISS FORTRAN code, but much of the
FORTRAN, all the C, and all the ksh has been replaced with Python and
/bin/sh.

The high-level step and sub-step architecture of GISTEMP is still
present in the code.  Lots of intermediate files are produced and then
immediately consumed; a lot of this style will be superseded in
future.

Some of the Python code is pretty clear, some is still very obscure.

The Python code is quite slow.  Optimizing for speed would be very
premature at the moment.

Some bugs have been reported to GISS and fixed in GISTEMP and in CCC.


2. DEPENDENCIES

You will need Fortran 9x and Python (2.4, or a more recent and compatible
version; 2.4 Python syntax is used so it will definitely not work on
earlier versions without modification), on a more-or-less Unix platform
(there are two /bin/sh shell scripts using Unix commands such as grep,
sed, sort, ln, mv, cp, rm, echo, etc).

We have observed that GNU Fortran 95 4.1.2 has broken rounding, and
strongly advise you not to use it.

At the moment we believe some of our code only works on little-endian
machines (e.g. ones with x86 processors), for essentially bogus
reasons.  This dependency will go away soon.


3. INSTALLATION

Unpack ccc-gistemp-x.x.x.tar.gz.


4. INPUT DATA

ccc-gistemp-x.x.x uses input data in the subdirectory input/ This can
be obtained from the originating organisations over the internet by
running the Python script code/preflight.py

5. RUNNING

code/run.sh  runs steps 0, 1, 2, 3, and 5.

code/STEP4_5/do_comb_step4.sh runs step 4, which updates the
sea-surface temperature file input/SSBX.HadR2.  This is optional.

Both these scripts use this directory structure:

ccc-gistemp-x.x.x/code/     Source code only
                 /config/   Configuration files
                 /bin/      Fortran compiled executables
                 /input/    Input data files
                 /work/     Intermediate data files
                 /log/      Log files
                 /result/   Final result files

Running the code should only write to the bin/ work/ log/ and result/
directories.  Before running code/run.sh, these directories can all be
deleted.


6. RESULTS

After running run.sh, the GISTEMP result files are all in the result/
directory.  A rough graphical output is made using the Google Chart API;
this file:

    result/google-chart.url

contains the URL of a chart showing the global annual mean surface
temperature anomaly, and if a network connexion can be made then this
file:

    result/GLB.Ts.ho2.GHCN.CL.PA.png

contains the PNG image downloaded from that URL.


A. REFERENCES

None.


B. DOCUMENT HISTORY

Most recent changes first:

2009-12-03 NB  Updated for transfer to GoogleCode project.
2008-09-19 DRJ Added PNG result.
2008-09-13 NB  Updated for CCC 0.1.0.
2008-09-12 NB  Updated for CCC 0.0.3.
2008-09-12 NB  Updated for CCC 0.0.2.
2008-09-11 NB  Updated for CCC 0.0.1.
2008-09-08 NB  Created.

C. COPYRIGHT AND LICENSE

This document is copyright (C) 2009 Ravenbrook Limited.  All rights
reserved.

Redistribution and use of this document in any form, with or without
modification, is permitted provided that redistributions of this
document retain the above copyright notice, this condition and the
following disclaimer.

THIS DOCUMENT IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
DOCUMENT, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

$Id$
