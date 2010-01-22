CLEAR CLIMATE CODE GISTEMP README FOR RELEASE 0.2.0

Nick Barnes, Clear Climate Code

$Date$


CONTENTS

  1. Introduction
  2. Installation
  A. References
  B. Document history
  C. Copyright and license

1. INTRODUCTION

This is release 0.2.0 of the Clear Climate Code GISTEMP project
(ccc-gistemp).

Clear Climate Code is reimplementing GISTEMP (the GISS surface
temperature analysis system) in Python, to make it clearer.

CCC GISTEMP release 0.2.0 is a release of CCC GISTEMP version 0.2.
The purpose of version 0.2 is: a first all-Python version.

Large amounts of the original GISS code was written in Fortran, some in
C, and some in ksh (and some in Python).  This code has all been
replaced with equivalent Python code.

The high-level step and sub-step architecture of GISTEMP is still
present in the code.  Lots of intermediate files are produced and then
immediately consumed; a lot of this style will be superseded in
future.

Some of the Python code is pretty clear, some is still very obscure.

The Python code is quite slow.  Optimizing for speed would be very
premature at the moment.

Some bugs have been reported to GISS and fixed in GISTEMP and in CCC.

URLs for further information:

http://clearclimatecode.org/ Clear Climate Code website and blog.
http://ccc-gistemp.googlecode.com/ ccc-gistemp code repository.


2. DEPENDENCIES

The only dependency is Python 2.5.2 (or a more recent Python version).
The code should run on OS X, FreeBSD, Windows, and probably a variety of
other Unix-like operating systems.  A network connection is required to
download the input files (which need only be done once), and to produce
an optional graph from the results.

Python may already be installed on your machine (for example, it comes
pre-installed on OS X), it may be possible to install it using your
operating system's package manager; for Windows you can download an
installer from http://www.python.org/download/ .  We recommend you use a
stable production release from the Python 2.x series (Python 3.x will
not work).


3. INSTALLATION

Unpack ccc-gistemp-0.2.0.tar.gz.


4. INPUT DATA

ccc-gistemp-0.2.0 uses input data in the subdirectory input/.  This
input data includes large files of temperature records from GHCN,
USHCN, and sea surface data, and small files of additional temperature
records and records and station tables from GISS.  ccc-gistemp-0.2.0
includes code (tool/preflight.py) to fetch this data from the
originating organisations over the internet.  It will not download a
file if it is already present in the input/ directory, so if you wish
to run ccc-gistemp with updated input data, you can delete the input
directory before you start.


5. RUNNING

To run CCC-GISTEMP:

python tool/run.py

To run only one step (0-5), give a -s <step> argument.  For instance:

python tool/run.py -s 1

We use this directory structure:

ccc-gistemp-x.x.x/code/     Source code only
                 /config/   Configuration files
                 /doc/      Internal developer documentation
                 /input/    Input data files
                 /log/      Log files
                 /tool/     Tools - sources other than the core of ccc-gistemp
                 /work/     Intermediate data files
                 /result/   Final result files

Running the code should write to the input/ directory when fetching
input data, but subsequently only write to the work/ log/ and result/
directories.  Before running tool/run.py, these directories can all be
deleted (if you wish, for example, to have a clean run).


6. RESULTS

After running run.py, the GISTEMP result files are all in the result/
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

2010-01-22 NB  Updated to reflect some code moving to tool/.
2010-01-11 NB  Updated to describe preflight better.
2010-01-06 DRJ Updated for our all-Python status.
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

$URL$
$Rev$
