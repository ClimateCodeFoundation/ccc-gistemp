#!/usr/bin/env python
# $URL$
# $Rev$
# Copyright (C) 2008 Ravenbrook Limited.
# David Jones
#
# preflight.py
#
# Pre-flight checks to avoid obvious mistakes that can be reasonably
# spotted early.  Mostly missing input files.  That sort of thing.

# http://www.python.org/doc/2.4.4/lib/os-process.html
import os
import sys

def checkit(log):
    """Perform pre-flight checks.  Messages are written to log, which
    could reasonably be something like sys.stderr.
    """

    def missing_input(name):
        """Complain about a missing input file."""

        log.write('MISSING: input/%s\n' % name)

    def missing_config(name):
        """Complain about a missing config file."""

        log.write('MISSING: config/%s\n' % name)

    def missing_big_files(list):
        """Check list (of names) and return mising ones.  .Z extension
        is also tried if supplied name is not found.
        """

        missing = []
        for name in list:
            if not any(input_ok(name+suffix) for suffix in ['','.Z','.gz']):
                missing_input(name)
                missing.append(name)
        return missing

    rc = 0

    step0 = """
        antarc1.list
        antarc1.txt
        antarc2.list
        antarc2.txt
        antarc3.list
        antarc3.txt
        t_hohenpeissenberg_200306.txt_as_received_July17_2003
        ushcn2.tbl
        ushcnV2_cmb.tbl
        """.split()
    for name in step0:
        if not input_ok(name):
            missing_input(name)
            rc = max(rc, 2)

    step0big = '9641C_200907_F52.avg v2.mean'.split()
    step5big = 'SBBX.HadR2'.split()
    big = step0big + step5big
    assert big
    missing = missing_big_files(big)
    if missing:
        log.write('Attempting to fetch missing files: %s\n' %
            ' '.join(missing))
        # Call fetch.py as if it were a program.
        import fetch
        fetch.main(argv=['fetch'] + missing)
        missing = missing_big_files(big)
        if missing:
            log.write("PROBLEM: Tried fetching missing files but it didn't work.\n")
            rc = max(rc, 2)
        else:
            log.write('OK: Fetching missing files looks like it worked\n')

    step1 = """
        mcdw.tbl
        ushcn2.tbl
        sumofday.tbl
        v2.inv
        """.split()
    step1_config = """
        combine_pieces_helena.in
        Ts.strange.RSU.list.IN
        Ts.discont.RS.alter.IN
        """.split()
    for name in step1:
        if not input_ok(name):
            missing_input(name)
            rc = max(rc, 2)
    for name in step1_config:
        if not config_ok(name):
            missing_config(name)
            rc = max(rc, 2)

    return rc

def input_ok(name):
    """Test if an input file is okay.  Checks it can open it for
    reading.  Return a true value if okay, false otherwise.
    """

    return file_readable('input/' + name)

def config_ok(name):
    return file_readable('config/' + name)

def file_readable(n):
    """Return True if file name n exists, can be opened for reading, and
    1 byte can be read (this avoids problems with 0 length files).
    """

    try:
        f = open(n, 'rb')
        if f.read(1) == '':
            return False
    except IOError:
        return False
    f.close()
    return True


if __name__ == '__main__':
    sys.exit(checkit(sys.stderr))
