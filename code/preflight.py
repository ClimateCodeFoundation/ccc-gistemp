#!/usr/bin/env python
# $Id$
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
            if not (input_ok(name) or input_ok(name + '.Z')):
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
        ushcn.tbl
        """.split()
    for name in step0:
        if not input_ok(name):
            missing_input(name)
            rc = max(rc, 2)

    step0big = 'v2.mean hcn_doe_mean_data'.split()
    assert step0big
    missing = missing_big_files(step0big)
    if len(missing) == len(step0big):
        log.write('Attempting to fetch missing files: %s\n' %
            ' '.join(missing))
        os.system('code/fetch.py')
        missing = missing_big_files(step0big)
        if missing:
            log.write("PROBLEM: Tried fetching missing files but it didn't work.\n")
            rc = max(rc, 2)
        else:
            log.write('OK: Fetching missing files looks like it worked\n')
    elif missing:
        log.write('PROBLEM: Some, but not all, of [%s] are missing.\n'
                  '         Please sort it out by hand.\n' %
                ' '.join(step0big))
        rc = max(rc, 2)

    step1 = """
        mcdw.tbl
        ushcn.tbl
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
    """Return True if file name n exists, and can be opened for reading."""

    try:
        f = open(n)
    except IOError:
        return False
    f.close()
    return True


if __name__ == '__main__':
    sys.exit(checkit(sys.stderr))
