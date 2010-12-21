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

    def missing_file(name):
        for suffix in ['', '.gz']:
            if input_ok(name+suffix):
                return False
        return True

    def missing_files(names):
        """Check list of names and return missing ones.  Compression
        extensions (.gz) are also tried if supplied name is not found.
        Each name can be a string or a list of strings; if a list of
        strings then any string will do, and the first string will be
        returned if none are present.
        """

        missing = []
        for item in names:
            if isinstance(item, str):
                if missing_file(item):
                    # Didn't find file, even trying compressed suffixes.
                    missing_input(item)
                    missing.append(item)
            else:
                assert isinstance(item, list)
                for name in item:
                    if not missing_file(name):
                        break
                else:
                    missing_input(item[0])
                    missing.append(item[0])
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
    step1 = """
        mcdw.tbl
        ushcn2.tbl
        sumofday.tbl
        v2.inv
        """.split()
    step4 = """
    oisstv2_mod4.clim.gz
    """.split()

    ushcn_alternatives = """
        ushcnv2
        9641C_201003_F52.avg
        9641C_201002_F52.avg
        9641C_200907_F52.avg
        """.split()

    step0big = 'v2.mean'.split() + [ushcn_alternatives]
    step5big = 'SBBX.HadR2'.split()
    big = step0big + step5big
    assert big
    all = step0 + step1 + step4 + big
    missing = missing_files(all)
    if missing:
        log.write('Attempting to fetch missing files: %s\n' %
            ' '.join(missing))
        # Call fetch.py as if it were a program.
        import fetch
        fetch.main(argv=['fetch'] + missing)
        missing = missing_files(all)
        if missing:
            log.write("PROBLEM: Tried fetching missing files but it didn't work.\n")
            rc = max(rc, 2)
        else:
            log.write('OK: Fetching missing files looks like it worked\n')

    return rc

def input_ok(name):
    """Test if an input file is okay.  Checks it can open it for
    reading.  Return a true value if okay, false otherwise.
    """

    return file_readable('input/' + name)

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
