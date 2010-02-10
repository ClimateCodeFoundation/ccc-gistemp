#!/usr/bin/env python
"""Script used to run step3 in isolation.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import sys
import getopt

import extend_path

import code.step3
import tool.giss_io

import tool.step2
from tool.fork import fork


# TODO: What should I do?
def do_step3(label='GHCN.CL.PA', radius=1200, audit=None, subbox=()):
    """Do STEP3 of the GISTEMP algorithm.  label specifies the common part
    of the filenames for the 6 input files.  The filenames that are
    actually used as inputs will be formed by prepending "Ts." to label
    and appending ".n" where n is an integer from 1 to 6.  radius is the
    RCRIT value used in the algorithm, measured in Kilometres.
    """

    # Label string including radius
    labelr = 'Ts.%(label)s.%(radius)d' % locals()

    # Open the input station record source
    path = 'work/Ts.%s' % (label)
    f = open(path, "rb")
    reader = tool.giss_io.StationReader(f, bos='<')

    # Create the output file writer.
    subbox_grid_output = open('work/SBBX1880.%s' % labelr, 'wb')
    subbox_output = tool.giss_io.SubboxWriter(subbox_grid_output,
            trimmed=False)

    for el in code.step3.step3(reader, subbox_output, radius=radius,
            year_begin=1880, audit=audit, subbox=subbox):
        pass


# TODO: What should I do?
def old_main(argv=None):
    if argv is None:
        argv = sys.argv
    audit = None
    subbox = ()
    opt,arg = getopt.getopt(argv[1:], 'a:s:')
    for o,v in opt:
        if o == '-a':
            audit = eval(v, {'__builtins__':None})
        if o == '-s':
            subbox = map(int, v.split(','))
    return do_step3(audit=audit, subbox=subbox)


def get_inputs(steps=(), save_work=True):
    if 2 in steps:
        source = tool.step2.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step2.get_outputs()
            source = fork(source, sink)
        return source
    f = open("work/Ts.GHCN.CL.PA", "rb")
    return tool.giss_io.StationReader(f, bos='<')


def get_step_iter(steps=(), save_work=True):
    return code.step3.step3(get_inputs(steps, save_work))


def get_outputs():
    f = open('work/SBBX1880.Ts.GHCN.CL.PA.1200', 'wb')
    return tool.giss_io.SubboxWriter(f, trimmed=False)


def main(argv=None):
    record_sink = get_outputs()
    record_source = get_step_iter()
    record_sink.add_meta(record_source.meta)

    for record in record_source:
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main(argv=sys.argv)
