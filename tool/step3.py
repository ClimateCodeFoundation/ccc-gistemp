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


def get_inputs(steps=(), save_work=True):
    if 2 in steps:
        source = tool.step2.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step2.get_outputs()
            source = fork(source, sink)
        return source
    f = open("work/Ts.GHCN.CL.PA", "rb")
    return tool.giss_io.StationReader(f, bos='<')


def get_step_iter(steps=(), save_work=True, audit=None, subbox=()):
    return code.step3.step3(get_inputs(steps, save_work),
            audit=audit, subbox=subbox)


def get_outputs():
    f = open('work/SBBX1880.Ts.GHCN.CL.PA.1200', 'wb')
    return tool.giss_io.SubboxWriter(f, trimmed=False)


def main(audit=None, subbox=()):
    record_sink = get_outputs()
    record_source = get_step_iter(audit=audit, subbox=subbox)

    for record in record_source:
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    audit = None
    subbox = ()
    opt,arg = getopt.getopt(sys.argv[1:], 'a:s:')
    for o,v in opt:
        if o == '-a':
            audit = eval(v, {'__builtins__':None})
        if o == '-s':
            subbox = map(int, v.split(','))
    main(audit=audit, subbox=subbox)
