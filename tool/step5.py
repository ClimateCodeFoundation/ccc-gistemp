#!/usr/bin/env python
"""Script used to run step5 in isolation.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import os
import sys

import extend_path

import code.step5
import tool.giss_io

import tool.step3
import tool.step4
from tool.fork import fork


def get_inputs(steps=(), save_work=True):
    outputs = []
    if 3 in steps:
        land_source = tool.step3.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step3.get_outputs()
            land_source = fork(land_source, sink)
    else:
        print "Reading from work/SBBX1880.Ts.GHCN.CL.PA.1200"
        land = open(os.path.join('work', 'SBBX1880.Ts.GHCN.CL.PA.1200'), 'rb')
        land_source = tool.giss_io.SubboxReader(land)

    if 4 in steps:
        ocean_source = tool.step4.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step4.get_outputs()
            ocean_source = fork(ocean_source, sink)
    else:
        print "Reading from work/SBBX.HadR2"
        ocean = open("work/SBBX.HadR2", "rb")
        ocean_source = tool.giss_io.SubboxReader(ocean)

    return land_source, ocean_source


def get_step_iter(steps=(), save_work=True):
    return code.step5.step5(get_inputs(steps, save_work))


def get_outputs():
    # Step 5 writes its output files directly.
    pass


def main(argv=None):
    record_sink = get_outputs()
    record_source = get_step_iter()

    return code.step5.step5(inputs=(land_reader, ocean_reader))


if __name__ == '__main__':
    main()

