#!/usr/bin/env python
# $URL$
# $Rev$

"""Module that knows how to run step 1.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import extend_path

import code.step1
import tool.giss_io

import tool.step0
from tool.fork import fork


def get_inputs(steps=(), save_work=True):
    if 0 in steps:
        source = tool.step0.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step0.get_outputs()
            source = fork(source, sink)
        return source
    return tool.giss_io.V2MeanReader("work/v2.mean_comb")


def get_step_iter(steps=(), save_work=True):
    return code.step1.step1(get_inputs(steps, save_work))


def get_outputs():
    return tool.giss_io.StationTsWriter("work/Ts.txt")


def main(argv=None):
    record_sink = get_outputs()
    for record in get_step_iter():
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main()
