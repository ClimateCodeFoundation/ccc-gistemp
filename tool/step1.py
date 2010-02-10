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


def get_inputs(steps=()):
    if 0 in steps:
        return tool.step0.get_step_iter(steps)
    return tool.giss_io.V2MeanReader("work/v2.mean_comb")


def get_step_iter(steps=()):
    return code.step1.step1(get_inputs(steps))


def get_outputs():
    return tool.giss_io.StationTsWriter("work/Ts.txt")


def main(argv=None):
    record_sink = get_outputs()
    for record in get_step_iter():
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main()
