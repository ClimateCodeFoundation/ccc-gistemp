#!/usr/bin/env python
"""Module that knows how to run step 2.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import os
import sys

import extend_path

import code.step2
import tool.giss_io

import tool.step1
from tool.fork import fork


def get_inputs(steps=(), save_work=True):
    if 1 in steps:
        source = tool.step1.get_step_iter(steps, save_work)
        if save_work:
            sink = tool.step1.get_outputs()
            source = fork(source, sink)
        return source
    return tool.giss_io.StationTsReader("work/Ts.txt")


def get_step_iter(steps=(), save_work=True):
    return code.step2.step2(get_inputs(steps, save_work))


def get_outputs():
    return tool.giss_io.StationRecordWriter(open("work/Ts.GHCN.CL.PA", "wb"),
            bos='<')


def main(argv=None):
    record_sink = get_outputs()
    record_source = get_step_iter()

    for record in record_source:
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main()

