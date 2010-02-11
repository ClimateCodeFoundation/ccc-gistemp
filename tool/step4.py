#!/usr/bin/env python
"""Script used to run step4 in isolation.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import extend_path

import code.step4
import tool.giss_io


def get_inputs(steps=(), save_work=True):
    f = open('input/SBBX.HadR2', 'rb')
    return tool.giss_io.SubboxReader(f)


def get_step_iter(steps=(), save_work=True):
    return code.step4.step4(get_inputs(steps))


def get_outputs():
    f = open('work/SBBX.HadR2', 'wb')
    return tool.giss_io.SubboxWriter(f)


def main(argv=None):
    record_sink = get_outputs()
    record_source = get_step_iter()

    for record in record_source:
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main(argv=sys.argv)
