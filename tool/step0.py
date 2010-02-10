#!/usr/bin/env python
"""Module that knows how to run step 0.

<+Detailed multiline documentation+>
"""
__docformat__ = "restructuredtext"

import itertools

import extend_path

import code.step0
import tool.giss_io


class Struct:
    pass


def get_inputs(steps=(), save_work=True):
    inputs = Struct()
    inputs.ushcn_source = tool.giss_io.USHCNReader(
            "input/9641C_200907_F52.avg", 'input/ushcn2.tbl')
    inputs.ghcn_source = tool.giss_io.V2MeanReader( "input/v2.mean")
    inputs.antarc_source = itertools.chain(
            tool.giss_io.AntarcticReader("input/antarc1.txt", 
                "input/antarc1.list", '8'),
            tool.giss_io.AntarcticReader("input/antarc3.txt",
                "input/antarc3.list", '9'),
            tool.giss_io.AustroAntarcticReader("input/antarc2.txt",
                "input/antarc2.list", '7'))
    inputs.hohenpeis_source = tool.giss_io.HohenpeissenbergReader(
            "input/t_hohenpeissenberg_200306.txt_as_received_July17_2003")

    return inputs


def get_step_iter(steps=(), save_work=True):
    return code.step0.step0(get_inputs())


def get_outputs():
    return tool.giss_io.V2MeanWriter("work/v2.mean_comb")


def main():
    record_sink = get_outputs()
    for record in get_step_iter():
        record_sink.add_record(record)
    record_sink.close()


if __name__ == '__main__':
    main()
