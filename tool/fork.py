#!/usr/bin/env python
"""Support for forking iterated record flows.

This contains the single generator `fork`.

"""
__docformat__ = "restructuredtext"


def fork(source, sink):
    """Generator: All records in source are added to the sink.

    This is used to pipeline the CCC steps, whilst also writing the
    intermediate work files. This generator yields all records in the `source`
    and also passes them to the `sink` using the sink's `add_record` method.

    :Param source:
        This must be an iterable of some sort.
    :Param sink:
        This must provide an ``add_record`` method, which can accept
        every recorded yielded th the `source`.

    """
    for record in source:
        sink.add_record(record)
        yield record



