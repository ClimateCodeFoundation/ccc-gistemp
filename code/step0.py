#! /usr/bin/env python
# step0.py
# $URL$
# $Rev$
#
# step1.py
#
# Nick Barnes and David Jones.
# Copyright (C) Ravenbrook Limited, 2008-2010.
# Copyright (C) Climate Code Fonudation, 2010-2012.
#
# BSD license, see license.txt

"""
Python code for the STEP0 part of the GISTEMP algorithm: combining
diverse inputs into a single dataset.
"""

# http://docs.python.org/release/2.4.4/lib/module-os.path.html
import os

def correct_Hohenpeissenberg(ghcn_records, hohenpeissenberg_dict):
    """Replace Hohenpeissenberg data from 1880 to 2002 in the GHCN
    dataset with the priv. comm. data."""
    print "Correct the GHCN Hohenpeissenberg record."

    # We expect the hohenpeissenberg_dict to contain a single record.
    # The assignment will raise an exception if this assumption fails.
    (key,hohenpeissenberg), = hohenpeissenberg_dict.items()
    del hohenpeissenberg_dict[key]
    cut = hohenpeissenberg.last_year + 1

    for record in ghcn_records.itervalues():
        if record.station_uid == hohenpeissenberg.station_uid:
            # Extract the data for the years more recent than the priv.
            # comm. data.
            new_data = []
            for year in range(cut, record.last_year + 1):
                new_data.extend(record.get_a_year(year))

            if record.uid == hohenpeissenberg.uid:
                # If the record has the same UID as Hohenpeissenberg then
                # replace the data with Hohenpeissenberg plus the recent
                # years.
                last_year = record.last_year
                record.set_series(hohenpeissenberg.first_month,
                                  hohenpeissenberg.series)
                for i, year in enumerate(range(cut, last_year + 1)):
                    record.add_year(year, new_data[i * 12:(i + 1) * 12])
            else:
                # For other records, just replace the data with the later
                # years.
                record.set_series(cut * 12 + 1, new_data)

def step0(input):
    """An iterator for Step 0.  Produces a stream of
    `giss_data.Series` instances.  *input* should be an instance that
    has an open() method.  input.open(x) is called for each data source x.
    (typically, this input object is made by the tool.io.step0_input()
    function).
    """

    # Read each data input into dictionary form.
    data = {}
    for source in input.sources:
        print "Load %s records" % source.upper()
        records = input.open(source)
        data[source] = dict((record.uid, record) for record in records)

    # If we're using GHCN and Hohenpeissenberg then we correct one with the other.
    if 'ghcn' in data and 'hohenpeissenberg' in data:
            correct_Hohenpeissenberg(data['ghcn'], data['hohenpeissenberg'])

    # Join all data sources together.
    records = {}
    for source in input.sources:
        records.update(data[source])

    # We sort here - as does GISTEMP - so that all the records for a
    # given 11-digit station ID are grouped together, ready for
    # combining in the next step.
    for _, record in sorted(records.iteritems()):
        if record:
            yield record
