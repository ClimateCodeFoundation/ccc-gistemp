#!/usr/bin/env python
# 
# $Id: //info.ravenbrook.com/project/ccc/master/tool/unbdb.py#1 $
# 
# convert a GISTEMP STEP1 .bdb file to a .txt file suitable for
# diffing.  This allows me to compare our (CCC) intermediate data
# products with those from original GISTEMP.
#
# Based closely on bdb_to_text.py (added sorting the ids)

import re, stationstring, bsddb
def unbdb(file):
    txt_file = file + '.txt'
    db_file = file + '.bdb'
    db = bsddb.hashopen(db_file, 'r')
    print "reading", db_file
    f = open(txt_file, 'w')
    print "creating", txt_file
    ids = db['IDS'].split()
    ids.sort()
    count = 0
    for id in ids:
        count = count + 1
        if count % 1000 == 0:
            print count
        s = db[id]
        st = stationstring.new(s)
        f.write(st.to_text(id))
    print count
    f.close()
    db.close()
