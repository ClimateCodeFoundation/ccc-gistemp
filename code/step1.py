#! /usr/bin/env python
# 
# step1.py
# 
# Python code reproducing the STEP1 part of the GISTEMP algorithm.
#
# Very much first draft, but does product the correct Ts.txt.
# Nick Barnes, Ravenbrook Limited, 2008-08-06
# 
# $Id: //info.ravenbrook.com/project/ccc/master/code/step1.py#8 $
#
# To run:
# $ ./step1.py
# $
# 
# Requires the following files in the input/ directory,
# from GISTEMP STEP1/input_files/:
# 
# mcdw.tbl
# ushcn.tbl
# sumofday.tbl
# v2.inv
#
# Do not be tempted to replace v2.inv with the apparently similar
# v2.temperature.inv file available from NOAA,
# ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.temperature.inv .  The
# NOAA file has been treated for GISTEMP's use by, for example, adding
# records corresponding to Antarctic stations that are not used in GHCN
# but are used in the GISTEMP analysis.  Step 1 (this step) expects to
# find a record in v2.inv for every station it has a time series for.
#
# Requires the following files in the config/ directory,
# from GISTEMP STEP1/input_files/:
# 
# combine_pieces_helena.in
# Ts.strange.RSU.list.IN
# Ts.discont.RS.alter.IN
#
# Also requires the existence of writeable work/ and log/ directories.

import math, sys, struct, bsddb, re, array

# progress counter, used when writing records into a database file.

def inc_count(count):
    count += 1
    if count % 100 == 0:
        print count
    return count

def FEQUALS(x, y):
  return abs(x - y) < 0.1

# From stationstring.py
#
# A StationString encapsulates a dict and a set of temperature records
# for a single station.  The set of temperature records is a 2D array
# of floating-point temperatures: data[m][y], where m ranges from 0 to
# 11 and y from 0 to N.  Missing values are recorded as 999.9.
#
# This "2D array" is actually a list of array.array values. For
# instance, data[2] is an array of all the values for March.
#
# The dict has various fields, including the following:
# 'begin'   : first year of data
# 'name'
# 'lat'     : (floating-point) latitude
# 'lon'     : (floating-point) longitude
# 'elevs'
# 'elevg'
# 'pop'
# 'ipop'
# 'topo'
# 'stveg',
# 'stloc'
# 'iloc'
# 'airstn'
# 'itowndis'
# 'grveg'
# 
# The BDB files contain StationStrings as serialized by the
# serialize() function.  When manipulating data for a station from a
# BDB file, it is converted to a StationString using the __init__()
# method.
# 
# The data and dict items are originally created (from the
# v2.mean_comb file and metadata files) by the from_lines() and
# v2_fill_dbm() functions.
#
# StationStrings are finally written out to a text file, at the end of
# step 1, as serialized by the to_text() method.
# 
# The rotation between year-major data (in the inputs and outputs
# of this step) and month-major data (everywhere inside this step) is
# done by from_lines() and to_text().

class StationString:
  def to_text(self, id):
    dict = self.dict
    data = self.data
    begin = dict['begin']
    text = ' %4d%5d%s%4s%-36s\n' % (math.floor(dict['lat']*10+0.5),
                                    math.floor(dict['lon']*10+0.5),
                                    id, dict['elevs'], dict['name'])
    lines = []
    for y in range(self.years):
      l = map(lambda m: "%5d" % math.floor(data[m][y] * 10.0 + 0.5), range(12))
      lines.append('%4d%s\n' % (y + begin, ''.join(l)))
    text = text + ''.join(lines)
    return text

  def __init__(self, sp):
    sep = sp.index('\n')

    self.dict = {}
    l = sp[:sep].split('\t')
    while l:
        k = l[0]
        v = l[1]
        l = l[2:]
        if k in ['begin', 'ipop']:
          v = int(v)
        elif k in ['lat', 'lon']:
          v = float(v)
        self.dict[k] = v

    l = sp[sep+1:].split(' ')
    years = int(len(l)/12)
    assert len(l) == years * 12
    l = map(lambda s: float(s) * 0.1, l)
    self.data = []
    while(l):
      self.data.append(array.array('d',l[:years]))
      l = l[years:]
    self.years = years

def serialize(dict, data):
    dict_l = []
    for (k,v) in dict.items():
      dict_l.append(k)
      dict_l.append(str(v))
    data_l = []
    for m in range(12):
      data_l.extend(map(lambda d: "%d" % math.floor(d*10.0 + 0.5), data[m]))
    return   '\t'.join(dict_l) + '\n' + ' '.join(data_l)

def from_lines(lines):
    begin = 9999
    end = 0
    for line in lines:
      year = int(line[12:16])
      if year < begin:
        begin = year
      if year > end:
        end = year
    years = end - begin + 1
    l = [999.9]*years
    data = map(lambda x: array.array('d',l), range(12))
    for line in lines:
      year = int(line[12:16])
      index = year-begin
      for m in range(12):
        datum = int(line[16+5*m:21+5*m])
        if datum == -9999:
          datum = 9999
        data[m][index] = datum * 0.1
    return (data, begin)

# Derived from the monthlydata.c python extension, which defined a lot
# of methods which are not called anywhere in GISTEMP.
# 
# The only thing it was used for is calculating a mean and set of
# annual anomalies.
# 
# However, the particular algorithm for computing these (via monthly
# and then seasonal means and anomalies, with particular thresholds
# for missing values at each step) must be maintained for bit-for-bit
# compatibility with GISTEMP.

def monthly_annual(data, bad):
        years = len(data[0])
        # compute monthly means and anomalies   
        monthly_means = array.array('d',[bad]*12)
        monthly_anoms = []
        for m in range(12):
            row = data[m]
            sum = 0.0
            count = 0
            for n in range(years):
                datum = row[n]
                if not FEQUALS(datum, bad):
                    sum += datum
                    count += 1
            if count > 0:
                monthly_means[m] = sum/count
            mean = monthly_means[m]
            monthly_anoms_row = array.array('d',[bad]*years)
            for n in range(years):
                datum = row[n]
                if (not FEQUALS(mean, bad) and not FEQUALS(datum, bad)):
                    monthly_anoms_row[n] = datum - mean
            monthly_anoms.append(monthly_anoms_row)

        seasonal_means = array.array('d',[bad, bad, bad, bad])
        seasonal_anoms = []
        # compute seasonal anomalies
        for s in range(4):
            if s == 0:
                months = [11, 0, 1]
            else:
                months = range(s*3-1, s*3+2)
            sum = 0.0
            count = 0
            for i in range(3):
                m = months[i]
                monthly_mean = monthly_means[m]
                if not FEQUALS(monthly_mean, bad):
                    sum += monthly_mean
                    count += 1
            if count > 1:
                seasonal_means[s] = sum / count
            seasonal_anoms_row = array.array('d', [bad] * years)
            for n in range(years):
                sum = 0.0
                count = 0
                for i in range(3):
                    m = months[i]
                    if m == 11:
                      if n == 0:
                        continue
                      else:
                        monthly_anom = monthly_anoms[m][n-1]
                    else:
                        monthly_anom = monthly_anoms[m][n]
                    if not FEQUALS(monthly_anom, bad):
                        sum += monthly_anom
                        count += 1
                if count > 1:
                    seasonal_anoms_row[n] = sum / count
            seasonal_anoms.append(seasonal_anoms_row)
        
        # annual mean and anomalies
        sum = 0.0
        count = 0
        for s in range(4):
            seasonal_mean = seasonal_means[s]
            if not FEQUALS(seasonal_mean, bad):
                sum += seasonal_mean
                count += 1
        if count > 2:
                annual_mean = sum / count
        else:
                annual_mean = bad
        annual_anoms = array.array('d',[bad]*years)
        for n in range(years):
            sum = 0.0
            count = 0
            for s in range(4):
                anom = seasonal_anoms[s][n]
                if not FEQUALS(anom, bad):
                    sum += anom
                    count += 1
            if count > 2:
                annual_anoms[n] = sum / count
        return (annual_mean, annual_anoms)

# From v2_to_bdb.py

def v2_fill_dbm(f, dbm, info, sources):
    ids = []
    count = 0
    line = f.readline()
    while line:
        lines = [line]
        id = line[:12]
        ids.append(id)
        last_id = id
        line = f.readline()
        while line:
            id = line[:12]
            if id != last_id:
                break
            lines.append(line)
            line = f.readline()
        data, begin = from_lines(lines)
        st_id = last_id[:-1]
        dict = info[st_id]
        dict['begin'] = begin
        dict['source'] = sources.get(last_id, 'UNKNOWN')
        mystring = serialize(dict, data)
        dbm[last_id] = mystring
        count = inc_count(count)
    print count
    dbm['IDS'] = ' '.join(ids)
    dbm['IBAD'] = '9999'
    dbm.close()

def v2_get_sources():
    sources = {}
    for source in ['MCDW', 'USHCN', 'SUMOFDAY']:
        for line in open('input/%s.tbl' % source.lower()):
            _, id, rec_no = line.split()
            sources[id + rec_no] = source
    return sources

def v2_get_info():
    keys = ('id', 'name', 'lat', 'lon', 'elevs',
            'elevg', 'pop', 'ipop', 'topo', 'stveg',
            'stloc', 'iloc', 'airstn', 'itowndis', 'grveg')
    info = {}
    for row in open('input/v2.inv', 'r'):
        values = struct.unpack('11sx30sx6sx7sx4s5sc5s2s2s2s2sc2s16sx', row[:101])
        # Would use struct.unpack_from if we knew we had Python 2.5
        hash = dict(zip(keys, values))
        id = hash['id']
        info[id] = hash
        del hash['id']
    return info

def v2_to_bdb(infile_name):
    bdb_name = 'work/' + infile_name + '.bdb'
    f = open('work/' + infile_name, 'r')
    print "reading " + infile_name
    info = v2_get_info()
    sources = v2_get_sources()
    dbm = bsddb.hashopen(bdb_name, 'n')
    print "writing " + bdb_name
    v2_fill_dbm(f, dbm, info, sources)

# from comb_records.py
# 
# Note that unlike the other phases, some of this code was shared with
# comb_pieces.py

BAD = '???'
MIN_OVERLAP = 4
comb_log = None

def average(new_sums, new_counts, new_data, years):
    for m in range(12):
        for n in range(years):
            count = new_counts[m][n]
            if count != 0:
                new_data[m][n] = new_sums[m][n] / count

def final_average(new_sums, new_wgts, new_data, years, begin):
    min, max = 9999, -9999
    for m in range(12):
        mon_min, mon_max = 9999, -9999
        wgts_row = new_wgts[m]
        for n in range(years):
            wgt = wgts_row[n]
            if wgt == 0:
                continue
            if mon_min > n:
                mon_min = n
            if mon_max < n:
                mon_max = n
        if min > mon_min:
            min = mon_min
        if max < mon_max:
            max = mon_max
    if min == 0 and max == years - 1:
        average(new_sums, new_wgts, new_data, years)
        return begin
    years = max - min + 1
    begin = begin + min
    end = begin + years - 1
    for m in range(12):
        new_sums[m] = new_sums[m][min: max + 1]
        new_wgts[m] = new_wgts[m][min: max + 1]
        new_data[m] = [BAD] * years
    average(new_sums, new_wgts, new_data, years)
    return begin

def add(new_sums, new_wgts, diff, begin, record):
    rec_begin = record['dict']['begin']
    rec_years = record['years']
    rec_data = record['data']
    for m in range(12):
        sums_row = new_sums[m]
        wgts_row = new_wgts[m]
        data_row = rec_data[m]
        for n in range(rec_years):
            datum = data_row[n]
            if abs(datum - BAD) < 0.1:
                continue
            year = n + rec_begin
            index = year - begin
            sums_row[index] = sums_row[index] + datum - diff
            wgts_row[index] = wgts_row[index] + 1

def get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records):
    end = begin + years - 1
    average(new_sums, new_wgts, new_data, years)
    ann_mean, ann_anoms = monthly_annual(new_data, BAD)
    overlap = 0
    length = 0
    for rec_id, record in records.items():
        rec_ann_anoms = record['ann_anoms']
        rec_ann_mean = record['ann_mean']
        rec_years = record['years']
        rec_begin = record['dict']['begin']
        sum = wgt = 0
        for n in range(rec_years):
            rec_anom = rec_ann_anoms[n]
            if abs(rec_anom - BAD) < 0.1:
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if abs(anom - BAD) < 0.1:
                continue
            wgt = wgt + 1
            sum = sum + (rec_ann_mean + rec_anom) - (ann_mean + anom)
        if wgt < MIN_OVERLAP:
            continue
        if wgt < overlap:
            continue
        overlap = wgt
        diff = sum / wgt
        best_id = rec_id
        best_record = record
    if overlap < MIN_OVERLAP:
        return 0, 0, BAD
    return best_record, best_id, diff

def combine(*args):
    new_sums, new_wgts, new_data, begin, years, records = args
    record, rec_id, diff = apply(get_longest_overlap, args)
    if abs(diff - BAD) < 0.1:
        comb_log.write("\tno other records okay\n")
        return 0
    del records[rec_id]
    add(new_sums, new_wgts, diff, begin, record)
    rec_begin = record['dict']['begin']
    comb_log.write("\t %s %d %d %f\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1, diff))
    return 1

def get_best(records):
    ranks = {'MCDW': 4, 'USHCN': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    rids = records.keys()
    rids.sort()
    for rec_id in rids:
        record = records[rec_id]
        source = record['dict']['source']
        length = record['length']
        rank = ranks[source]
        if rank > best:
            best = rank
            best_rec = record
            best_id = rec_id
        elif length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    if best > 1:
        return best_rec, best_id
    return longest_rec, longest_id

def get_records(old_db, rec_ids):
    res = {}
    min, max = 9999, -9999
    for rec_id in rec_ids:
        s = old_db[rec_id]
        st = StationString(s)
        assert st.years == len(st.data[0])
        begin = st.dict['begin']
        end = begin + st.years - 1
        if min > begin:
            min = begin
        if max < end:
            max = end
        ann_mean, ann_anoms = monthly_annual(st.data, BAD)
        length = 0
        for anom in ann_anoms:
            if anom != BAD:
                length = length + 1
        res[rec_id] = {'dict': st.dict, 'data': st.data,
                       'string': s, 'years': st.years,
                       'length': length, 'ann_anoms': ann_anoms,
                       'ann_mean': ann_mean}
    return res, min, max

def records_get_new_data(record, begin, years):
    sums, wgts, data = [None] * 12, [None] * 12, [None] *12
    rec_data = record['data']
    rec_begin, rec_years = record['dict']['begin'], record['years']
    rec_end = rec_begin + rec_years - 1
    for m in range(12):
        sums_row, wgts_row, data[m] = [0.0] * years, [0] * years, [BAD] * years
        rec_row = rec_data[m]
        for n in range(rec_years):
            datum = rec_row[n]
            if abs(datum - BAD) < 0.1:
                continue
            index = n + rec_begin - begin
            sums_row[index] = datum
            wgts_row[index] = 1
        sums[m] = sums_row
        wgts[m] = wgts_row
    return sums, wgts, data

def records_fill_new_db(ids, rec_id_dict, old_db, new_db):
    count = 0
    new_ids = []
    for id in ids:
        rec_ids = rec_id_dict[id]
        comb_log.write('%s\n' % id)
        while 1:
            records, begin, end = get_records(old_db, rec_ids)
            if len(records) == 1:
                rec_id, record = records.items()[0]
                new_db[rec_id] = record['string']
                new_ids.append(rec_id)
                count = inc_count(count)
                break
            years = end - begin + 1
            record, rec_id = get_best(records)
            rec_dict = record['dict']
            source = rec_dict['source']
            del records[rec_id]
            new_sums, new_wgts, new_data = records_get_new_data(record, begin, years)
            new_dict = rec_dict.copy()
            new_dict['begin'] = begin
            comb_log.write ("\t%s %s %s -- %s\n" % (rec_id, begin, end,source))
            ok_flag = 1
            while ok_flag and records.keys():
                ok_flag = combine(new_sums, new_wgts, new_data,
                                  begin, years, records)
            begin = final_average(new_sums, new_wgts, new_data, years, begin)
            new_dict['begin'] = begin
            s = serialize(new_dict, new_data)
            new_db[rec_id] = s
            new_ids.append(rec_id)
            count = inc_count(count)
            rec_ids = records.keys()
            if not rec_ids:
                break
    print count
    new_db['IDS'] = ' '.join(new_ids)

def get_ids(old_db):
    rec_ids = old_db['IDS'].split()
    rec_id_dict = {}
    for rec_id in rec_ids:
        id = rec_id[:-1]
        if not rec_id_dict.has_key(id):
            rec_id_dict[id] = []
        rec_id_dict[id].append(rec_id)
    ids = rec_id_dict.keys()
    ids.sort()
    return ids, rec_id_dict

def comb_records(stemname):
    global comb_log
    comb_log = open('log/comb.log','w')
    old_bdb_name = 'work/' + stemname + '.bdb'
    new_bdb_name = 'work/' + stemname + '.combined.bdb'
    print "opening db file '%s'" % old_bdb_name
    old_db = bsddb.hashopen(old_bdb_name, 'r')
    print "creating db file '%s'" % new_bdb_name
    new_db = bsddb.hashopen(new_bdb_name, 'n')
    global BAD
    BAD = int(old_db["IBAD"]) / 10.0
    BAD = 0.1 * int(BAD * 10.0)
    new_db['IBAD'] = old_db['IBAD']
    ids, rec_id_dict = get_ids(old_db)
    records_fill_new_db(ids, rec_id_dict, old_db, new_db)
    comb_log.close()

# From comb_pieces.py.

BAD = '???'
MIN_MID_YEARS = 5            # MIN_MID_YEARS years closest to ymid
BUCKET_RADIUS = 10
BIG_X = 1.                     # "we may then use a fraction X of ..."
GET_ALL = 1
VERBOSE = 0
pieces_log = None

HELENA_IN_NAME = 'combine_pieces_helena.in'

def sigma(list, bad):
    sos = sum = count = 0
    for x in list:
        if x == bad: continue
        sos = sos + x * x
        sum = sum + x
        count = count + 1
    if count == 0:
        return bad
    else:
        mean = sum/count
        return math.sqrt(sos / count - mean * mean)

def pieces_get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records):
    end = begin + years - 1
    average(new_sums, new_wgts, new_data, years)
    ann_mean, ann_anoms = monthly_annual(new_data, BAD)
    overlap = 0
    length = 0
    for rec_id, record in records.items():
        rec_ann_anoms = record['ann_anoms']
        rec_years = record['years']
        rec_begin = record['dict']['begin']
        sum = wgt = 0
        for n in range(rec_years):
            rec_anom = rec_ann_anoms[n]
            if abs(rec_anom - BAD) < 0.1:
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if abs(anom - BAD) < 0.1:
                continue
            wgt = wgt + 1
            sum = sum + rec_anom - anom
        if wgt < overlap:
            continue
        overlap = wgt
        # diff = sum / wgt
        best_id = rec_id
        best_record = record
    ignore_me = "ignore_me"
    return best_record, best_id, ignore_me

def get_actual_endpoints(new_wgts, begin, years):
    tmp_begin = 9999
    tmp_end = 0
    for m in range(12):
        wgts_row = new_wgts[m]
        for n in range(years):
            if wgts_row[n] == 0:
                continue
            year = n + begin
            if year < tmp_begin:
                tmp_begin = year
            if tmp_end < year:
                tmp_end = year
    return tmp_begin, tmp_end

def find_quintuples(new_sums, new_wgts, new_data,
                    begin, years, record, rec_begin,
                    new_id=None, rec_id=None):
    end = begin + years - 1

    rec_begin = record['dict']['begin']
    rec_end = rec_begin + record['years'] - 1

    actual_begin, actual_end = get_actual_endpoints(new_wgts, begin, years)

    max_begin = max(actual_begin, rec_begin)
    min_end = min(actual_end, rec_end)
    middle_year = int(.5 * (max_begin + min_end) + 0.5)
    pieces_log.write("max begin: %s\tmin end: %s\n" % (max_begin, min_end))

    average(new_sums, new_wgts, new_data, years)
    ann_mean, sublist1 = monthly_annual(new_data, BAD)
    ann_std_dev = sigma(sublist1, BAD)
    pieces_log.write("ann_std_dev = %s\n" % ann_std_dev)
    offset1 = (middle_year - begin)
    len1 = len(sublist1)

    sublist2 = record['ann_anoms']
    rec_ann_mean = record['ann_mean']
    offset2 = (middle_year - rec_begin)
    len2 = len(sublist2)

    ov_success = 0
    okay_flag = 0
    for rad in range(1, BUCKET_RADIUS + 1):
        if rad > 1:
            if VERBOSE:
                pieces_log.write("incrementing bucket radius: %s -> %s\n" % (rad - 1, rad))
        count1 = sum1 = 0
        count2 = sum2 = 0
        for i in range(0, rad + 1):
            for sign in [-1, 1]:
                if sign == 1 and i == 0:
                    continue
                index1 = i * sign + offset1
                index2 = i * sign + offset2
                if index1 < 0 or index1 >= len1:
                    anom1 = BAD
                else:
                    anom1 = sublist1[index1]
                if index2 < 0 or index2 >= len2:
                    anom2 = BAD
                else:
                    anom2 = sublist2[index2]
                if VERBOSE:
                    rec_mean = anom2 + rec_ann_mean
                    mean = anom1 + ann_mean
                    pieces_log.write ("*%s\t%s\t%s\t%s\t%s*\n" % 
                                      (middle_year, index1 + begin,
                                       index2 + rec_begin,
                                       mean, rec_mean))
                if anom1 != BAD and (GET_ALL or count1 < MIN_MID_YEARS):
                    sum1 = sum1 + anom1 + ann_mean
                    count1 = count1 + 1
                if anom2 != BAD and (GET_ALL or count2 < MIN_MID_YEARS):
                    sum2 = sum2 + anom2 + rec_ann_mean
                    count2 = count2 + 1
        if count1 >= MIN_MID_YEARS and count2 >= MIN_MID_YEARS:
            pieces_log.write("overlap success: %s %s\n" % (new_id, rec_id))
            ov_success = 1
            avg1 = sum1 / float(count1)
            avg2 = sum2 / float(count2)
            diff = abs(avg1 - avg2)
            pieces_log.write("diff = %s\n" % diff)
            if diff < BIG_X * ann_std_dev:
                okay_flag = 1
                pieces_log.write("combination success: %s %s\n" % (new_id, rec_id))
            else:
                pieces_log.write("combination failure: %s %s\n" % (new_id, rec_id))
            break
    if not ov_success:
        pieces_log.write("overlap failure: %s %s\n" % (new_id, rec_id))
    pieces_log.write("counts: %s\n" % ((count1, count2),))
    return okay_flag

def pieces_combine(*args, **kwargs):
    new_sums, new_wgts, new_data, begin, years, records = args
    new_id = kwargs['new_id']
    end = begin + years - 1

    record, rec_id, ignore_me = apply(pieces_get_longest_overlap, args)
    rec_begin = record['dict']['begin']
    rec_end = rec_begin + record['years'] - 1

    pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))

    is_okay = find_quintuples(new_sums, new_wgts, new_data,
                              begin, years, record, rec_begin,
                              new_id=new_id, rec_id=rec_id)

    if not is_okay:
        pieces_log.write("\t***no other pieces okay***\n")
        return 0
    else:
        del records[rec_id]
        add(new_sums, new_wgts, 0.0, begin, record)
        pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))
        return 1

def get_longest(records):
    longest = 0
    for rec_id, record in records.items():
        length = record['length']
        if length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    return longest_rec, longest_id

def pieces_fill_new_db(ids, rec_id_dict, old_db, new_db, helena_ds):
    new_ids = []
    count = 0
    for id in ids:
        rec_ids = rec_id_dict[id]
        pieces_log.write('%s\n' % id)
        while 1:
            records, begin, end = get_records(old_db, rec_ids)

            if len(records) == 1:
                rec_id, record = records.items()[0]
                new_db[rec_id] = record['string']
                count = inc_count(count)
                new_ids.append(rec_id)
                break

            if helena_ds.has_key(id):
                id1, id2, this_year, month, summand = helena_ds[id]
                dict = records[id1]['dict']
                data = records[id1]['data']
                begin = dict['begin']
                years = records[id1]['years']
                for n in range(years):
                    year = n + begin
                    if year > this_year:
                        break
                    for m in range(12):
                        if year == this_year and m > month:
                            break
                        datum = data[m][n]
                        if datum == BAD:
                            continue
                        data[m][n] = datum + summand
                del helena_ds[id]
                ann_mean, ann_anoms = monthly_annual(data, BAD)
                records[id1]['ann_anoms'] = ann_anoms
                records[id1]['ann_mean'] = ann_mean

            record, rec_id = get_longest(records)
            years = end - begin + 1
            rec_dict = record['dict']
            source = rec_dict['source']
            del records[rec_id]
            new_sums, new_wgts, new_data = records_get_new_data(record, begin, years)
            new_dict = rec_dict.copy()
            new_dict['begin'] = begin
            pieces_log.write("\t%s %s %s -- %s\n" % (rec_id, begin, begin + years -1,source))
            ok_flag = 1
            while ok_flag and records.keys():
                ok_flag = pieces_combine(new_sums, new_wgts, new_data,
                                  begin, years, records, new_id=rec_id)
            begin = final_average(new_sums, new_wgts, new_data, years, begin)
            new_dict['begin'] = begin
            s = serialize(new_dict, new_data)
            new_db[rec_id] = s
            count = inc_count(count)
            new_ids.append(rec_id)
            rec_ids = records.keys()
            if not rec_ids:
                break
    print count
    new_db['IDS'] = ' '.join(new_ids)

def comb_pieces(stemname):
    global pieces_log
    pieces_log = open('log/pieces.log','w')
    old_name = 'work/' + stemname + '.bdb'
    new_name = 'work/' + stemname + '.pieces.bdb'
    print "opening db file '%s'" % old_name
    old_db = bsddb.hashopen(old_name, 'r')
    print "creating db file '%s'" % new_name
    new_db = bsddb.hashopen(new_name, 'n')
    global BAD
    BAD = int(old_db["IBAD"]) / 10.0
    BAD = 0.1 * int(BAD * 10.0)
    comb_records.BAD = BAD
    new_db['IBAD'] = old_db['IBAD']
    ids, rec_id_dict = get_ids(old_db)

    try:
        helena_in = open('config/' + HELENA_IN_NAME, 'r')
        helena_ls = helena_in.readlines()
        helena_in.close()
        print "successfully read '%s'" % HELENA_IN_NAME
    except IOError:
        helena_ls = []
    helena_ds = {}
    for x in helena_ls:
        id1, id2, year, month, summand = x.split()
        year = int(year)
        month = int(month)
        summand = float(summand)
        id = id1[:-1]
        helena_ds[id] = id1, id2, year, month, summand

    pieces_fill_new_db(ids, rec_id_dict, old_db, new_db, helena_ds)
    pieces_log.close()

# From drop_strange.py
#
# Removes records as specified in config/Ts.strange.RSU.list.IN. Each
# line in that file begins with a 12-digit station ID - actually the
# tuple (country-code, WMO station, modifier, duplicate) - and ends
# with either yyyy/mm, specifying a month datum to omit or with
# xxxx-yyyy, specifying years to omit.  xxxx can be 0, meaning from
# the beginning. yyyy can be 9999, meaning to the end.
#
# I have looked through this code; it's very straightforward.  It
# could be very much more concise, and quicker, but it's not too bad.
# 
# Nick Barnes 2008-09-18.

BAD_YEAR = 9999
BAD = None

TS_CHANGES_NAME = "Ts.strange.RSU.list.IN"

def drop_get_new_data(new_begin, new_years, old_begin, old_years, old_data,
                 year1, year2):
    assert new_begin >= old_begin
    assert new_begin + new_years <= old_begin + old_years

    new_data = [None] * 12
    for m in range(12):
        new_monthly = []
        old_monthly = old_data[m]
        for n in range(new_years):
            year = n + new_begin
            if year1 <= year <= year2:
                datum = BAD
            else:
                datum = old_monthly[year - old_begin]
            new_monthly.append(datum)
        new_data[m] = new_monthly
    return new_data

def remove_middle(dict, data, years, begin, year1, year2):
    end = begin + years - 1
    if year1 <= begin:
        new_begin = year2 + 1
    else:
        new_begin = begin
    if year2 > end:
        year2 = end
    new_end = end
    new_years = new_end - new_begin + 1

    new_data = drop_get_new_data(new_begin, new_years, begin, years, data,
                            year1, year2)
    new_dict = dict.copy()
    new_dict['begin'] = new_begin

    return serialize(new_dict, new_data)

def implement_changes(changes_dict, old_bdb, new_bdb):
    new_ids = []
    count = 0
    for id, el in changes_dict.items():
        s = None
        for a, x, y in el:
            if x == 0 and y == BAD_YEAR:
                continue
            if s is None:
                s = old_bdb[id]
            st = StationString(s)
            begin = st.dict['begin']
            if a == "month":
                year = x
                month = y
                st.data[month - 1][year - begin] = BAD
                new_s = serialize(st.dict, st.data)
            else:
                assert a == "period"
                new_s = remove_middle(st.dict, st.data, st.years, begin, x, y)
            s = new_s
        if s is not None:
            new_bdb[id] = s
            count = inc_count(count)
            new_ids.append(id)

    old_ids = old_bdb['IDS'].split()
    for id in old_ids:
        if changes_dict.has_key(id):
            continue
        new_bdb[id] = old_bdb[id]
        count = inc_count(count)
        new_ids.append(id)

    print count
    new_ids.sort()
    new_bdb['IDS'] = ' '.join(new_ids)
    new_bdb['IBAD'] = old_bdb['IBAD']

def get_changes_dict():
    dict = {}

    f = open('config/' + TS_CHANGES_NAME, 'r')
    print "reading", TS_CHANGES_NAME

    for line in f.readlines():
        split_line = line.split()
        id = split_line[0]
        try:
            year1, year2 = map(int, split_line[-1].split("-"))
            val = ("period", year1, year2)
        except ValueError:
            year, month = map(int, split_line[-1].split("/"))
            val = ("month", year, month)
        try:
            dict[id].append(val)
        except KeyError:
            dict[id] = [val]
    return dict

def drop_strange(in_name):
    changes_dict = get_changes_dict()

    OLD_BDB_NAME = 'work/' + in_name + ".bdb"
    NEW_BDB_NAME = 'work/' + in_name + ".strange.bdb"
    old_bdb = bsddb.hashopen(OLD_BDB_NAME, 'r')
    print "reading", OLD_BDB_NAME

    global BAD
    BAD = int(old_bdb["IBAD"]) / 10.0

    new_bdb = bsddb.hashopen(NEW_BDB_NAME, "c" or "n")
    print "writing", NEW_BDB_NAME

    implement_changes(changes_dict, old_bdb, new_bdb)

# From alter_discont.py
#
# Modifies records as specified in config/Ts.discont.RS.alter.IN.
# Each line in that file has a 12 digit station ID, a month, a year,
# and a floating-point temperature delta.
# This phase adds the delta to every datum for that station prior to
# the specified month.
#
# I have looked through this code; it's tolerable. It could be very
# much more concise, and quicker, but it's not too bad.  Also it could
# be combined with other phases very easily.

# 
# Nick Barnes 2008-09-18.

BAD = None
IN_FILE_NAME = "Ts.discont.RS.alter.IN"

def alter_discont(in_name):
    IN_BDB_NAME = 'work/' + in_name + ".bdb"
    OUT_BDB_NAME = 'work/' + in_name + ".alter.bdb"
    in_bdb = bsddb.hashopen(IN_BDB_NAME, "r")
    print "reading", IN_BDB_NAME
    out_bdb = bsddb.hashopen(OUT_BDB_NAME, "c")
    print "writing", OUT_BDB_NAME

    alter_dict = get_alter_dict()

    global BAD
    BAD = int(in_bdb["IBAD"]) / 10.0
    BAD = 0.1 * int(BAD * 10.0)
    out_bdb["IBAD"] = in_bdb["IBAD"]

    alter_fill_new_bdb(in_bdb, out_bdb, alter_dict)

def get_alter_dict():
    f = open('config/' + IN_FILE_NAME)
    print "reading", IN_FILE_NAME
    dict = {}
    for line in f:
        id, month, year, num = line.split()
        dict[id] = [int(month), int(year), float(num)]
    return dict

def alter_fill_new_bdb(in_bdb, out_bdb, alter_dict):
    ids = in_bdb["IDS"].split()
    count = 0
    for id, (a_month, a_year, a_num) in alter_dict.items():
        s = in_bdb[id]
        st = StationString(s)
        data = st.data
        years = len(data[0])
        begin = st.dict["begin"]
        for year in range(begin, a_year + 1):
            index = year - begin
            for m in range(12):
                if year == a_year and m > a_month - 2:
                    continue
                datum = data[m][index]
                if datum == BAD:
                    continue
                data[m][index] = datum + a_num
        s = serialize(st.dict, data)
        out_bdb[id] = s
        count = inc_count(count)
    for id in ids:
        if alter_dict.has_key(id):
            continue
        out_bdb[id] = in_bdb[id]
        count = inc_count(count)
    print count
    out_bdb["IDS"] = ' '.join(ids)

# From bdb_to_text.py

def bdb_to_text(db_file, txt_file):
    db = bsddb.hashopen('work/' + db_file, 'r')
    print "reading", db_file
    f = open('work/' + txt_file, 'w')
    print "creating", txt_file
    ids = db['IDS'].split()
    count = 0
    for id in ids:
        count = inc_count(count)
        s = db[id]
        st = StationString(s)
        f.write(st.to_text(id))
    print count
    f.close()
    db.close()

def main():
  v2_to_bdb('v2.mean_comb')
  comb_records('v2.mean_comb')
  comb_pieces('v2.mean_comb.combined')
  drop_strange('v2.mean_comb.combined.pieces')
  alter_discont('v2.mean_comb.combined.pieces.strange')
  bdb_to_text('v2.mean_comb.combined.pieces.strange.alter.bdb', 'Ts.txt')

if __name__ == '__main__':
  main()

# Notes on STEP1 algorithm:
#
# The step is driven by do_comb_step1.sh, a ksh script
# 
# 1. v2_to_bdb.py creates v2.mean_comb.bdb containing the same data as v2.mean_comb.
# The records in each BDB are as follows:
# 
# IDS -> space-separated list of station IDs
# BAD -> '9999'
# id  -> station data (format defined by stationstring.serialize)
# 
# 2. comb_records.py creates v2.mean_comb.combined.bdb and comb.log
#
# 3. comb_pieces.py creates v2.mean_comb.combined.pieces.bdb and piece.log
#
# 4. drop_strange.py creates v2.mean_comb.combined.pieces.strange.bdb
#
# 5. alter_discont.py creates  v2.mean_comb.combined.pieces.strange.alter.bdb
#
# 6. bdb_to_text.py creates v2.mean_comb.combined.pieces.strange.alter.txt
# which is then renamed to to_next_step/Ts.txt
