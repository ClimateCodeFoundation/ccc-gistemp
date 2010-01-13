#!/usr/bin/env python
# step3.py
# 
# Python code reproducing the STEP3 part of the GISTEMP algorithm.
#
# Work in progress.
# David Jones, Ravenbrook Limited, 2008-08-06
#
# The code is derived from the Fortran GISTEMP code.
#
# Python notes:
#
# Much of the binary IO uses the struct module.  A quick review will
# help.
#
# I feel obliged to tell you that the Python expression N * list gives a
# new list that conists of list repeated N times.  This is used every
# now and then.
#
# We also make use of list slice assignment: a[2:4] = range(2)
# 
# $Id: //info.ravenbrook.com/project/ccc/master/code/step3.py#20 $

# Ravenbrook
import earth # required for radius.
# Ravenbrook
import eqarea
# Ravenbrook
import fort

# http://www.python.org/doc/2.3.5/lib/module-array.html
from array import array
# http://www.python.org/doc/2.3.5/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.3.5/lib/module-math.html
import math
# http://www.python.org/doc/2.3.5/lib/module-struct.html
import struct
# http://www.python.org/doc/2.3.5/lib/module-sys.html
import sys
# http://www.python.org/doc/2.3.5/lib/module-warnings.html
from warnings import warn

# The width of a standard word according to Python's struct module.
w = len(struct.pack('=I', 0))

class Stations:
    """A representation of a collection of station data.  Concretely a
    station is a temperature recording station situated at some point on
    the Earth; abstractly each station is a series of numeric values
    with some associated metadata, including location.

    Each instance has the following public members:
    kq - the KQ quantity from the header record.
    mavg - the MAVG flag from the header record.
    monm - maximum length of any time record.
    recsize - INFO(5) from the header record.  Probably not much use.
    yrbeg - YRBEG from the header record.  The year of first data.
    bad - the XBAD value used in the data.
    trace - INFO(8) from the header record.
    scale - the output scale derived from KQ.  Equivalent to SCL in the
        Fortran code.
    km - number of records per year, derived from MAVG.
    """

    def __init__(self, infile):
        """infile should be a sequence of file objects.  Each file
        object should be opened in binary onto a trimmed NCAR file,
        identical in format to the input files that GISTEMP's
        to.SBBXgrid.f program expects.  Normally there will be 6 zonal
        files, but this is not required; all input files are processed.
        """

        self.yrbeg = None
        self.kq = None
        self.monm = None
        self.bad = None
        self.header = []

        self.station = []

        for f in infile:
            ff = fort.File(f)

            r = ff.readline()
            self.header.append(r)
            self.title = r[9*w:]
            # First record contains various INFO items.
            # When referring to comments in to.SBBXgrid.f recall that Fortran
            # arrays are typically indexed from 1 onwards; Python from 0
            # onwards.  Therefore INFO(1) corresponds to a[0]
            a = struct.unpack('9i', r[:9*w])
            mfirst = a[0]
            mlast = a[8]
            kq = a[1]
            self.kq = self.kq or kq
            mavg = a[2]
            monm = a[3]
            self.monm = self.monm or monm
            recsize = a[4]
            yrbeg = a[5]
            # Check all YRBEGs are the same
            if self.yrbeg and yrbeg != self.yrbeg :
                warn(('File %s has a YRBEG of %d, different from the ' +
                      'YRBEG of %d on the first file.') %
                      (f.name, yrbeg, self.yrbeg))
            self.yrbeg = self.yrbeg or yrbeg
            bad = a[6]
            self.bad = self.bad or bad
            trace = a[7]

            for r in ff:
                self.station.append(Station(r, mfirst, mlast, bad=bad))
                mfirst,mlast = struct.unpack('2I', r[-2*w:])

        # This computation of scale isn't clearly correct, for example,
        # are KQ values bigger than 20 really supposed to be silently
        # accepted?  But it is what the Fortran code does.
        self.scale = 1.0 # SCL in Fortran code
        # Scale values, see to.SBBXgrid.f line 107
        # Python array has a leading initial 0.0 entry to make the
        # indexing work.
        scale = [0.0] + 3*[0.1] + 9*[0.0] + [0.1] + 7*[0.0]
        # Compute scale, to.SBBXgrid.f lines 133 and following
        if self.kq < len(scale) :
            self.scale = scale[self.kq]
        if self.scale == 0 :
            raise 'Program not ready for quantity %d' % stations.kq

        # Derive KM
        # to.SBBXgrid.f lines 140 and following
        self.km = 1
        if mavg == 6:
            self.km = 12
        if mavg == 7:
            self.km = 4

    def all(self):
        """An iterator that yields every station."""

        for st in self.station:
            yield st

    def inbox(self, lats, latn, longw, longe):
        """An iterator that yields every station within the box bounded
        by the lines of latitude lats (to the south), latn (to the
        north), and the meridians at longw (to the west), and longe (to
        the east).

        In order to accommodate boxes that overlap the meridian at -180
        it is permissible for either longw to be < -180 or for longe to
        be > +180.

        For stations exactly on the boundary the
        "lower-left" rule is used.  Stations are returned if they lie on
        the southern boundary or the western boundary except for
        corners; stations lying exactly on a corner are only returned if
        it is the south-west corner.  Note therefore that to include a
        station situated exactly at the North Pole a latn value slightly
        larger than 90 should be used.
        """

        assert lats <= latn
        assert longw <= longe

        for st in self.station:
            lat = st.lat
            lon = st.lon

            # See to.SBBXgrid.f lines 213 and following
            if lon > longe:
                lon -= 360
            elif lon < longw:
                lon += 360

            if (lats < lat < latn) and (longw < lon < longe):
                yield st
            if lats == lat and longw <= lon < longe:
                yield st
            if longw == lon and lats <= lat < latn:
                yield st

    def incircle(self, arc, lat, lon):
        """Returns all stations within a certain distance of a point of
        interest.  See the function incircle for details, the entire
        station list is used as the iterable."""

        return incircle(self.all(), arc, lat, lon)

def incircle(iterable, arc, lat, lon):
    """An iterator that filters iterable (the argument) and yields every
    station with a certain distance of the point of interest given by
    lat and lon (in degrees).

    A station is returned if the great circle arc between it and the point
    of interest is less than arc radians (using angles is slightly
    odd, but makes it independent of sphere size).

    A pair (2-tuple) (t,station) is returned where t is the
    weight assigned to the station based on its distance from
    the point.  0 <= t <= 1.  The GISTEMP algorithm is used: t
    is documented as 1-(d/arc) where d is the angle, in radians,
    made by the arc between the station and the point of interest,
    but in fact d is approximated by the chord length.

    This is essentially a filter; the stations that are returned are in
    the same order in which they appear in iterable.

    The pair is returned with the weight first so that a list of the
    pairs (for example obtained from
    list(incircle(stations.all(), 0.19, 51.5, 0))) can be sorted by
    weight with a simple call to its sort method, without requiring
    a compare or key method.

    From Python 2.4 onwards a sorted list (of pairs) can be got
    with: sorted(stations.incircle(stations.all(), 0.19, 51.5, 0)) .
    """

    # Warning: lat,lon in degrees; arc in radians!

    cosarc = math.cos(arc)
    coslat = math.cos(lat*math.pi/180)
    sinlat = math.sin(lat*math.pi/180)
    coslon = math.cos(lon*math.pi/180)
    sinlon = math.sin(lon*math.pi/180)

    for st in iterable:
        # A possible improvement in speed (which the corresponding
        # Fortran code does) would be to store the trig values of
        # the station location in the station object.
        sinlats = math.sin(st.lat*math.pi/180)
        coslats = math.cos(st.lat*math.pi/180)
        sinlons = math.sin(st.lon*math.pi/180)
        coslons = math.cos(st.lon*math.pi/180)

        # Cosine of angle subtended by arc between 2 points on a
        # unit sphere is the vector dot product.  See to.SBBXgrid.f
        # line 272.
        cosd = (sinlats*sinlat +
            coslats*coslat*(coslons*coslon + sinlons*sinlon))
        if cosd > cosarc:
            # Calculate weight based on chord length
            # See to.SBBXgrid.f line 285
            d = math.sqrt(2*(1-cosd))
            d /= arc
            yield (1-d, st)


class Station:
    """A representation of a time series of data from a station (and
    associated metadata."""

    def __init__(self, record, m0, ml, bad=9999):
        """Create a station series from its representation as a Fortran
        binary record.  m0 is the month corresponding to the first data
        item; ml is the month corresponding to the last data item.  ml
        is not used except to check that the lengh of the record is
        correct.  bad is the data value corresponding to a bad entry
        (normally stored in the file's metadata); bad is used to
        calculate the ngood member, the number of good data values."""

        # Length of record trail, the non-variable part, in bytes.
        ltrail = 15*w

        # The month where the data series begins.  This value is simply
        # taken from the input file, but conventionally the value 1
        # indicates January in the year YRBEG (which is indicated in the
        # file header).
        self.m0 = m0

        trail = record[-ltrail:]
        lat,lon,id,height = struct.unpack('2iIi', trail[:4*w])
        id = '%09d' % id
        # Latitude in degrees
        self.lat = lat * 0.1
        # Longitude in degrees
        self.lon = lon * 0.1
        # Height in metres.  Not expected to be used.
        self.height = height
        if len(record[:-ltrail]) != w*(ml-m0+1):
            warn(('Station ID %s has suspect record length. ' +
                'mfirst=%s mlast=%d record-length=%d\n') %
                (id, m0, ml, len(record)), SyntaxWarning)
        name = trail[4*w:-2*w]
        # raw name field (36 characters, including the final 6 used as
        # metadata)
        self.rawname = name
        # Some metadata is stored in the name field, *sigh*
        meta = name[-6:]
        name = name[:-6]
        # 3 digit country code
        cc = meta[3:6]
        meta = meta[0:3]
        # Prepend country code to station ID.
        id = cc + id
        # Station ID as a 12-digit string
        self.id = id

        data = record[:-ltrail]
        n = len(data)//w
        # The time series as an array.array
        self.data = array('i', struct.unpack('%di' % n, data))
        self.ngood = len(self.data) - self.data.count(bad)
        # We derive from m0 a beginning and end pair that correspond to
        # the range that the station data occupies in the AVG array used
        # when combining (the AVG array goes from the first month to the
        # last month)
        self.databeg = m0 - 1
        self.dataend = self.databeg + len(self.data)

def sort(l, cmp):
    """Sort the list l (in place) according to the comparison function
    cmp.  The comparison function, cmp(x, y), should return something
    less than 0 when x < y, 0 when x == y, something greater than 0 when
    x > y.  The sort is ascending in the following sense:  For the sorted
    list: When i < j, cmp(l[i], l[j]) <= 0.

    This sort is not stable.  In fact it is a deliberate emulation of
    the O(n**2) sort implemented in the SORT subroutine of
    to.SBBXgrid.f.  This is necessary to achieve results that are as
    close as possible to the GISS code.
    """

    # See to.SBBXgrid.f lines 605 and following

    for n in range(len(l)-1):
        nlmax = n
        for nn in range(n+1, len(l)):
            if cmp(l[nn], l[nlmax]) < 0:
                nlmax = nn
        # swap items at n and nlmax
        t = l[nlmax]
        l[nlmax] = l[n]
        l[n] = t
    return

def subbox_grid(infile,
        subbox_out, box_out, station_use_out, radius=1200,
        year_begin=1880, base_year=(1951,1980), audit=None):
    """(This is the equivalent of to.SBBXgrid.f) Convert the input files,
    infile[0] through infile[5], into gridded datasets which are output
    on the file (-like) object subbox_out and box_out.  station_use_out
    should also be a file object, it is used to record which records are
    used for which box (region).

    radius specifies a radius in kilometres, it is
    equivalent to the RCRIT value in the Fortran code.

    year_begin specifies the earliest year to include in the output.

    base_year is a pair specifying the base (reference) period in
    (first_year, last_year).  The last year is included in the base
    period.
    """

    log = sys.stdout

    # Helper functions (rather big ones).  These are separate in the
    # Fortran.

    # Equivalent of the Fortran subroutine CMBINE
    # Parameters as per Fortran except that km is bound lexically, and
    # nsm is returned directly.
    def combine(avg, wt, dnew, nf1, nl1, wt1, wtm, id, NOVRLP=20):
        """Run the GISTEMP combining algorithm.  This combines the data
        in the dnew array into the avg array.

        Each of the argument avg, wt, dnew is a linear array that is
        divided into "years" by considering each contiguous segment of
        km elements a year.  Only data for years in range(nf1, nl1) are
        considered and combined.  Note that range(nf1, nl1) includes nf1
        but excludes nl1 (and that this differs from the Fortran
        convention).
        
        Each month (or other subdivision, such as season, according to
        km) of the year is considered separately.  For the set of times
        where both avg and dnew have data the mean difference (a bias)
        is computed.  If there are fewer than NOVRLP years in common the
        data (for that month of the year) are not combined.  The bias is
        subtracted from the dnew record and it is point-wise combined
        into avg according to the station weight, wt1, and the exist
        weight for avg.

        id is an identifier used only when diagnostics are issued
        (when combining stations it expected to be the station ID).
        """

        # See to.SBBXgrid.f lines 519 and following

        # The Fortran code handles the arrays by just assuming them to
        # be 2D arrays of shape (*,KM).  Sadly Python array handling
        # just isn't that convenient, so look out for repeated uses of
        # "[k+km*n]" instead.

        nsm = 0
        missed = km
        missing = [True]*km
        for k in range(km):
            sumn = 0    # Sum of data in dnew
            sum = 0     # Sum of data in avg
            ncom = 0    # Number of years where both dnew and avg are valid
            for n in range(nf1, nl1):
                kn = k+km*n     # CSE for array index
                # Could specify that arguments are array.array and use
                # array.count(BAD) and sum, instead of this loop.
                if avg[kn] >= XBAD or dnew[kn] >= XBAD:
                    continue
                ncom += 1
                sum += avg[kn]
                sumn += dnew[kn]
            if ncom < NOVRLP:
                continue
            biask = float(sum-sumn)/ncom
            # Find mean bias
            wtmnew = wtm[k]+wt1
            # Note: bias is lexically captured
            bias[k] = float(wtm[k]*bias[k]+wt1*biask)/wtmnew
            wtm[k]=wtmnew
            # Update period of valid data, averages and weights
            for n in range(nf1, nl1):
                kn = k+km*n     # CSE for array index
                if dnew[kn] >= XBAD:
                    continue
                wtnew = wt[kn] + wt1
                avg[kn] = float(wt[kn]*avg[kn] + wt1*(dnew[kn]+biask))/wtnew
                wt[kn] = wtnew
                nsm += 1
            missed -= 1
            missing[k] = False
        if False and missed > 0:
            print "Unused data - ID/SUBBOX,WT", id, wt1, missing
        return nsm

    # Equivalent to Fortran subroutine TAVG.
    # See to.SBBXgrid.f lines 563 and following.
    def tavg(data, nyrs, base, limit, nr, nc, deflt=0.0):
        """tavg computes the time averages (separately for each calendar
        month if km=12) over the base period (year base to limit) and
        saves them in bias. In case of no data, the average is set to
        deflt if nr=0 or computed over the whole period if nr>0.

        Similarly to combine() data is treated as a linear array divided
        into years by considering contiguous chunks of km elements.

        Note: the Python convention for base and limit is used, the base
        period consists of the years starting at base and running up to,
        but including, the year limit.
        """

        missed = km
        len = km*[0]    # Warning: shadows builtin "len"
        for k in range(km):
            bias[k] = deflt
            sum = 0.0
            m = 0
            for n in range(base, limit):
                kn = k+km*n     # CSE for array index
                if data[kn] >= XBAD:
                    continue
                m += 1
                sum += data[kn]
            len[k] = m
            if m == 0:
                continue
            bias[k] = float(sum)/float(m)
            missed -= 1
        if nr*missed == 0:
            return
        # Base period is data free (for at least one month); use bias
        # with respect to whole series.
        for k in range(km):
            if len[k] > 0:
                continue
            print "No data in base period - MONTH,NR,NC", k, nr, nc
            sum = 0.0
            m = 0
            for n in range(nyrs):
                kn = k+km*n     # CSE for array index
                if data[kn] >= XBAD:
                    continue
                m += 1
                sum += data[kn]
            if m == 0:
                continue
            bias[k] = float(sum)/float(m)
        return

    # In Python you can't write to lexically scoped variables in an
    # output scope.  As a hack workaround we make the buffer a list and
    # use output_line_buffer[0]
    output_line_buffer = [None]

    def buffer_line_out(out, line, count=None):
        """Buffers a line of output so that the output file can be
        written in trimmed format.  See header comment of trimSBBX.f for
        details.  out is the (Fortran binary) file on which to output;
        line is the record to output except for the first word; count is
        the number of data in this record, it is actually output at the
        beginning of the _previous_ record.
        """

        # The assertion only fire after flush_out has been called (see
        # below); once flush_out has been called, no-one should be
        # calling this function, so this assertion provides a check that
        # no-one is.
        assert output_line_buffer[0] is not IOError

        # Only the first record output (the header row) is allowed to
        # not specify a count.
        assert count or not output_line_buffer[0]
        if output_line_buffer[0] is not None:
            out.writeline(struct.pack('>i', count) +
                          output_line_buffer[0])
        output_line_buffer[0] = line

    def flush_out(out):
        """Flush the last line of the output file."""

        out.writeline(struct.pack('>i',1) + output_line_buffer[0])
        # inhibit further output, see buffer_line_out
        output_line_buffer[0] = IOError
    
    def swrite(out, dmin):
        """Write binary record onto (the binary fortran file opened for
        writing) out."""

        # The Fortran code uses NSTNS which is actually an alias, via a
        # COMMON/DIAG/ declaration, to NSTCMB in the main program.  Gak.
        # See to.SBBXgrid.f lines 103 and 634
        record = struct.pack('>4i', *latlon)
        record += struct.pack('>2if', nstcmb, nstmns, dmin)
        record += struct.pack('>%df' % monm, *avg[-monm:])
        buffer_line_out(out, record, monm)

    # trimmed input files will be assumed.

    # Parameters.  These are a mixture of computed constants, and
    # parameters that influence the behaviour of the algorithm.  Many of
    # the parameters are tempting to change, but some of the code might
    # accidentally rely on particular values.

    assert len(infile) == 6
    assert radius > 0

    radius = float(radius)

    # Fortran version of output file.
    # must be big-endian because it's a GISTEMP data product.
    subboxf = fort.File(subbox_out, bos='>')

    # number of boxes
    nbox = 80
    # number of subboxes within each box
    nsubbox = 100

    # Much of the Fortran code assumes that various "parameters" have
    # particular fixed values (probably accidentally).  I don't trust the
    # Python code to avoid similar assumptions.  So assert the
    # "parameter" values here.  Just in case anyone tries changing them.
    # Note to people reading comment becuse the assert fired:  Please
    # don't assume that the code will "just work" when you change one of
    # the parameter values.  It's supposed to, but it might not.
    assert nbox == 80
    assert nsubbox == 100

    # area of subbox in squared kilometres
    # Recall area of sphere is 4*pi*(r**2)
    # Should equivalent to to.SBBXgrid.f line 113 (but with a different
    # precision)
    km2persubbox = (4*math.pi*earth.radius**2) / (nbox * nsubbox)

    # Critical radius as an angle of arc
    arc = radius / earth.radius

    stations = Stations(infile)

    # Output series cannot begin earlier than input series.
    # See to.SBBXgrid.f line 151
    if year_begin < stations.yrbeg:
        year_begin = stations.yrbeg

    BAD = stations.bad
    XBAD = float(BAD)
    km = stations.km
    # Number of months to _output_
    monm = (stations.yrbeg + (stations.monm//km) - year_begin) * km

    # Output header
    # See to.SBBXgrid.f lines 160 and following
    info = list(struct.unpack('8I', stations.header[0][:8*w]))
    info[0] = 1
    info[3] = monm
    info[4] = info[3]+7
    info[5] = year_begin
    # The edits to the title are done with it as a list of characters.
    # It's just easier that way, as slice assignment can be used.
    title = [' ']*80
    title[:len(stations.title)] = stations.title
    title[20:40] = 'ANOM (C)  CR     KM '
    # Equivalent to to.SBBXgrid.f lines 172 and 115 to 119
    title[33:37] = '%4d' % int(radius)
    if stations.kq == 2:
        title[25:29] = '(mm)'
    title[46:58] = '%4d-present' % year_begin
    # Convert title back to traditional string.
    title = ''.join(title)
    print title
    buffer_line_out(subboxf, struct.pack('>7I', *info[1:]) + title)

    regions = list(eqarea.gridsub())
    if audit:
        # Restrict processing to region containing single gridbox of
        # interest.
        regions=[regions[audit//100]]

    for region in regions:
        box = region[0]
        # Extend box according to the algorithm used in subroutine
        # GRIDEA.  See to.SBBXgrid.f lines 484 and following.
        # (that is, by half a box east and west and by a sufficient
        # amount north and south to incorporate radius).
        extent = [0]*4  # Just a dummy list of the right size
        extent[2] = box[2] - 0.5*(box[3]-box[2])
        extent[3] = box[3] + 0.5*(box[3]-box[2])
        if box[0] <= -90 or box[1] >= 90:
            # polar
            extent[2] = -180.0
            extent[3] = +180.0
        # In the GISS Fortran code the box is extended north and south
        # by DDLAT degrees.  Note the computation of DDLAT (on line 463
        # of to.SBBXgrid.f), it is the angle (in degrees) subtended by
        # an arc of length RCRIT (on a spherical earth's surface).
        # We already have that in radians in the variable arc.
        ddlat = arc*180/math.pi
        extent[0] = box[0] - ddlat
        extent[1] = box[1] + ddlat

        regionstations = list(stations.inbox(*extent))
        if False and audit:
            print 'UNSORTED LIST'
            for station in regionstations:
                print station.id
        # Descending sort by number of good records
        # Sadly we cannot use Python's sort method here, we must use an
        # emulation of the GISS Fortran.
        sort(regionstations, lambda x,y: y.ngood-x.ngood)
        if True and audit:
            print 'SORTED LIST'
            for station in regionstations:
                print station.id

        subboxes = list(region[1])
        if audit:
            # Select single subbox for audit
            subboxes = [subboxes[audit % 100]]

        # Used to generate the "subbox at" rows in the log.
        lastcentre = (None, None)
        for subbox in subboxes:
            # Convert latitude longitude to integer 100ths for swrite.
            latlon = map(lambda x:int(round(100*x)), subbox)
            # Select and weight stations
            centre = eqarea.centre(subbox)
            if centre[0] != lastcentre[0]:
                log.write("\nsubbox at %+05.1f" % centre[0])
            log.write('%+06.1f' % centre[1])
            log.flush()
            lastcentre = centre
            # Of possible stations in this region, filter for those with
            # radius of subbox centre.  Note that it is important that
            # the ordering within the regionstations list is retained
            # (in order to match the GISS code).
            # station is a list (wt, station) pairs:
            station = list(incircle(regionstations, arc, *centre))
            # Split list of pairs into pair of lists
            wti = map(lambda x: x[0], station)
            station = map(lambda x: x[1], station)
            nstmns = 0
            nstcmb = 0
            # Combine data.  See to.SBBXgrid.f lines 301 to 376
            # Note: stations.monm is equivalent to MONM0 in the Fortran code.
            # It is the maximum length of any station record.
            wt = [0.0]*stations.monm
            avg = [XBAD]*stations.monm
            if len(station) == 0:
                swrite(subboxf, dmin=XBAD)
                log.write('\nNo stations for center %+05.1f%+06.1f\n' %
                  (centre))
                continue
            # Initialise data with first station
            # See to.SBBXgrid.f lines 315 to 330
            nstmns = station[0].ngood
            nstcmb = 1
            id0 = [station[0].id]
            if audit :
                print "%+06.2f%+06.2f%+07.2f%+07.2f %s %r" % (subbox[0], subbox[1],
                    subbox[2], subbox[3], station[0].id, wti[0])
            # :todo: increment use count here. NUSEID
            wmax = wti[0]
            wtm = km * [wti[0]]
            bias = km * [0.0]
            offset = station[0].databeg
            avg[offset:station[0].dataend] = station[0].data
            a = station[0].data # just a temporary
            for i in range(len(a)):
                if a[i] < BAD :
                    wt[i + offset] = wti[0]
            # Add in the remaining stations
            # See to.SBBXgrid.f lines 331 and following
            for i in range(1,len(station)):
                dnew = [XBAD]*stations.monm
                dnew[station[i].databeg:station[i].dataend] = station[i].data
                # index, 0-based, of first year with data
                nf1 = station[i].databeg // km
                # one more than the index of last year with data
                nl1 = 1 + (station[i].dataend-1) // km
                if audit :
                    print "ADD %s %r" % (station[i].id, wti[i])
                nsm = combine(avg, wt, dnew, nf1, nl1, wti[i], wtm, station[i].id)
                nstmns += nsm
                if nsm == 0:
                    continue
                nstcmb += 1
                id0.append(station[i].id)
                # :todo: increment use count here
                if wmax < wti[i]:
                    wmax = wti[i]
            # See to.SBBXgrid.f line 354
            # This is conditional in the Fortran code, but in practice
            # NFB is always bigger than 0, so TAVG always gets called.
            tavg(avg, len(avg)//km, base_year[0]-stations.yrbeg,
                # :todo: remove dummy 99s
                base_year[1]-stations.yrbeg+1, 99, 99)
            # Subtract BIAS, then scale and write the result to disk
            m = 0
            for y in range(stations.monm//km):
                for k in range(km):
                    if avg[m] < XBAD:
                        avg[m] = stations.scale*(avg[m]-bias[k])
                    m += 1
                    # :todo: increment LENC
            swrite(subboxf, dmin=radius*(1-wmax))
    flush_out(subboxf)


def step3(label='GHCN.CL.PA', radius=1200, audit=None):
    """Do STEP3 of the GISTEMP algorithm.  label specifies the common part
    of the filenames for the 6 input files.  The filenames that are
    actually used as inputs will be formed by prepending "Ts." to label
    and appending ".n" where n is an integer from 1 to 6.  radius is the
    RCRIT value used in the algorithm, measured in Kilometres.
    """

    # Label string including radius
    labelr = 'Ts.%(label)s.%(radius)d' % locals()

    subbox_grid_output = open('work/SBBX1880.%s' % labelr, 'wb')
    box_grid_output = open('work/BX.%s' % labelr, 'wb')
    station_use_output = open('work/statn.use.%s' % labelr, 'wb')

    # Open the 6 input files
    infile = map(lambda n: open('work/Ts.%s.%d' % (label,n+1), 'rb'), range(6))

    subbox_grid(infile,
        subbox_grid_output, box_grid_output, station_use_output,
        radius=radius, audit=audit)

    # if [[ $rad -eq 1200 ]] ; then ./zonav $label ; fi
    # ./trimSBBX SBBX1880.Ts.${label}.$rad


# Mostly for debugging.
# Can be expensive.  about 50 CPU seconds on drj's MacBook.
# Note: checksubbox(0.06) emits one warning, for the Amundsen--Scott
# station at the South Pole.
def checksubboxsmall(r, stations=None, grid=None):
    """Check that for the grid, which defaults to that returned by
    eqarea.grid8k(), every station in each grid box is within r of
    the box centre.  stations should be a Stations instance, and
    defaults to fst().

    r is given as an arc-length in radians.

    This check is useful because the GISTEMP code assigns a positive
    weight to stations that are in a subbox but not within the
    critical distance; we currently dispense with these stations
    altogether, only considering stations within the critical
    distance, but we expect there to be no stations outside the
    critical distance but in the subbox.  In other words, we expect
    that r is large with respect to the subbox.  This actually
    checks this is so.
    """

    # Ravenbrook
    import eqarea

    if grid is None:
        grid = eqarea.grid8k()
    if stations is None:
        stations = fst()

    for box in grid:
        # All stations in box
        l = list(stations.inbox(*box))
        # Filtered by the radius r
        c = list(incircle(l, r, *eqarea.centre(box)))
        if len(c) != len(l):
            warn('box with centre %s has stations outside distance %f' %
                (eqarea.centre(box), r))

# Mostly for debugging and development.  Avoid public use.
def fst(prefix='Ts.GHCN.CL.PA.'):
    """Load station data from file matching the glob "prefix?" and
    return a Stations instance.
    """

    # http://www.python.org/doc/2.3.5/lib/module-glob.html
    import glob

    fl = map(lambda n: open(n, 'rb'), glob.glob(prefix + '?'))
    return Stations(fl)

# Strictly utility, not part of the public interface.
def singlefloat(x):
    """Convert 32-bit (unsigned) integer to single-precision floating
    point.  Actually, a Python float (hence a double) of the same value
    as the single-precision float is returned.
    """

    return struct.unpack('f', struct.pack('I', x))[0]

def main(argv=None):
    if argv is None:
        argv = sys.argv
    audit = None
    opt,arg = getopt.getopt(argv[1:], 'a:')
    for o,v in opt:
        if o == '-a':
            audit = eval(v, {'__builtins__':None})
    return step3(audit=audit)

if __name__ == '__main__':
    main()
