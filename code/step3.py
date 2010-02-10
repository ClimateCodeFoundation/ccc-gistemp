#!/usr/bin/env python
# $URL$
# $Rev$
#
# step3.py
#
# David Jones, Ravenbrook Limited, 2008-08-06

""" 
Python code reproducing the STEP3 part of the GISTEMP algorithm.

python code/step3.py [-a N] [-s X,Y,Z]

-a N
  Restrict to subbox N and producing audit information

-s X,Y,Z
  Restrict to the list of subboxes

Work in progress.

The code is derived from the Fortran GISTEMP code.

Python notes:

Much of the binary IO uses the struct module.  A quick review will
help.

I feel obliged to tell you that the Python expression N * list gives a
new list that conists of list repeated N times.  This is used every
now and then.

We also make use of list slice assignment: a[2:4] = range(2)
"""

# Clear Climate Code
import earth # required for radius.
# Clear Climate Code
import eqarea
# Clear Climate Code
import giss_data

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

    def __init__(self, infile, km, meta):
        """infile should be a sequence of file objects.  Each file
        object should be opened in binary onto a trimmed NCAR file,
        identical in format to the input files that GISTEMP's
        to.SBBXgrid.f program expects.
        
        Historically (and in emulation of the GISS GISTEMP code) there
        were be 6 zonal files, but now only one is usually supplied.
        The "for f infile" loop is a relic, but still kept.
        """
        (self.reader, ) = infile
        self.meta = meta
        self.station_records = None

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

        self.km = km

    def __getattr__(self, name):
        return getattr(self.meta, name)

    def all(self):
        """An iterator that yields every station record."""
        if self.station_records is not None:
            for record in self.station_records:
                yield record
            return

        station_records = []
        for record in self.reader:
            station_records.append(record)
            yield record
        self.station_records = station_records

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

        for record in self.all():
            st = record.station
            lat = st.lat_fixed_1
            lon = st.lon_fixed_1

            # See to.SBBXgrid.f lines 213 and following
            if lon > longe:
                lon -= 360
            elif lon < longw:
                lon += 360

            if (lats < lat < latn) and (longw < lon < longe):
                yield record
            if lats == lat and longw <= lon < longe:
                yield record
            if longw == lon and lats <= lat < latn:
                yield record


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

    for record in iterable:
        st = record.station
        s_lat, s_lon = st.lat_fixed_1, st.lon_fixed_1
        # A possible improvement in speed (which the corresponding
        # Fortran code does) would be to store the trig values of
        # the station location in the station object.
        sinlats = math.sin(s_lat*math.pi/180)
        coslats = math.cos(s_lat*math.pi/180)
        sinlons = math.sin(s_lon*math.pi/180)
        coslons = math.cos(s_lon*math.pi/180)

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
            yield (1-d, record)


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


# Equivalent of the Fortran subroutine CMBINE
# Parameters as per Fortran except that km is bound lexically, and
# nsm is returned directly.
def combine(km, XBAD, bias, avg, wt, dnew, nf1, nl1, wt1, wtm, id,
        NOVRLP=20):
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

    from itertools import izip

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
        for a,n in izip(avg[nf1*km+k: nl1*km: km],
                       dnew[nf1*km+k: nl1*km: km]):
            if a >= XBAD or n >= XBAD:
                continue
            ncom += 1
            sum += a
            sumn += n
        if ncom < NOVRLP:
            continue
        biask = float(sum-sumn)/ncom
        # Find mean bias
        wtmnew = wtm[k]+wt1
        # Note: bias is lexically captured
        bias[k] = float(wtm[k]*bias[k]+wt1*biask)/wtmnew
        wtm[k]=wtmnew
        # Update period of valid data, averages and weights
        for kn in range(nf1*km+k, nl1*km, km):
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


# TODO: This is an almost duplicate of the function in step5.
#       Make the code common.

# Equivalent to Fortran subroutine TAVG.
# See to.SBBXgrid.f lines 563 and following.
def tavg(km, XBAD, bias, data, nyrs, base, limit, nr, nc, deflt=0.0):
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


def iter_subbox_grid(stations, monm, radius=1200,
        base_year=(1951,1980), audit=None,
        just_subbox=()):
    """(This is the equivalent of to.SBBXgrid.f) Convert the input file,
    infile, into gridded datasets which are output on the file (-like)
    object box_out.

    radius specifies a radius in kilometres, it is
    equivalent to the RCRIT value in the Fortran code.

    year_begin specifies the earliest year to include in the output.

    base_year is a pair specifying the base (reference) period in
    (first_year, last_year).  The last year is included in the base
    period.
    """

    log = sys.stdout

    # trimmed input files will be assumed.

    # Parameters.  These are a mixture of computed constants, and
    # parameters that influence the behaviour of the algorithm.  Many of
    # the parameters are tempting to change, but some of the code might
    # accidentally rely on particular values.

    assert radius > 0

    radius = float(radius)

    # Critical radius as an angle of arc
    arc = radius / earth.radius

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

    BAD = giss_data.MISSING
    XBAD = giss_data.XMISSING
    km = stations.km

    regions = list(eqarea.gridsub())
    if audit:
        just_subbox = [audit]
    if just_subbox:
        # Restrict processing to region containing single gridbox of
        # interest.
        n = list(set(map(lambda x: x//100, just_subbox)))
        if len(n) > 1:
            raise Exception('Sorry, only works within one region currently.')
        regions = map(lambda x: regions[x], n)

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
        if audit:
            print 'UNSORTED LIST'
            for station in regionstations:
                print station.id
        # Descending sort by number of good records
        # Sadly we cannot use Python's sort method here, we must use an
        # emulation of the GISS Fortran.
        sort(regionstations, lambda x,y: y.good_count - x.good_count)
        if audit:
            print 'SORTED LIST'
            for station in regionstations:
                print station.id

        subboxes = list(region[1])
        if just_subbox:
            # Select subboxes for processing.
            l = map(lambda x: x % 100, just_subbox)
            subboxes = map(lambda x: subboxes[x], l)

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
            # Of possible stations in this region, filter for those within
            # radius of subbox centre.  Note that it is important that
            # the ordering within the regionstations list is retained
            # (in order to match the GISS code).
            # *station* is a list of (wt, station) pairs:
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
                box_obj = giss_data.SubboxRecord(
                    lat_S=latlon[0], lat_N=latlon[1], lon_W=latlon[2],
                    lon_E=latlon[3], stations=nstcmb, station_months=nstmns,
                    d=XBAD, series=avg)
                log.write('*')
                log.flush()
                yield box_obj
                continue
            # Initialise data with first station
            # See to.SBBXgrid.f lines 315 to 330
            nstmns = station[0].good_count
            nstcmb = 1
            id0 = [station[0].station.uid]
            if audit:
                print "%+06.2f%+06.2f%+07.2f%+07.2f %s %r" % (subbox[0], subbox[1],
                    subbox[2], subbox[3], station[0].station.uid, wti[0])
            # :todo: increment use count here. NUSEID
            wmax = wti[0]
            wtm = km * [wti[0]]
            bias = km * [0.0]
            #offset = station[0].databeg PAO
            r = station[0]
            offset = r.rel_first_month - 1
            avg[offset:offset + len(r.series_as_tenths)] = r.series_as_tenths
            a = r.series_as_tenths # just a temporary
            for i in range(len(a)):
                if a[i] < BAD :
                    wt[i + offset] = wti[0]
            # Add in the remaining stations
            # See to.SBBXgrid.f lines 331 and following
            for i in range(1,len(station)):
                # TODO: A StationMethod method to produce a padded data series
                #       would be good here. Hence we could just do:
                #           dnew = station[i].padded_series(stations.monm)
                dnew = [XBAD]*stations.monm
                aa, bb = station[i].rel_first_month, station[i].rel_last_month
                dnew[aa - 1:bb] = station[i].series_as_tenths
                # index, 0-based, of first year with data
                #nf1 = station[i].databeg // km
                # TODO: I think first_year/last_year (+1) should work here.
                nf1 = (station[i].rel_first_month - 1) // km
                # one more than the index of last year with data
                #nl1 = 1 + (station[i].dataend-1) // km
                nl1 = 1 + station[i].rel_last_month // km
                if audit:
                    print "ADD %s %r" % (station[i].id, wti[i])
                nsm = combine(km, XBAD, bias, avg, wt, dnew, nf1, nl1,
                        wti[i], wtm, station[i].uid)
                nstmns += nsm
                if nsm == 0:
                    continue
                nstcmb += 1
                # TODO: What is id0 used for?
                id0.append(station[i].uid)
                # :todo: increment use count here
                if wmax < wti[i]:
                    wmax = wti[i]
            # See to.SBBXgrid.f line 354
            # This is conditional in the Fortran code, but in practice
            # NFB is always bigger than 0, so TAVG always gets called.
            tavg(km, XBAD, bias, avg, len(avg)//km, base_year[0]-stations.yrbeg,
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
            box_obj = giss_data.SubboxRecord(n=monm,
                    lat_S=latlon[0], lat_N=latlon[1], lon_W=latlon[2],
                    lon_E=latlon[3], stations=nstcmb, station_months=nstmns,
                    d=radius*(1-wmax), series=avg[:monm])
            yield box_obj
    print >>log


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

    # Clear Climate Code
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
def fst(file='Ts.GHCN.CL.PA'):
    """Load station data from file and
    return a Stations instance.
    """

    return Stations([file])


class Step3Iterator(object):
    def __init__(self, record_source, radius=1200, year_begin=1880,
            audit=None, subbox=()):
        # Ensure record_source is a true generator, so we can pop off the
        # leading meta data before handing the station records to the next
        # phase.
        record_source = iter(record_source)
        m = record_source.next()
        self.meta = giss_data.SubboxMetaData(m.mo1, m.kq, m.mavg, m.monm,
                m.monm + 7, m.yrbeg, m.missing_flag, m.precipitation_flag,
                m.title)

        # Only monthly data supported, indicated by self.meta.mavg == 6, hence
        # km = 12. Load all the station records into a `Stations` instance.
        # TODO: The names `Stations` is slightly misleading.
        assert self.meta.mavg == 6
        km = 12
        stations = Stations([record_source], km, self.meta)

        # Output series cannot begin earlier than input series.
        year_begin = max(year_begin, self.meta.yrbeg)
        monm = (self.meta.yrbeg + (self.meta.monm // km) - year_begin) * km

        units = '(C)'
        if stations.kq == 2:
            units = '(mm)'
        title = "%20.20s ANOM %-4s CR %4dKM %s-present" % (self.meta.title,
                units, radius, year_begin)
        self.meta.mo1 = 1
        self.meta.title = title.ljust(80)

        self.box_source = iter_subbox_grid(stations, monm,
                radius=radius, base_year=(1951,1980), audit=None,
                just_subbox=subbox)

    def __iter__(self):
        return self._it()

    def _it(self):
        yield self.meta
        for box in self.box_source:
            yield box


def step3(record_source, radius=1200, year_begin=1880, audit=None, subbox=()):
    """Step 3 of the GISS processing.

    """
    return Step3Iterator(record_source,
            radius=radius, year_begin=year_begin, audit=audit, subbox=subbox)
