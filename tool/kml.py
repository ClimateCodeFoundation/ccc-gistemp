#!/usr/bin/env python
# $URL$
# $Rev$

"""
Convert v2.inv file to a KML file for browsing in Google Earth.

Usage:
python tool/kml.py [v2.temperature.inv] [> ghcn.kml]

If no input file is specified, it defaults to using input/v2.inv.
"""


# Images for the icons that denote the "urban flag" in the v2.mean file;
# this flag is either R (rural), S (semi-urban), U (urban).
urbanicon = dict(zip('RSU', """
http://google-maps-icons.googlecode.com/files/forest.png
http://google-maps-icons.googlecode.com/files/home.png
http://google-maps-icons.googlecode.com/files/factory.png
""".split()))

def xmlquote(x):
    import re

    def it(s):
        if s == '&':
            return '&amp;'
        if s == '<':
            return '&lt;'

    return re.sub('[&<]', it, x)

def kml(inp, out):
    """Convert the input *inp*, which is assumed to be (an open file
    object on) the NOAA file v2.temperature.inv file or one in an
    equivalent format, to a KML file.  The output is written to *out*.

    The KML file can be opened with Google Earth.
    """

    import re

    # Most of the KML knowledge was scraped from
    # http://code.google.com/apis/kml/documentation/kml_tut.html
    # and
    # http://code.google.com/apis/kml/documentation/kmlreference.html
    print '<?xml version="1.0" encoding="UTF-8"?>'
    print '<kml xmlns="http://www.opengis.net/kml/2.2">'
    print '<Document>'
    for code,url in urbanicon.items():
        print """
<Style id="%(code)s">
  <IconStyle><Icon>
    <href>%(url)s</href>
  </Icon></IconStyle>
</Style>""" % locals()

    for l in inp:
        # See ftp.ncdc.noaa.gov/v2.read.inv.f
        stationid = l[0:11]
        # A 30 character name field.  It often contains a country and
        # lots of internal spaces.  We reduce all multiples of SPACE
        # character to two SPACE characters.
        name = l[12:42]
        name = name.strip()
        name = re.sub(' {2,}', '  ', name)
        name = xmlquote(name)

        # v2.inv has north then east, KML has east then north.
        north = l[43:49]
        east = l[50:57]
        north = north.strip()
        east = east.strip()
        # used to determine icon (see urbanicon, above)
        urbanflag = l[67]
        print '<Placemark>'
        print '<name>%s</name>' % name
        print '<description>%s</description>' % stationid
        print '<Point><coordinates>'
        print "%(east)s,%(north)s" % locals()
        print '</coordinates></Point>'
        print '<styleUrl>#%s</styleUrl>' % urbanflag
        print '</Placemark>'
    print '</Document>'
    print '</kml>'

def main():
    import sys

    if len(sys.argv) > 1:
        f = open(sys.argv[1])
    else:
        f = open('input/v2.inv')
    return kml(f, sys.stdout)

if __name__ == '__main__':
    main()
