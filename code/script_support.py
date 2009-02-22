#!/usr/bin/env python
"""Some useful support for Python scripts.

TODO: Describle more.
"""
__docformat__ = "restructuredtext"

import sys
import os

def parseIntArg(idx, args):
    """Parse ``args[idx]`` as an integer.

    If the argument is not a valid integer then the program exits with
    code == 1.
    """
    try:
        return int(args[idx])
    except ValueError:
        sys.exit("Argument %d must be an integer" % (idx + 1,))


def makeParser(usage):
    """Create an optparse.OptionParser with common options.

    This is a simple factory that creates an ``optparse.OptionParser`` instance
    and adds some common options.
    
    Nothing fancy, but it does provide for a level of consistency.

    The standard options are:

    -h, --help
        Used to display a short usage message. This is standard in the
        ``optparse.OptionParser``.
    --man
        Used to display detailed help, including the usage message.
    --verbose=level, -v
        Allow the ``verbose`` to be set.
    """
    import optparse
    parser = optparse.OptionParser(usage)

    parser.add_option("--man", action="store_true",
        help="Print a simple man page.")
    parser.add_option("-v", action="count", dest="verbose",
        help="Increase verbose level")
    parser.add_option("--no-psyco", action="store_true",
        help="Do no use psyco")
    parser.add_option("--verbose", action="store", type="int",
        metavar="<N>", default=0,
        help="Set verbose level to <N>")

    return parser


def parseArgs(parser, doc, argRange=(0, 0)):
    """Wrapper for parsing command line options.
 
    This invokes ``parser.parse_args()`` and then handles some of the common options.

    :Param parser:
        The ``optparse.OptionParser`` instance, which should typically have
        been created using `makeParser`.
    :Param doc:
        The documentation for the ``--man`` option. Typically this is the calling
        programs doc-string.
    :Param argRange:
        A tuple of ``(min, max)`` defining the minimum allowed and maximum allowed
        number of command line arguments. This range is inclusive, so (for example)
        if you require exeactly 2 arguments, use (2, 2).

    :Return:
        The tuple ``(options, args)``, as returned by ``parser.parse_args()``.
        This function may not return, for example if the ``--help`` option is
        used on the command line.
    """
    options, args = parser.parse_args()
    a, b = argRange
    if not a <= len(args) <= b:
        parser.error("Wrong number of arguments")
    if options.man:
        print doc
        print "-" * 50
        parser.print_help()
        sys.exit(0)

    return options, args


def countReport(interval, fmt="%d ", f=None):
    """Generator:

    A count generator with he side effect of writing the count every ``interal``
    iterations:

    :Param interval:
        The reporting interval.
    :Param f:
        The file to report the count to. If not set or ``None`` then sys.stderr
        is used.
    :Param fmt:
        Format string for the report. The report will be ``fmt % count``. No
        newline is automatically added.
    """
    count = 0
    while True:
        count += 1
        if count % interval == 0:
            f = f or sys.stderr
            f.write(fmt % count)
            try:
                f.flush()
            except NameError:
                pass
        yield count


def enablePysco(progFile, *functions):
    def psycoFilt(co):
        if not os.path.basename(co.co_filename) == os.path.basename(progFile):
            return False
        # print >>sys.stderr, "Psyco", co.co_name
        return True

    try:
        import psyco
        psyco.setfilter(psycoFilt)
        for f in functions:
            psyco.bind(f)
    except ImportError:
        pass


