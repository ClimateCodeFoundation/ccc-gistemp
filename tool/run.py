#!/usr/bin/env python
# $URL$
# $Rev$
# 
# run.py -- run steps of the GISTEMP algorithm
#
# Gareth Rees, 2009-12-08

"""run.py [options] -- run steps of the GISTEMP algorithm.
Options:
   --help         Print this text.
   --steps=STEPS  Specify which steps to run, as a comma-separated list of
                  numbers from 0 to 5.  For example, --steps=2,3,5
                  The steps are run in the order you specify.
                  If this option is omitted, run all steps in order.
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os
# http://docs.python.org/release/2.4.4/lib/module-re.html
import re
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys

# Clear Climate Code
import extend_path
import giss_io

class Fatal(Exception):
    pass

# :todo: remove me
# Record the original standard output so we can log to it; in steps 2 and 5
# we'll be changing the value of sys.stdout before calling other modules that
# use "print" to generate their output.
logfile = sys.stdout

def log(msg):
    print >>logfile, msg

def mkdir(path):
    """mkdir(PATH): create the directory PATH, and all intermediate-level
    directories needed to contain it, unless it already exists."""
    if not os.path.isdir(path):
        log("... creating directory %s" % path)
        os.makedirs(path)

# Each of the run_stepN functions below takes a data object, its input,
# and produces a data object, its output.  Ordinarily the data objects
# are iterators, either produced from the previous step, or an iterator
# that feeds from a file.
def run_step0(data):
    from code import step0
    if data is None:
        data = giss_io.step0_input()
    result = step0.step0(data)
    return giss_io.step0_output(result)

def run_step1(data):
    from code import step1
    if data is None:
        data = giss_io.step1_input()
    result = step1.step1(data)
    return giss_io.step1_output(result)

def run_step2(data):
    from code import step2
    if data is None:
        data = giss_io.step2_input()
    result = step2.step2(data)
    return giss_io.step2_output(result)

def run_step3(data):
    from code import step3
    if data is None:
        data = giss_io.step3_input()
    result = step3.step3(data)
    return giss_io.step3_output(result)

def run_step3c(data):
    """An alternative to Step 3 that reads (copies) the output file
    created by the ordinary Step 3.  Effectively using the data produced
    by Step 3 without re-running it."""
    if data:
        raise Fatal("Expect to run 3c first in pipeline.")
    return giss_io.step3c_input()

def run_step4(data):
    from code import step4
    data = giss_io.step4_input(data) 
    result = step4.step4(data)
    return giss_io.step4_output(result)

def run_step5(data):
    from code import step5
    if data is None:
        data = giss_io.step5_input()
    result = step5.step5(data)
    giss_io.step5_output(result)
    return vischeck(result)

def vischeck(data):
    # Suck data through pipeline.
    for _ in data:
        pass
    log("... running vischeck")
    import vischeck
    vischeck.chartit(
      [open(os.path.join('result', 'GLB.Ts.ho2.GHCN.CL.PA.txt'))],
      out = open(os.path.join('result', 'google-chart.url'), 'w'))

    log("See result/google-chart.url")
    yield "vischeck completed"

def parse_steps(steps):
    """Parse the -s, steps, option.  Produces a list of strings."""
    steps = steps.strip()
    if not steps:
        return range(0, 6)
    result = set()
    for part in steps.split(','):
        # Part can be integer number with an optional letter suffix...
        if re.match(r'^\d+[a-z]?$', part):
            result.add(part)
        else:
            # Or a range in the form '1-3'.
            try:
                l,r = part.split('-')
                result.update(str(s) for s in range(int(l), int(r)+1))
            except ValueError:
                # Expect to catch both
                # "ValueError: too many values to unpack" when the split
                # produces too many values ("1-3-"), and
                # "ValueError: invalid literal for int() with base 10: 'a'"
                # when int fails ("1,a")
                raise Fatal("Can't understand steps argument.")

    return list(sorted(result))


def parse_options(arglist):
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)

    parser.add_option("-s", "--steps", action="store", metavar="S[,S]",
            default="",
            help="Select range of steps to run")
    parser.add_option('-p', '--parameter', action='store', default='',
        help="Redefine parameter from parameters.py during run")
    parser.add_option("--no-work_files", "--suppress-work-files",
            action="store_false", default=True, dest="save_work",
            help="Do not save intermediate files in the work sub-directory")
    options, args = parser.parse_args(arglist)
    if len(args) != 0:
        parser.error("Unexpected arguments")

    options.steps = parse_steps(options.steps)
    return options, args

def update_parameters(parm):
    """Take a parameter string from the command line and update the
    parameters.py module."""

    if not parm:
        return

    from code import parameters

    parm = parm.split(';')
    for p in parm:
        try:
            key,value = p.split('=', 1)
        except ValueError:
            raise Fatal("Can't understand parameter option: %r" % p)
        if not hasattr(parameters, key):
            print "Ignoring unknown parameter %r" % key
            continue
        # Coerce value, a string, to the same type as the existing parameter
        # value.  That works nicely for strings, ints, and floats...
        x = getattr(parameters, key)
        # ... but we need a hack for bool.
        if type(x) == bool:
            try:
                value = ['false', 'true'].index(value.lower())
            except ValueError:
                raise Fatal("Boolean parameter %r must be True or False"
                  % key)
            # Now value is 0 or 1 and the default case will correctly
            # coerce it.
        value = type(x)(value)
        setattr(parameters, key, value)


def main(argv=None):
    import sys
    # http://www.python.org/doc/2.4.4/lib/module-time.html
    import time

    if argv is None:
        argv = sys.argv
    options, args = parse_options(argv[1:])

    update_parameters(options.parameter)

    step_list = options.steps
    try:
        rootdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if os.getcwd() != rootdir:
            raise Fatal("The GISTEMP procedure must be run from the root "
                        "directory of the project.\nPlease change directory "
                        "to %s and try again." % rootdir)

        # Carry out preflight checks and fetch missing files.
        import preflight
        preflight.checkit(sys.stderr)

        # Create all the temporary directories we're going to use.
        for d in ['log', 'result', 'work']:
            mkdir(d)

        step_fn = {
            '0': run_step0,
            '1': run_step1,
            '2': run_step2,
            '3': run_step3,
            '3c': run_step3c,
            '4': run_step4,
            '5': run_step5,
        }
        
        # Record start time now, and ending times for each step.
        start_time = time.time()

        cannot = [s for s in step_list if not step_fn.has_key(s)]
        if cannot:
            raise Fatal("Can't run steps %s" % str(cannot))

        # Create a message for stdout.
        if len(step_list) == 1:
            logit = "STEP %s" % step_list[0]
        else:
            assert len(step_list) >= 2
            try:
                t = [str(s) for s in range(step_list[0], step_list[-1]+1)]
            except:
                t = []
            if step_list == t:
                logit = "STEPS %s to %s" % (step_list[0], step_list[-1])
            else:
                logit = "STEPS %s" % ', '.join(step_list)
        log("====> %s  ====" % logit)
        data = None
        for step in step_list:
            data = step_fn[step](data)
        # Consume the data in whatever the last step was, in order to
        # write its output, and hence suck data through the whole
        # pipeline.
        for _ in data:
            pass

        end_time = time.time()
        log("====> Timing Summary ====")
        log("Run took %.1f seconds" % (end_time - start_time))

        return 0
    except Fatal, err:
        sys.stderr.write(str(err))
        sys.stderr.write('\n')
        return 2


if __name__ == '__main__':
    sys.exit(main())
