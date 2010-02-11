#!/usr/bin/env python
# $URL$
# $Id: run.py 159 2010-01-12 10:39:07Z drj@pobox.com $
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
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys

# Extend sys.path to support access to all the clear climate code.
import extend_path
import tool.step0
import tool.step1
import tool.step2
import tool.step3
import tool.step4
import tool.step5

import code.step5


class Fatal(Exception):
    def __init__(self, msg):
        self.msg = msg

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


def run_step(source, sink):
    for record in source:
        sink.add_record(record)
    sink.close()


def run_step0(options):
    source = tool.step0.get_step_iter()
    sink = tool.step0.get_outputs()
    run_step(source, sink)


def run_step1(options):
    source = tool.step1.get_step_iter(options.steps, options.save_work)
    sink = tool.step1.get_outputs()
    run_step(source, sink)


def run_step2(options):
    source = tool.step2.get_step_iter(options.steps, options.save_work)
    sink = tool.step2.get_outputs()
    run_step(source, sink)


def run_step3(options):
    source = tool.step3.get_step_iter(options.steps, options.save_work)
    sink = tool.step3.get_outputs()
    run_step(source, sink)


def run_step4(options):
    source = tool.step4.get_step_iter(options.steps, options.save_work)
    sink = tool.step4.get_outputs()
    run_step(source, sink)


def run_step5(options):
    inputs = tool.step5.get_inputs(options.steps, options.save_work)

    # Save standard output so we can restore it when we're done with this step.
    old_stdout = sys.stdout

    code.step5.step5(inputs)

    log("... running vischeck")
    import vischeck
    sys.stdout = open(os.path.join('result', 'google-chart.url'), 'w')
    vischeck.main(['', os.path.join('result', 'GLB.Ts.ho2.GHCN.CL.PA.txt')])
    sys.stdout.close()

    log("See result/google-chart.url")

    # Restore standard output.
    sys.stdout = old_stdout


def parse_steps(steps):
    steps = steps.strip()
    if not steps:
        return range(0, 6)
    try:
        step = int(steps)
        return [step]

    except (TypeError, ValueError):
        try:
            a, b = [int(v) for v in steps.split(",")]
            return range(a, b + 1)

        except ValueError:
            sys.exit("Do not understand steps arg %r" % steps)


def parse_options():
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)

    parser.add_option("-s", "--steps", action="store", metavar="S[,S]",
            default="",
            help="Select range of steps to run")
    parser.add_option("--no-work_files", "--suppress-work-files",
            action="store_false", default=True, dest="save_work",
            help="Do not save intermediate files in the work sub-directory")
    options, args = parser.parse_args()
    if len(args) != 0:
        parser.error("Unexpected arguments")

    options.steps = parse_steps(options.steps)
    return options, args


def main(options, args):
    # http://www.python.org/doc/2.4.4/lib/module-time.html
    import time

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
            0: run_step0,
            1: run_step1,
            2: run_step2,
            3: run_step3,
            4: run_step4,
            5: run_step5,
        }
        
        # Record start time now, and ending times for each step.
        start_time = time.time()
        s = max(step_list)

        if not step_fn.has_key(s):
            raise Fatal("Can't run step %d" % s)

        log("====> STEPS %d to %d  ====" % (step_list[0], step_list[-1]))
        # As a special case, if we are running a sequence of steps ending with step4,
        # we need to run step3 and step4 separately. This is beacuse the steps
        # join up like this::
        #
        #         0 ---> 1 ---> 2 ---> 3 -.
        #                                  |--> 5
        #                              4 -'
        #
        if s == 4 and 3 in step_list:
            step_fn[3](options)
            step_fn[4](options)
        else:
            step_fn[s](options)

        end_time = time.time()
        log("====> Timing Summary ====")
        log("Run took %.1f seconds" % (end_time - start_time))

        return 0
    except Fatal, err:
        sys.stderr.write(err.msg)
        sys.stderr.write('\n')
        return 2


if __name__ == '__main__':
    options, args = parse_options()
    sys.exit(main(options, args))
