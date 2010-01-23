#!/usr/bin/env python
"""A script to run CCC code on 'official' test data and compare with GISS
results.

This script temporarily replaces the ``input`` directory with a copy of the
``ccc-gistemp-test-2009-12-28/input`` directory. Then it executes all steps of
the CCC code, using the `run.py` module and finally used the `compare_results`
module to generate the comparison web page ``index.html``.

When this script exits, it restores the ``input`` directory to its former
state.

"""
__docformat__ = "restructuredtext"


import os
import sys
import shutil

# Make sure we can import from other CCC packages.
import extend_path

import ref_test_data
import run
import compare_results
from code import script_support


# Dummy options for when used as a module.
options = script_support.options


#: The path of a temporary directory used to save backups of files
#: that get modified while the script is running.
safe_dir = "__SAFE__"

_safe_exists_msg  = """\
The %(safe)s directory already exists. This may be because a previous run of
this program failed badly.

Check your ``input`` directory restoring
it from %(safe)s if necessary. When your input directory is restored, remove
the %(safe)s directory.

This program is designed to create %(safe)s to save files while it runs.
Normally the %(safe)s directory is deleted when the program exits.
"""


_readme_txt = """This was created byt the script tool/regression.py.

Normally this get automatically removed when the script finishes.
"""


#: A list of functions that need to be called in order to undo
#: changes made during the test run.
undo_actions = []


#{ General helper functions
def mkdir(path):
    """Make a directory, if it does not already exist.

    All missing intermediate directories are also created.

    :Param path:
        The name of the directory to create.
    """
    if not os.path.isdir(path):
        os.makedirs(path)


#{ Clean up functions
def rm_safe():
    """Cleanup function: Remove the `safe_dir` directory.

    """
    try:
        os.unlink(os.path.join(safe_dir, "README.txt"))
    except OSError:
        pass # Can happen if interrupted at just the tight moment.

    sys.stdout.write("Removing %s directory...\n" % safe_dir)
    try:
        os.rmdir(safe_dir)
    except OSError, exc:
        sys.stderr.write("BAD CLEANUP: Could not remove %s\n" % safe_dir)
        sys.stderr.write("             %s\n" % exc)


def restore_inputs():
    """Cleanup function: Restore the input directory to its old state.

    """
    sys.stdout.write("Restoring input directory...\n")
    source = os.path.join(safe_dir, "input")
    try:
        shutil.move(source, "input")
    except (OSError, shutil.Error), exc:
        sys.sdterr.write("BAD CLEANUP: Failed to move %s to %s\n" % (
            source, "input", dest))
        sys.sdterr.write("             %s\n" % exc)


def remove_test_input():
    """Cleanup function: Remove test version of input directory.

    """
    sys.stdout.write("Removing test input directory...\n")
    try:
        shutil.rmtree("input")
    except (OSError, shutil.Error), exc:
        sys.sdterr.write("BAD CLEANUP: Failed to remove input\n")
        sys.sdterr.write("             %s\n" % exc)


def clean_up():
    """Run just before program exit to undo all changes.

    Uses the `undo_actions` list to undo all the changes to the working
    directory.

    """
    for action in reversed(undo_actions):
        try:
            action()
        except Exception, exc:
            sys.stderr.write("BAD CLEANUP: Call to %s failed\n"
                    % action.func_name)
            sys.stderr.write("             %s\n" % exc)


#{ Test execution functions
def install_inputs():
    """Install the 'official' input files in the ``input`` directory.

    The original ``input`` directory is saved in the `safe_dir`.

    """
    dest = os.path.join(safe_dir, "input")
    sys.stdout.write("Moving directory %r to %r...\n" % ("input", dest))
    try:
        shutil.move("input", dest)
    except (OSError, shutil.Error), exc:
        sys.sdterr.write("Failed to move %r to %r\n" % ("input", dest))
        sys.sdterr.write("    %s\n" % exc)
        return 1
    undo_actions.append(restore_inputs)

    source = os.path.join(ref_test_data.test_data_dir, "input")
    sys.stdout.write("Copying directory %r to %r...\n" % (source, "input"))
    try:
        shutil.copytree(source, "input")
    except (OSError, shutil.Error), exc:
        sys.sdterr.write("Failed to move %r to %r\n" % (source, "input"))
        sys.sdterr.write("    %s\n" % exc)
        return 1
    undo_actions.append(remove_test_input)

    return 0


def run_test():
    """Run the regression test on the 'official' test data and check results.

    This is invoked once the `safe_dir` has been set up. This will install the
    'official' test input data, downloading and unpacking the files if
    necessary. Then the CCC code is executed, followed by the compare_results
    code.

    """
    # First off, we need the 'official' test data downloaded and unpacked
    # if necessary.
    if ref_test_data.install_and_check_test_files() != 0:
        return 1
    if install_inputs() != 0:
        return 1

    if not options.skip_run:
        #ret = run.main([sys.argv, "--steps", "3"])
        sys.stdout.write("Executing CCC code...\n")
        ret = run.main([sys.argv])
        if ret != 0:
            return 1

    test_result_dir = os.path.join(ref_test_data.test_data_dir, "result")
    sys.stdout.write("Comparing CCC results with 'official' results...\n")
    ret = compare_results.main([sys.argv[0], "result", test_result_dir])
    if ret != 0:
        return 1


#{ Main 
def main(args=[]):
    """Prepare for and run the regression test.

    The this module's documentation for an outline of what this does.

    :Param args:
        The command line arguments, parsed using ``optparse``.

    """
    # Do initial checks and setup. If this fails then the program exits without
    # returning.
    if os.path.exists(safe_dir):
        sys.stderr.write(_safe_exists_msg % {"safe": safe_dir})
        sys.exit(1)
    try:
        os.mkdir(safe_dir)
        undo_actions.append(rm_safe)
    except OSError, exc:
        sys.exit(str(exc))

    try:
        f = open(os.path.join(safe_dir, "README.txt"), "w")
        f.write(_readme_txt)
        f.close
        run_test()
    finally:
        clean_up()


def run_as_script():
    """Invoke this module as a script.

    This is invoked when this module is run as a script. This expects
    ``sys.argv`` to be correctly set up.

    """
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    parser.add_option("--skip-run", action="store_true",
            help="Do not do the run phase, skip to the compare.")

    global options
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)


if __name__ == "__main__":
    sys.exit(run_as_script())
