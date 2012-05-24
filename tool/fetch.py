#!/usr/bin/env python
# $URL$
# $Rev$
#
# fetch.py
#
# David Jones and Nick Barnes, Climate Code Foundation.
# Copyright (C) 2008-2010 Ravenbrook Limited.
# Copyright (C) 2011 Climate Code Foundation.

"""
fetch.py [--help] [--list] [--force] [--store <dir>] [--config <file>] [pattern] ...

Script to fetch (download from the internet) the inputs required for
the ccc-gistemp program.

The groups, bundles, bundle members, and individual files fetchable
are defined by the configuration file specified by --config <file>
(default 'config/sources').

Everything fetched ends up in the directory specified by --store <dir>
(default 'input').

If no such arguments are given, the default is to fetch those files in
the default group.  In the config file provided with ccc-gistemp, that
means the files required for normal ccc-gistemp operation.

Any arguments are treated as regular expressions and matched against
groups, bundles, files, or bundle members (in that order).  The first
matching item for each argument is fetched.  

Unless --force is set, no file that already exists is created.

--list lists all things that can be fetched.

The config file syntax is as follows:

    Comments begin '#' and run to the end of the line.
    Every line begins with a keyword, followed by a colon.
    Keywords are not case-sensitive.

    A fetchable item is either:

        file: <url> [<local filename>]

    which denotes a fetchable item which is also a source dataset, or

        bundle: <url> [<local filename>]

    which denotes a 'bundle': a file which can be unpacked into
    a number of files, one or more of which may be source datasets,
    identified thusly:

        member: <pattern> [<local filename>]

    <pattern> is a regular expression matching the tail of a pathname
    within the most-recently described bundle.

    The system works out for itself how to unpack a bundle.  These may
    be based on the bundle's name or contents: you shouldn't have to
    worry about it.

    <local filename> in each of the above is an optional name to give
    the fetched item or extracted member.  If absent, the system uses
    a filename derived from the fetched item or extracted member.

    <url> may be any ftp:// or http:// URL.  It may also be of this
    form:

       ftpmatch://<site>/<path>/<pattern>

    In which case the directory <path> on the FTP site <site> is
    searched for filenames matching <pattern> and the last such file
    is fetched.  This 'feature' was developed for USHCN version 2
    datasets.

    All the contents of this file may be divided into disjoint
    'groups'.  Each group is named.  Groups are introduced with group
    lines:

      group: <group name>

    The default group name is the empty string.
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://docs.python.org/release/2.4.4/lib/module-os.html
import os
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys
# http://www.python.org/doc/2.4.4/lib/module-urllib.html
import urllib

import itertools
import re

# http://www.python.org/doc/2.4.4/lib/module-tarfile.html
# Conditionally import our modified tarfile for Python 2.4.x, 2.5, and
# 2.5.1
if sys.version_info[:3] <= (2, 5, 1):
    import ccc_tarfile as tarfile
else:
    import tarfile

# Same for zipfile.  We need the open method of a ZipFile object, but
# that's only on Python 2.6 and above.  ccc_zipfile is a copy of the
# Python zipfile module from 2.6 and happily it works on Python 2.4.
if sys.version_info[:2] <= (2, 5):
    import ccc_zipfile as zipfile
else:
    import zipfile

class Fetcher(object):
    def __init__(self, **kwargs):
        self.force = kwargs.pop('force', False)
        self.output = kwargs.pop('output', sys.stdout)
        self.prefix = kwargs.pop('prefix', 'input/')
        self.config_file = kwargs.pop('config_file', 'config/sources')
        self.requests = kwargs.pop('requests', None)

    def fetch(self):
        (bundles, files) = self.find_requests(self.requests)
        for url, local in files:
            self.fetch_one(url, local)
        for ((url, local), members) in bundles.items():
            self.fetch_one(url, local, members=members)

    def make_prefix(self):
        try:
            os.makedirs(self.prefix)
        except OSError:
            # Expected if the directories already exist.
            pass

    def key_lines(self):
        comment_re = re.compile(r'((.*?[^\\])??)#')
        key_re = re.compile('^([a-zA-Z_]+)\s*:\s*(.*)$')
        for (no, l) in itertools.izip(itertools.count(1), open(self.config_file)):
            m = comment_re.match(l)
            if m:
                bare = m.group(1)
            else:
                bare = l
            bare = bare.strip()
            # ignore blank lines
            if len(bare) == 0:
                continue
            m = key_re.match(bare)
            if m:
                yield (no, m.groups())
            else:
                raise Error("%s:%d: malformed line '%s'" % (self.config_file, no, l.strip()))

    def read_config(self):
        valid_keys=dict(group  = re.compile(r'^\s*(.*?)\s*$'),
                        file   = re.compile(r'^([^\s]+)(\s+.*)?\s*$'),
                        bundle = re.compile(r'^([^\s]+)(\s+.*)?\s*$'),
                        member = re.compile(r'^([^\s]+)(\s+.*)?\s*$'))
        group=''
        config={'': dict(files = [], bundles = {})}
        for (no, (k,v)) in self.key_lines():
            k = k.lower()
            if k not in valid_keys:
                raise Error("%s:%d: unknown key '%s'" % (filename, no, k))
            m = valid_keys[k].match(v)
            if not m:
                raise Error("%s:%d: malformed '%s' line" %(filename, no, k))

            # 'bundle' only persists over 'member' lines.
            if k != 'member':
                bundle = None

            if k == 'group':
                group = m.group(1)
                config[group] = dict(files=[], bundles={})
            elif k == 'file':
                config[group]['files'].append(m.groups())
                pattern = m.group(1)
                local = m.group(2)
            elif k == 'bundle':
                bundle = m.groups()
                members = []
                config[group]['bundles'][bundle] = members
                pattern = m.group(1)
                local = m.group(2)
            elif k == 'member':
                if bundle is None:
                    raise Error("%s:%d: 'member' line with no bundle." %(filename, no))
                config[group]['bundles'][bundle].append(m.groups())
                pattern = m.group(1)
                local = m.group(2)
        return config

    def list_things(self):
        """List the things that we know how to fetch."""

        config = self.read_config()
        group_names = config.keys()
        group_names.sort()
        for g in group_names:
            if g == '':
                self.output.write("Default group: \n")
            else:
                self.output.write("Group '%s':\n" % g)
            bs = config[g]['bundles'].items()
            bs.sort()
            for ((pattern, local), members) in bs:
                self.output.write("  bundle '%s':\n" % pattern)
                if local:
                    self.output.write("   (read to '%s')\n" % local)
                for (pattern, local) in members:
                    self.output.write("    member '%s'\n" % pattern)
                    if local:
                        self.output.write("    (read to '%s')\n" % local)
            fs = config[g]['files']
            fs.sort()
            for (pattern, local) in fs:
                self.output.write("  file '%s'\n" % pattern)
                if local:
                    self.output.write("   (read to '%s')\n" % local)

    def find_requests(self, requests):
        config = self.read_config()
        bundles = {}
        files = []
        def add(fs, bs):
            for f in fs:
                files.append(f)
            for (b,ms) in bs.items():
                bundles[b] = bundles.get(b,[]) + ms
        if not requests:
            requests=['']
        for request in requests:
            if request in config:
                add(config[request]['files'], config[request]['bundles'])
                requests.remove(request)
        for request in requests:
            for group_name in config.keys():
                if re.search(request, group_name):
                    self.output.write("No group named '%s', using '%s' instead.\n"
                                      % (request, group_name))
                    add(config[group_name]['files'], config[group_name]['bundles'])
                    requests.remove(request)
        for request in requests:
            for dict in config.values():
                for (b,ms) in dict['bundles'].items():
                    (pattern, local) = b
                    if re.search(request, pattern) or (local is not None and re.search(request, local)):
                        self.output.write("No group matching '%s',\n"
                                          "    using bundle '%s:%s' instead.\n"
                                          % (request, pattern, local))
                        add([], {(pattern, local): ms})
                        requests.remove(request)
        for request in requests:
            for dict in config.values():
                for (pattern, local) in dict['files']:
                    if re.search(request, pattern) or (local is not None and re.search(request, local)):
                        self.output.write("No group or bundle matching '%s',\n"
                                          "    using file '%s:%s' instead.\n"
                                          % (request, pattern, local))
                        add([(pattern, local)], {})
                        requests.remove(request)
        for request in requests:
            for dict in config.values():
                for (b,ms) in dict['bundles'].items():
                    for (pattern, local) in ms:
                        if re.search(request, pattern) or (local is not None and re.search(request, local)):
                            self.output.write("No group or bundle matching '%s',\n"
                                              "    using member '%s:%s'\n"
                                              "    of bundle '%s:%s' instead.\n"
                                              % (request, pattern, local, b[0], b[1]))
                            add([], {b: [(pattern, local)]})
                            requests.remove(request)
        if requests:
            raise Error("Don't know how to fetch these items: %s" % requests)
        return (bundles, files)

    def fetch_one(self, url, local, members=[]):
        m = re.match('([a-z]+)://([^/]+)/(.*/)([^/]+)$', url)
        if m is None:
            raise Error("Malformed URL '%s'" % url)
        protocol = m.group(1)
        if protocol in 'http ftp'.split():
            self.fetch_url(url, local, members)
        elif protocol == 'ftpmatch':
            host = m.group(2)
            path = m.group(3)
            pattern = m.group(4)
            self.ftpmatch(host, path, pattern, local, members)
        else:
            raise Error("Unknown protocol '%s' in URL '%s'" % (protocol, url))

    def fetch_url(self, url, local, members):
        import os
        if local is None:
            local=url.split('/')[-1]
        name = os.path.join(self.prefix, local.strip())
        if os.path.exists(name) and not self.force:
            self.output.write("%s already exists.\n" % name)
        else:
            self.make_prefix()
            self.output.write("Fetching %s to %s\n" % (url, name))
            urllib.urlretrieve(url, name, progress_hook(self.output))
            self.output.write('\n')
            self.output.flush()
            if not os.path.exists(name):
                raise Error("Fetching %s to %s failed." % (url, name))
        if os.path.getsize(name) == 0:
            raise Error("%s is empty." % name)
        if members:
            self.extract(name, members)

    def ftpmatch(self, host, path, pattern, local, members):
        regexp = re.compile(pattern)
        # http://www.python.org/doc/2.4.4/lib/module-ftplib.html
        import ftplib

        remote = ftplib.FTP(host, 'ftp', 'info@climatecode.org')
        remote.cwd(path)
        dir = remote.nlst()
        good = filter(regexp.match, dir)
        good.sort()
        if not good:
            raise Error("Could not find any file matching '%s' at ftp://%s/%s" % (pattern, host, path))
        remotename = good[-1]
        path = path.strip('/')
        self.fetch_url('ftp://%s/%s/%s' % (host, path, remotename), local, members)

    def extract(self, name, members):
        exts = name.split('.')
        if exts[-1] in 'gz bz bz2'.split():
            exts = exts[:-1]
        if exts[-1] in 'tar tgz tbz tbz2'.split():
            self.extract_tar(name, members)
        elif exts[-1] in 'zip'.split():
            self.extract_zip(name, members)

    def extract_tar(self, name, members):
        # Could figure out compression type here, and pass it in to
        # tarfile.open, but apparently these days the tarfile module
        # does it for us.

        # The first argument, an empty string, is a dummy which works around
        # a bug in Python 2.5.1.  See
        # http://code.google.com/p/ccc-gistemp/issues/detail?id=26
        tar = tarfile.open('', mode='r', fileobj=open(name,'r'))
        for info in tar:
            # would like to use 'any', but that requires Python 2.5
            matches = [member for member in members if re.search(member[0]+'$', info.name)]
            if matches:
                if len(matches) > 1:
                    raise Error("Multiple patterns match '%s': %s" % (info.name, matches))
                members.remove(matches[0])
                local = matches[0][1]
                if local is None:
                    local = info.name.split('/')[-1]
                local = os.path.join(self.prefix, local.strip())
                if os.path.exists(local) and not self.force:
                    self.output.write("  ... %s already exists.\n" % local)
                else:
                    self.make_prefix()
                    out = open(local, 'wb')
                    self.output.write("  ... %s from %s.\n" % (local, info.name))
                    # The following used to be simply
                    # ``out.writelines(tar.extractfile(info))``, but the Python2.4
                    # tarfile.py does not provide iteration support.
                    member = tar.extractfile(info)
                    while True:
                        buf = member.read(4096)
                        if not buf:
                            break
                        out.write(buf)
        if members:
            raise Error("Couldn't find these members in '%s': %s" % (name, [member[0] for member in members]))

    def extract_zip(self, name, members):
        z = zipfile.ZipFile(name)
        for entry in z.namelist():
            matches = [member for member in members if re.search(member[0]+'$', entry)]
            if matches:
                if len(matches) > 1:
                    raise Error("Multiple patterns match '%s': %s" % (entry, matches))
                members.remove(matches[0])
                local = matches[0][1]
                if local is None:
                    local = entry.split('/')[-1]
                local = os.path.join(self.prefix, local.strip())
                if os.path.exists(local) and not self.force:
                    self.output.write("  ... %s already exists.\n" % local)
                else:
                    self.make_prefix()
                    # Only works for text files.
                    out = open(local, 'w')
                    self.output.write("  ... %s from %s.\n" % (local, entry))
                    src = z.open(entry)
                    while True:
                        s = src.read(4096)
                        if not s:
                            break
                        out.write(s)
                    out.close()
                    src.close()
        if members:
            raise Error("Couldn't find these members in '%s': %s" % (name, [member[0] for member in members]))

def progress_hook(out):
    """Return a progress hook function, suitable for passing to
    urllib.retrieve, that writes to the file object *out*.
    """

    def it(n, bs, ts):
        got = n*bs
        if ts < 0:
            outof = ''
        else:
            # On the last block n*bs can exceed ts, so we clamp it
            # to avoid awkward questions.
            got = min(got, ts)
            outof = '/%d [%d%%]' % (ts, 100*got//ts)
        out.write("\r  %d%s" % (got, outof))
        out.flush()
    return it

class Error(Exception):
    """Some sort of problem with fetch."""

# Guido's main, http://www.artima.com/weblogs/viewpost.jsp?thread=4829
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    write_list = False
    kwargs = dict()
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "", ["help", "list", "force", "store=", "config="])
            for o,a in opts:
                if o in ('--help',):
                    print __doc__
                    return 0
                if o == '--list':
                    write_list = True
                if o == '--force': 
                    kwargs.update(force=True)
                if o == '--config':
                    kwargs.update(config_file=a)
                if o == '--store':
                    kwargs.update(prefix=a)
        except getopt.error, msg:
             raise Usage(msg)
        fetcher = Fetcher(**kwargs)
        if write_list:
            fetcher.list_things()
        else:
            fetcher.fetch()
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    return 0

if __name__ == "__main__":
    sys.exit(main())
