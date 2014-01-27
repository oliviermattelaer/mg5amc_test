#! /usr/bin/env python

################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

"""This is a simple script to create a release for MadGraph5_aMC@NLO, based
on the latest Bazaar commit of the present version. It performs the
following actions:

1. bzr branch the present directory to a new directory
   MadGraph5_vVERSION

4. Create the automatic documentation in the apidoc directory -> Now tar.gz

5. Remove the .bzr directory

6. tar the MadGraph5_vVERSION directory.
"""

import sys

if not sys.version_info[0] == 2 or sys.version_info[1] < 6:
    sys.exit('MadGraph5_aMC@NLO works only with python 2.6 or later (but not python 3.X).\n\
               Please upgrate your version of python.')

import glob
import optparse
import logging
import logging.config
import time

import os
import os.path as path
import re
import shutil
import subprocess
import urllib

from datetime import date

# Get the parent directory (mg root) of the script real path (bin)
# and add it to the current PYTHONPATH

root_path = path.split(path.dirname(path.realpath( __file__ )))[0]
sys.path.append(root_path)

import madgraph.various.misc as misc
from madgraph import MG5DIR

# Write out nice usage message if called with -h or --help
usage = "usage: %prog [options] [FILE] "
parser = optparse.OptionParser(usage=usage)
parser.add_option("-l", "--logging", default='INFO',
                  help="logging level (DEBUG|INFO|WARNING|ERROR|CRITICAL) [%default]")
(options, args) = parser.parse_args()
if len(args) == 0:
    args = ''

# Set logging level according to the logging level given by options
logging.basicConfig(level=vars(logging)[options.logging],
                    format="%(message)s")

# 0. check that all modification are committed in this directory
#    and that the date/UpdateNote are up-to-date
diff_result = subprocess.Popen(["bzr", "diff"], stdout=subprocess.PIPE).communicate()[0] 

if diff_result:
    logging.warning("Directory is not up-to-date. The release follow the last committed version.")
    answer = raw_input('Do you want to continue anyway? (y/n)')
    if answer != 'y':
        exit()

release_date = date.fromtimestamp(time.time())
for line in file(os.path.join(MG5DIR,'VERSION')):
    if 'version' in line:
        logging.info(line)
        version = line.rsplit('=')[1].strip()
    if 'date' in line:
        if not str(release_date.year) in line or not str(release_date.month) in line or \
                                                           not str(release_date.day) in line:
            logging.warning("WARNING: The release time information is : %s" % line)
            answer = raw_input('Do you want to continue anyway? (y/n)')
            if answer != 'y':
                exit()

Update_note = file(os.path.join(MG5DIR,'UpdateNotes.txt')).read()
if version not in Update_note:
    logging.warning("WARNING: version number %s is not found in \'UpdateNotes.txt\'" % version)
    answer = raw_input('Do you want to continue anyway? (y/n)')
    if answer != 'y':
        exit()

# 1. Adding the file .revision used for future auto-update.
# Provide this only if version is not beta/tmp/...
pattern = re.compile(r'''[\d.]+$''')
if pattern.match(version):
    #valid version format
    # Get current revision number:
    p = subprocess.Popen(['bzr', 'revno'], stdout=subprocess.PIPE)
    rev_nb = p.stdout.read().strip()
    logging.info('find %s for the revision number -> starting point for the auto-update' % rev_nb)  
else:
    logging.warning("WARNING: version number %s is not in format A.B.C,\n" % version +\
         "in consequence the automatic update of the code will be deactivated" )
    answer = raw_input('Do you want to continue anyway? (y/n)')
    if answer != 'y':
        exit()
    rev_nb=None

# checking that the rev_nb is in a reasonable range compare to the old one.
if rev_nb:
    rev_nb_i = int(rev_nb)
    try:
        filetext = urllib.urlopen('http://madgraph.phys.ucl.ac.be/mg5_build_nb')
        web_version = int(filetext.read().strip())            
    except (ValueError, IOError):
        logging.warning("WARNING: impossible to detect the version number on the web")
        answer = raw_input('Do you want to continue anyway? (y/n)')
        if answer != 'y':
            exit()
        web_version = -1
    else:
        logging.info('version on the web is %s' % web_version)
    if web_version +1 == rev_nb_i or web_version == -1:
        pass # this is perfect
    elif rev_nb_i in [web_version+i for i in range(1,4)]:
        logging.warning("WARNING: current version on the web is %s" % web_version)
        logging.warning("Please check that this (small difference) is expected.")
        answer = raw_input('Do you want to continue anyway? (y/n)')
        if answer != 'y':
            exit()
    elif web_version < rev_nb_i:
        logging.warning("CRITICAL: current version on the web is %s" % web_version)
        logging.warning("This is a very large difference. Indicating a wrong manipulation.")
        logging.warning("and can creates trouble for the auto-update.")
        answer = raw_input('Do you want to continue anyway? (y/n)')
        if answer != 'y':
            exit()
    else:
        logging.warning("CRITICAL: current version on the web is %s" % web_version)
        logging.warning("This FORBIDS any auto-update for this version.")
        rev_nb=None
        answer = raw_input('Do you want to continue anyway? (y/n)')
        if answer != 'y':
            exit()                        
# 1. bzr branch the present directory to a new directory
#    MadGraph5_vVERSION

filepath = "MG5_aMC_v" + misc.get_pkg_info()['version'].replace(".", "_")
filename = "MG5_aMC_v" + misc.get_pkg_info()['version'] + ".tar.gz"
if path.exists(filepath):
    logging.info("Removing existing directory " + filepath)
    shutil.rmtree(filepath)

logging.info("Branching " + MG5DIR + " to directory " + filepath)
status = subprocess.call(['bzr', 'branch', MG5DIR, filepath])
if status:
    logging.error("bzr branch failed. Script stopped")
    exit()

# 1. Remove the .bzr directory and clean bin directory file,
#    take care of README files.

shutil.rmtree(path.join(filepath, '.bzr'))
for data in glob.glob(path.join(filepath, 'bin', '*')):
    if not data.endswith('mg5') and not data.endswith('mg5_aMC'):
        os.remove(data)
os.remove(path.join(filepath, 'README.developer'))
shutil.move(path.join(filepath, 'README.release'), path.join(filepath, 'README'))


# 1. Add information for the auto-update
if rev_nb:
    fsock = open(os.path.join(filepath,'input','.autoupdate'),'w')
    fsock.write("version_nb   %s\n" % rev_nb)
    fsock.write("last_check   %s\n" % int(time.time()))
    fsock.close()
    
# 1. Copy the .mg5_configuration_default.txt to it's default path
shutil.copy(path.join(filepath, 'input','.mg5_configuration_default.txt'), 
            path.join(filepath, 'input','mg5_configuration.txt'))
shutil.copy(path.join(filepath, 'input','proc_card_default.dat'), 
            path.join(filepath, 'proc_card.dat'))


# 1.1 Change the trapfpe.c code to an empty file
os.remove(path.join(filepath,'Template','NLO','SubProcesses','trapfpe.c'))
create_empty = open(path.join(filepath,'Template','NLO','SubProcesses','trapfpe.c'),'w')
create_empty.close()

# 2. Create the automatic documentation in the apidoc directory
try:
    status1 = subprocess.call(['epydoc', '--html', '-o', 'apidoc',
                               'madgraph', 'aloha',
                               os.path.join('models', '*.py')], cwd = filepath)
except:
    logging.error("Error while trying to run epydoc. Do you have it installed?")
    logging.error("Execution cancelled.")
    sys.exit()

if status1:
    logging.error('Non-0 exit code %d from epydoc. Please check output.' % \
                 status)
    sys.exit()
# tarring the apidoc directory
status2 = subprocess.call(['tar', 'czf', 'doc.tgz', 'apidoc'], cwd=filepath)

if status2:
    logging.error('Non-0 exit code %d from tar. Please check result.' % \
                 status)
    sys.exit()
else:
    # remove the apidoc file.
    shutil.rmtree(os.path.join(filepath,'apidoc'))

# 3. tar the MadGraph5_vVERSION directory.

logging.info("Create the tar file " + filename)
# clean all the pyc
os.system("cd %s;find . -name '*.pyc' -delete" % filepath)
status2 = subprocess.call(['tar', 'czf', filename, filepath])

if status2:
    logging.error('Non-0 exit code %d from tar. Please check result.' % \
                 status)
    sys.exit()

logging.info("Running tests on directory %s", filepath)


logging.config.fileConfig(os.path.join(root_path,'tests','.mg5_logging.conf'))
logging.root.setLevel(eval('logging.CRITICAL'))
for name in logging.Logger.manager.loggerDict.keys():
    logging.getLogger(name).setLevel(eval('logging.CRITICAL'))
logging.getLogger('cmdprint').setLevel(eval('logging.CRITICAL'))
logging.getLogger('tutorial').setLevel(eval('logging.CRITICAL'))

# Change path to use now only the directory comming from bzr
sys.path.insert(0, os.path.realpath(filepath))
import tests.test_manager as test_manager

# reload from the bzr directory the element loaded here (otherwise it's 
#mixes the path for the tests
import madgraph
reload(madgraph)
import madgraph.various
reload(madgraph.various)
import madgraph.various.misc
reload(madgraph.various.misc)


test_results = test_manager.run(package=os.path.join('tests',
                                                     'unit_tests'))

if test_results.errors:
    logging.error("Removing %s and quitting..." % filename)
    os.remove(filename)
    exit()


a_test_results = test_manager.run(package=os.path.join('tests',
                                                       'acceptance_tests'),
                                  )
# Set logging level according to the logging level given by options
logging.basicConfig(level=vars(logging)[options.logging],
                    format="%(message)s")
logging.root.setLevel(vars(logging)[options.logging])

if a_test_results.errors:
    logging.error("Removing %s and quitting..." % filename)
    os.remove(filename)
    exit()

p_test_results = test_manager.run(['test_short_.*'],
                                  re_opt=0,
                                  package=os.path.join('tests',
                                                       'parallel_tests')
                                  )

if not test_results.wasSuccessful():
    logging.error("Failed %d unit tests, please check!" % \
                    (len(test_results.errors) + len(test_results.failures)))

if not a_test_results.wasSuccessful():
    logging.error("Failed %d acceptance tests, please check!" % \
                  (len(a_test_results.errors) + len(a_test_results.failures)))

if not p_test_results.wasSuccessful():
    logging.error("Failed %d parallel tests, please check!" % \
                  (len(p_test_results.errors) + len(p_test_results.failures)))

if p_test_results.errors:
    logging.error("Removing %s and quitting..." % filename)
    os.remove(filename)
    exit()


try:
    os.remove("%s.asc" % filename)
except:
    pass

try:
    status1 = subprocess.call(['gpg', '--armor', '--sign', '--detach-sig',
                               filename])
    if status1 == 0:
        logging.info("gpg signature file " + filename + ".asc created")
except:
    logging.warning("Call to gpg to create signature file failed. " +\
                    "Please install and run\n" + \
                    "gpg --armor --sign --detach-sig " + filename)


if not a_test_results.failures and not test_results.failures and not p_test_results.failures:
    logging.info("All good. Removing temporary %s directory." % filepath)
    shutil.rmtree(filepath)
else:
    logging.error("Some failures - please check before using release file")
    if p_test_results.failures:
        logging.error('This include discrepancy in parallel test please be carefull')

logging.info("Thanks for creating a release.")
