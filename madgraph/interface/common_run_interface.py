################################################################################
#
# Copyright (c) 2011 The MadGraph5_aMC@NLO Development team and Contributors
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
"""A user friendly command line interface to access MadGraph5_aMC@NLO features.
   Uses the cmd package for command interpretation and tab completion.
"""
from __future__ import division

import ast
import atexit
import cmath
import cmd
import glob
import logging
import math
import optparse
import os
import pydoc
import random
import re
import shutil
import signal
import stat
import subprocess
import sys
import time
import traceback
import glob
import sets

try:
    import readline
    GNU_SPLITTING = ('GNU' in readline.__doc__)
except:
    GNU_SPLITTING = True

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
root_path = os.path.split(root_path)[0]
sys.path.insert(0, os.path.join(root_path,'bin'))

# usefull shortcut
pjoin = os.path.join
# Special logger for the Cmd Interface
logger = logging.getLogger('madgraph.stdout') # -> stdout
logger_stderr = logging.getLogger('madgraph.stderr') # ->stderr


try:
    import madgraph
except ImportError:    
    # import from madevent directory
    import internal.extended_cmd as cmd
    import internal.banner as banner_mod
    import internal.shower_card as shower_card_mod
    import internal.misc as misc
    import internal.cluster as cluster
    import internal.check_param_card as check_param_card
    import internal.files as files
    import internal.histograms as histograms
    import internal.save_load_object as save_load_object
    import internal.gen_crossxhtml as gen_crossxhtml
    from internal import InvalidCmd, MadGraph5Error
    MADEVENT=True    
else:
    # import from madgraph directory
    import madgraph.interface.extended_cmd as cmd
    import madgraph.various.banner as banner_mod
    import madgraph.various.shower_card as shower_card_mod
    import madgraph.various.misc as misc
    import madgraph.iolibs.files as files
    import madgraph.various.cluster as cluster
    import madgraph.iolibs.save_load_object as save_load_object
    import madgraph.madevent.gen_crossxhtml as gen_crossxhtml
    import models.check_param_card as check_param_card
    import madgraph.various.histograms as histograms
    
    from madgraph import InvalidCmd, MadGraph5Error, MG5DIR
    MADEVENT=False

#===============================================================================
# HelpToCmd
#===============================================================================
class HelpToCmd(object):
    """ The Series of help routins in common between amcatnlo_run and
    madevent interface"""

    def help_treatcards(self):
        logger.info("syntax: treatcards [param|run] [--output_dir=] [--param_card=] [--run_card=]")
        logger.info("-- create the .inc files containing the cards information." )

    def help_set(self):
        logger.info("syntax: set %s argument" % "|".join(self._set_options))
        logger.info("-- set options")
        logger.info("   stdout_level DEBUG|INFO|WARNING|ERROR|CRITICAL")
        logger.info("     change the default level for printed information")
        logger.info("   timeout VALUE")
        logger.info("      (default 20) Seconds allowed to answer questions.")
        logger.info("      Note that pressing tab always stops the timer.")
        logger.info("   cluster_temp_path PATH")
        logger.info("      (default None) Allow to perform the run in PATH directory")
        logger.info("      This allow to not run on the central disk. This is not used")
        logger.info("      by condor cluster (since condor has it's own way to prevent it).")

    def help_plot(self):
        logger.info("syntax: help [RUN] [%s] [-f]" % '|'.join(self._plot_mode))
        logger.info("-- create the plot for the RUN (current run by default)")
        logger.info("     at the different stage of the event generation")
        logger.info("     Note than more than one mode can be specified in the same command.")
        logger.info("   This require to have MadAnalysis and td require. By default")
        logger.info("     if those programs are installed correctly, the creation")
        logger.info("     will be performed automaticaly during the event generation.")
        logger.info("   -f options: answer all question by default.")

    def help_compute_widths(self):
        logger.info("syntax: compute_widths Particle [Particles] [OPTIONS]")
        logger.info("-- Compute the widths for the particles specified.")
        logger.info("   By default, this takes the current param_card and overwrites it.") 
        logger.info("   Precision allows to define when to include three/four/... body decays (LO).")
        logger.info("   If this number is an integer then all N-body decay will be included.")
        logger.info("  Various options:\n")
        logger.info("  --body_decay=X: Parameter to control the precision of the computation")
        logger.info("        if X is an integer, we compute all channels up to X-body decay.")
        logger.info("        if X <1, then we stop when the estimated error is lower than X.")
        logger.info("        if X >1 BUT not an integer, then we X = N + M, with M <1 and N an integer")
        logger.info("              We then either stop at the N-body decay or when the estimated error is lower than M.")
        logger.info("        default: 4.0025")
        logger.info("  --min_br=X: All channel which are estimated below this value will not be integrated numerically.")
        logger.info("        default: precision (decimal part of the body_decay options) divided by four")
        logger.info("  --precision_channel=X: requested numerical precision for each channel")
        logger.info("        default: 0.01")
        logger.info("  --path=X: path for param_card")
        logger.info("        default: take value from the model")
        logger.info("  --output=X: path where to write the resulting card. ")
        logger.info("        default: overwrite input file. If no input file, write it in the model directory")
        logger.info("  --nlo: Compute NLO width [if the model support it]")

    def help_shower(self):
        logger.info("syntax: shower [shower_name] [shower_options]")
        logger.info("-- This is equivalent to running '[shower_name] [shower_options]'")

    def help_pgs(self):
        logger.info("syntax: pgs [RUN] [--run_options]")
        logger.info("-- run pgs on RUN (current one by default)")
        self.run_options_help([('-f','answer all question by default'),
                               ('--tag=', 'define the tag for the pgs run'),
                               ('--no_default', 'not run if pgs_card not present')])

    def help_delphes(self):
        logger.info("syntax: delphes [RUN] [--run_options]")
        logger.info("-- run delphes on RUN (current one by default)")
        self.run_options_help([('-f','answer all question by default'),
                               ('--tag=', 'define the tag for the delphes run'),
                               ('--no_default', 'not run if delphes_card not present')])

    def help_decay_events(self, skip_syntax=False):
        if not skip_syntax:
            logger.info("syntax: decay_events [RUN]")
        logger.info("This functionality allows for the decay of resonances")
        logger.info("in a .lhe file, keeping track of the spin correlation effets.")
        logger.info("BE AWARE OF THE CURRENT LIMITATIONS:")
        logger.info("  (1) Only a succession of 2 body decay are currently allowed")



class CheckValidForCmd(object):
    """ The Series of check routines in common between amcatnlo_run and
    madevent interface"""

    def check_set(self, args):
        """ check the validity of the line"""


        if len(args) < 2:
            if len(args)==1 and "=" in args[0]:
                args[:] = args[0].split("=",1)
            else:
                self.help_set()
                raise self.InvalidCmd('set needs an option and an argument')

        if args[0] not in self._set_options + self.options.keys():
            self.help_set()
            raise self.InvalidCmd('Possible options for set are %s' % \
                                  (self._set_options+self.options.keys()))

        if args[0] in ['stdout_level']:
            if args[1] not in ['DEBUG','INFO','WARNING','ERROR','CRITICAL'] \
                                                       and not args[1].isdigit():
                raise self.InvalidCmd('output_level needs ' + \
                                      'a valid level')

        if args[0] in ['timeout']:
            if not args[1].isdigit():
                raise self.InvalidCmd('timeout values should be a integer')

    def check_compute_widths(self, args):
        """check that the model is loadable and check that the format is of the
        type: PART PATH --output=PATH -f --precision=N
        return the model.
        """
        
        # Check that MG5 directory is present .
        if MADEVENT and not self.options['mg5_path']:
            raise self.InvalidCmd, '''The automatic computations of widths requires that MG5 is installed on the system.
            You can install it and set his path in ./Cards/me5_configuration.txt'''
        elif MADEVENT:
            sys.path.append(self.options['mg5_path'])
        try:
            import models.model_reader as model_reader
            import models.import_ufo as import_ufo
        except ImportError:
            raise self.ConfigurationError, '''Can\'t load MG5.
            The variable mg5_path should not be correctly configure.'''
        
        
        ufo_path = pjoin(self.me_dir,'bin','internal', 'ufomodel')
        # Import model
        if not MADEVENT:
            modelname = self.find_model_name()
            #restrict_file = None
            #if os.path.exists(pjoin(ufo_path, 'restrict_default.dat')):
            #    restrict_file = pjoin(ufo_path, 'restrict_default.dat')
            
            force_CMS = self.mother and self.mother.options['complex_mass_scheme']
            model = import_ufo.import_model(modelname, decay=True, 
                                   restrict=True, complex_mass_scheme=force_CMS)
        else:
            #pattern for checking complex mass scheme.
            has_cms = re.compile(r'''set\s+complex_mass_scheme\s*(True|T|1|true|$|;)''')
            force_CMS =  has_cms.search(open(pjoin(self.me_dir,'Cards',
                                                   'proc_card_mg5.dat')).read())
            model = import_ufo.import_model(pjoin(self.me_dir,'bin','internal',
                         'ufomodel'), decay=True, complex_mass_scheme=force_CMS)
            
#        if not hasattr(model.get('particles')[0], 'partial_widths'):
#            raise self.InvalidCmd, 'The UFO model does not include partial widths information. Impossible to compute widths automatically'
            
        # check if the name are passed to default MG5
        if '-modelname' not in open(pjoin(self.me_dir,'Cards','proc_card_mg5.dat')).read():
            model.pass_particles_name_in_mg_default()        
        model = model_reader.ModelReader(model)
        particles_name = dict([(p.get('name'), p.get('pdg_code'))
                                               for p in model.get('particles')])
        particles_name.update(dict([(p.get('antiname'), p.get('pdg_code'))
                                               for p in model.get('particles')]))        
        
        output = {'model': model, 'force': False, 'output': None, 
                  'path':None, 'particles': set(), 'body_decay':4.0025,
                  'min_br':None, 'precision_channel':0.01}
        for arg in args:
            if arg.startswith('--output='):
                output_path = arg.split('=',1)[1]
                if not os.path.exists(output_path):
                    raise self.InvalidCmd, 'Invalid Path for the output. Please retry.'
                if not os.path.isfile(output_path):
                    output_path = pjoin(output_path, 'param_card.dat')
                output['output'] = output_path       
            elif arg == '-f':
                output['force'] = True
            elif os.path.isfile(arg):
                ftype = self.detect_card_type(arg)
                if ftype != 'param_card.dat':
                    raise self.InvalidCmd , '%s is not a valid param_card.' % arg
                output['path'] = arg
            elif arg.startswith('--path='):
                arg = arg.split('=',1)[1]
                ftype = self.detect_card_type(arg)
                if ftype != 'param_card.dat':
                    raise self.InvalidCmd , '%s is not a valid param_card.' % arg
                output['path'] = arg
            elif arg.startswith('--'):
                if "=" in arg:
                    name, value = arg.split('=',1)
                    try:
                        value = float(value)
                    except Exception:
                        raise self.InvalidCmd, '--%s requires integer or a float' % name
                    output[name[2:]] = float(value)
                elif arg == "--nlo":
                    output["nlo"] = True
            elif arg in particles_name:
                # should be a particles
                output['particles'].add(particles_name[arg])
            elif arg.isdigit() and int(arg) in particles_name.values():
                output['particles'].add(ast.literal_eval(arg))
            elif arg == 'all':
                output['particles'] = set(['all'])
            else:
                self.help_compute_widths()
                raise self.InvalidCmd, '%s is not a valid argument for compute_widths' % arg
        if self.force:
            output['force'] = True

        if not output['particles']:
            raise self.InvalidCmd, '''This routines requires at least one particle in order to compute
            the related width'''
            
        if output['output'] is None:
            output['output'] = output['path']

        return output

    def check_delphes(self, arg, nodefault=False):
        """Check the argument for pythia command
        syntax: delphes [NAME] 
        Note that other option are already remove at this point
        """
        
        # If not pythia-pgs path
        if not self.options['delphes_path']:
            logger.info('Retry to read configuration file to find delphes path')
            self.set_configuration()
      
        if not self.options['delphes_path']:
            error_msg = 'No valid Delphes path set.\n'
            error_msg += 'Please use the set command to define the path and retry.\n'
            error_msg += 'You can also define it in the configuration file.\n'
            raise self.InvalidCmd(error_msg)  

        tag = [a for a in arg if a.startswith('--tag=')]
        if tag: 
            arg.remove(tag[0])
            tag = tag[0][6:]
            
                  
        if len(arg) == 0 and not self.run_name:
            if self.results.lastrun:
                arg.insert(0, self.results.lastrun)
            else:
                raise self.InvalidCmd('No run name currently define. Please add this information.')             
        
        if len(arg) == 1 and self.run_name == arg[0]:
            arg.pop(0)

        filepath = None        
        if not len(arg):
            prev_tag = self.set_run_name(self.run_name, tag, 'delphes')
            paths = [pjoin(self.me_dir,'Events',self.run_name, '%(tag)s_pythia_events.hep.gz'),
                     pjoin(self.me_dir,'Events',self.run_name, '%(tag)s_pythia8_events.hepmc.gz'),
                     pjoin(self.me_dir,'Events',self.run_name, '%(tag)s_pythia_events.hep'),
                     pjoin(self.me_dir,'Events',self.run_name, '%(tag)s_pythia8_events.hepmc'),
                     pjoin(self.me_dir,'Events','pythia_events.hep'),
                     pjoin(self.me_dir,'Events','pythia_events.hepmc'),
                     pjoin(self.me_dir,'Events','pythia8_events.hep.gz'),
                     pjoin(self.me_dir,'Events','pythia8_events.hepmc.gz')
                     ]
            for p in paths:
                if os.path.exists(p % {'tag': prev_tag}):
                    filepath = p % {'tag': prev_tag}
                    break
            else:
                a = raw_input("NO INPUT")          
                if nodefault:
                    return False
                else:
                    self.help_pgs()
                    raise self.InvalidCmd('''No file file pythia_events.* currently available
            Please specify a valid run_name''')
        
        if len(arg) == 1:
            prev_tag = self.set_run_name(arg[0], tag, 'delphes')
            if os.path.exists(pjoin(self.me_dir,'Events',self.run_name, '%s_pythia_events.hep.gz' % prev_tag)):            
                filepath = pjoin(self.me_dir,'Events',self.run_name, '%s_pythia_events.hep.gz' % prev_tag)
            elif os.path.exists(pjoin(self.me_dir,'Events',self.run_name, '%s_pythia8_events.hepmc.gz' % prev_tag)):
                filepath = pjoin(self.me_dir,'Events',self.run_name, '%s_pythia8_events.hepmc.gz' % prev_tag)
            elif os.path.exists(pjoin(self.me_dir,'Events',self.run_name, '%s_pythia_events.hep' % prev_tag)):            
                filepath = pjoin(self.me_dir,'Events',self.run_name, '%s_pythia_events.hep.gz' % prev_tag)
            elif os.path.exists(pjoin(self.me_dir,'Events',self.run_name, '%s_pythia8_events.hepmc' % prev_tag)):
                filepath = pjoin(self.me_dir,'Events',self.run_name, '%s_pythia8_events.hepmc.gz' % prev_tag)
            else:                
                raise self.InvalidCmd('No events file corresponding to %s run with tag %s.:%s '\
                    % (self.run_name, prev_tag, 
                       pjoin(self.me_dir,'Events',self.run_name, '%s_pythia_events.hep.gz' % prev_tag)))
        else:
            if tag:
                self.run_card['run_tag'] = tag
            self.set_run_name(self.run_name, tag, 'delphes')
            
        return filepath               


    
    



    def check_open(self, args):
        """ check the validity of the line """

        if len(args) != 1:
            self.help_open()
            raise self.InvalidCmd('OPEN command requires exactly one argument')

        if args[0].startswith('./'):
            if not os.path.isfile(args[0]):
                raise self.InvalidCmd('%s: not such file' % args[0])
            return True

        # if special : create the path.
        if not self.me_dir:
            if not os.path.isfile(args[0]):
                self.help_open()
                raise self.InvalidCmd('No MadEvent path defined. Unable to associate this name to a file')
            else:
                return True

        path = self.me_dir
        if os.path.isfile(os.path.join(path,args[0])):
            args[0] = os.path.join(path,args[0])
        elif os.path.isfile(os.path.join(path,'Cards',args[0])):
            args[0] = os.path.join(path,'Cards',args[0])
        elif os.path.isfile(os.path.join(path,'HTML',args[0])):
            args[0] = os.path.join(path,'HTML',args[0])
        # special for card with _default define: copy the default and open it
        elif '_card.dat' in args[0]:
            name = args[0].replace('_card.dat','_card_default.dat')
            if os.path.isfile(os.path.join(path,'Cards', name)):
                files.cp(os.path.join(path,'Cards', name), os.path.join(path,'Cards', args[0]))
                args[0] = os.path.join(path,'Cards', args[0])
            else:
                raise self.InvalidCmd('No default path for this file')
        elif not os.path.isfile(args[0]):
            raise self.InvalidCmd('No default path for this file')

    def check_treatcards(self, args):
        """check that treatcards arguments are valid
           [param|run|all] [--output_dir=] [--param_card=] [--run_card=]
        """

        opt = {'output_dir':pjoin(self.me_dir,'Source'),
               'param_card':pjoin(self.me_dir,'Cards','param_card.dat'),
               'run_card':pjoin(self.me_dir,'Cards','run_card.dat')}
        mode = 'all'
        for arg in args:
            if arg.startswith('--') and '=' in arg:
                key,value =arg[2:].split('=',1)
                if not key in opt:
                    self.help_treatcards()
                    raise self.InvalidCmd('Invalid option for treatcards command:%s ' \
                                          % key)
                if key in ['param_card', 'run_card']:
                    if os.path.isfile(value):
                        card_name = self.detect_card_type(value)
                        if card_name != key:
                            raise self.InvalidCmd('Format for input file detected as %s while expecting %s'
                                                  % (card_name, key))
                        opt[key] = value
                    elif os.path.isfile(pjoin(self.me_dir,value)):
                        card_name = self.detect_card_type(pjoin(self.me_dir,value))
                        if card_name != key:
                            raise self.InvalidCmd('Format for input file detected as %s while expecting %s'
                                                  % (card_name, key))
                        opt[key] = value
                    else:
                        raise self.InvalidCmd('No such file: %s ' % value)
                elif key in ['output_dir']:
                    if os.path.isdir(value):
                        opt[key] = value
                    elif os.path.isdir(pjoin(self.me_dir,value)):
                        opt[key] = pjoin(self.me_dir, value)
                    else:
                        raise self.InvalidCmd('No such directory: %s' % value)
            elif arg in ['MadLoop','param','run','all']:
                mode = arg
            else:
                self.help_treatcards()
                raise self.InvalidCmd('Unvalid argument %s' % arg)

        return mode, opt

    def check_decay_events(self,args):
        """Check the argument for decay_events command
        syntax is "decay_events [NAME]"
        Note that other option are already remove at this point
        """

        opts = []
        if '-from_cards' in args:
            args.remove('-from_cards')
            opts.append('-from_cards')

        if len(args) == 0:
            if self.run_name:
                args.insert(0, self.run_name)
            elif self.results.lastrun:
                args.insert(0, self.results.lastrun)
            else:
                raise self.InvalidCmd('No run name currently defined. Please add this information.')
                return

        if args[0] != self.run_name:
            self.set_run_name(args[0])

        args[0] = self.get_events_path(args[0])

        args += opts

    def check_check_events(self,args):
        """Check the argument for decay_events command
        syntax is "decay_events [NAME]"
        Note that other option are already remove at this point
        """

        if len(args) == 0:
            if self.run_name:
                args.insert(0, self.run_name)
            elif self.results.lastrun:
                args.insert(0, self.results.lastrun)
            else:
                raise self.InvalidCmd('No run name currently defined. Please add this information.')
                return
        
        if args[0] and os.path.isfile(args[0]):
            pass
        else:
            if args[0] != self.run_name:
                self.set_run_name(args[0], allow_new_tag=False)
    
            args[0] = self.get_events_path(args[0])


    def get_events_path(self, run_name):
        """return the path to the output events
        """


        if self.mode == 'madevent':
            possible_path = [
                pjoin(self.me_dir,'Events', run_name, 'unweighted_events.lhe.gz'),
                pjoin(self.me_dir,'Events', run_name, 'unweighted_events.lhe')]
        else:
            possible_path = [
                           pjoin(self.me_dir,'Events', run_name, 'events.lhe.gz'),
                           pjoin(self.me_dir,'Events', run_name, 'events.lhe')]

        for path in possible_path:
            if os.path.exists(path):
                correct_path = path
                break
        else:
            raise self.InvalidCmd('No events file corresponding to %s run. ' % run_name)
        return correct_path



class MadEventAlreadyRunning(InvalidCmd):
    pass
class AlreadyRunning(MadEventAlreadyRunning):
    pass

#===============================================================================
# CommonRunCmd
#===============================================================================
class CommonRunCmd(HelpToCmd, CheckValidForCmd, cmd.Cmd):

    debug_output = 'ME5_debug'
    helporder = ['Main Commands', 'Documented commands', 'Require MG5 directory',
                   'Advanced commands']

    # The three options categories are treated on a different footage when a
    # set/save configuration occur. current value are kept in self.options
    options_configuration = {'pythia8_path': './pythia8',
                       'hwpp_path': './herwigPP',
                       'thepeg_path': './thepeg',
                       'hepmc_path': './hepmc',
                       'madanalysis_path': './MadAnalysis',
                       'madanalysis5_path': './HEPTools/madanalysis5',
                       'pythia-pgs_path':'./pythia-pgs',
                       'td_path':'./td',
                       'delphes_path':'./Delphes',
                       'exrootanalysis_path':'./ExRootAnalysis',
                       'syscalc_path': './SysCalc',
                       'lhapdf': 'lhapdf-config',
                       'timeout': 60,
                       'f2py_compiler':None,
                       'web_browser':None,
                       'eps_viewer':None,
                       'text_editor':None,
                       'fortran_compiler':None,
                       'cpp_compiler': None,
                       'auto_update':7,
                       'cluster_type': 'condor',
                       'cluster_status_update': (600, 30),
                       'cluster_nb_retry':1,
                       'cluster_local_path': "/cvmfs/cp3.uclouvain.be/madgraph/",
                       'cluster_retry_wait':300}

    options_madgraph= {'stdout_level':None}

    options_madevent = {'automatic_html_opening':True,
                        'notification_center':True,
                         'run_mode':2,
                         'cluster_queue':None,
                         'cluster_time':None,
                         'cluster_size':100,
                         'cluster_memory':None,
                         'nb_core': None,
                         'cluster_temp_path':None}


    def __init__(self, me_dir, options, *args, **opts):
        """common"""

        self.force_run = False # this flag force the run even if RunWeb is present
        if 'force_run' in opts and opts['force_run']:
            self.force_run = True
            del opts['force_run']

        cmd.Cmd.__init__(self, *args, **opts)
        # Define current MadEvent directory
        if me_dir is None and MADEVENT:
            me_dir = root_path

        self.me_dir = me_dir
        self.options = options
        
        self.param_card_iterator = [] #an placeholder containing a generator of paramcard for scanning

        # usefull shortcut
        self.status = pjoin(self.me_dir, 'status')
        self.error =  pjoin(self.me_dir, 'error')
        self.dirbin = pjoin(self.me_dir, 'bin', 'internal')

        # Check that the directory is not currently running_in_idle
        if not self.force_run:
            if os.path.exists(pjoin(me_dir,'RunWeb')): 
                message = '''Another instance of the program is currently running.
                (for this exact same directory) Please wait that this is instance is 
                closed. If no instance is running, you can delete the file
                %s and try again.''' % pjoin(me_dir,'RunWeb')
                raise AlreadyRunning, message
            else:
                pid = os.getpid()
                fsock = open(pjoin(me_dir,'RunWeb'),'w')
                fsock.write(`pid`)
                fsock.close()
    
                misc.Popen([os.path.relpath(pjoin(self.dirbin, 'gen_cardhtml-pl'), me_dir)],
                            cwd=me_dir)

        self.to_store = []
        self.run_name = None
        self.run_tag = None
        self.banner = None
        # Load the configuration file
        self.set_configuration()
        self.configure_run_mode(self.options['run_mode'])

        # Define self.proc_characteristics
        self.get_characteristics()
        
        if not  self.proc_characteristics['ninitial']:
            # Get number of initial states
            nexternal = open(pjoin(self.me_dir,'Source','nexternal.inc')).read()
            found = re.search("PARAMETER\s*\(NINCOMING=(\d)\)", nexternal)
            self.ninitial = int(found.group(1))
        else:
            self.ninitial = self.proc_characteristics['ninitial']


    ############################################################################
    def split_arg(self, line, error=False):
        """split argument and remove run_options"""

        args = cmd.Cmd.split_arg(line)
        for arg in args[:]:
            if not arg.startswith('-'):
                continue
            elif arg == '-c':
                self.configure_run_mode(1)
            elif arg == '-m':
                self.configure_run_mode(2)
            elif arg == '-f':
                self.force = True
            elif not arg.startswith('--'):
                if error:
                    raise self.InvalidCmd('%s argument cannot start with - symbol' % arg)
                else:
                    continue
            elif arg.startswith('--cluster'):
                self.configure_run_mode(1)
            elif arg.startswith('--multicore'):
                self.configure_run_mode(2)
            elif arg.startswith('--nb_core'):
                self.options['nb_core'] = int(arg.split('=',1)[1])
                self.configure_run_mode(2)
            elif arg.startswith('--web'):
                self.pass_in_web_mode()
                self.configure_run_mode(1)
            else:
                continue
            args.remove(arg)

        return args

    @misc.multiple_try(nb_try=5, sleep=2)
    def load_results_db(self):
        """load the current results status"""
        
        # load the current status of the directory
        if os.path.exists(pjoin(self.me_dir,'HTML','results.pkl')):
            try:
                self.results = save_load_object.load_from_file(pjoin(self.me_dir,'HTML','results.pkl'))
            except Exception:
                #the pickle fail -> need to recreate the library
                model = self.find_model_name()
                process = self.process # define in find_model_name
                self.results = gen_crossxhtml.AllResults(model, process, self.me_dir)
                self.results.resetall(self.me_dir)
            else:
                try:                                
                    self.results.resetall(self.me_dir)
                except Exception, error:
                    logger.debug(error)
                    # Maybe the format was updated -> try fresh
                    model = self.find_model_name()
                    process = self.process # define in find_model_name
                    self.results = gen_crossxhtml.AllResults(model, process, self.me_dir)
                    self.results.resetall(self.me_dir)
                    self.last_mode = ''
            try:
                self.last_mode = self.results[self.results.lastrun][-1]['run_mode']
            except:
                self.results.resetall(self.me_dir)
                self.last_mode = ''

        else:
            model = self.find_model_name()
            process = self.process # define in find_model_name
            self.results = gen_crossxhtml.AllResults(model, process, self.me_dir)
            self.results.resetall(self.me_dir)
            self.last_mode=''

        return self.results

    ############################################################################
    def do_treatcards(self, line, amcatnlo=False):
        """Advanced commands: create .inc files from param_card.dat/run_card.dat"""

        keepwidth = False
        if '--keepwidth' in line:
            keepwidth = True
            line = line.replace('--keepwidth', '')
        args = self.split_arg(line)
        mode,  opt  = self.check_treatcards(args)

        if mode in ['run', 'all']:
            if not hasattr(self, 'run_card'):
                if amcatnlo:
                    run_card = banner_mod.RunCardNLO(opt['run_card'])
                else:
                    run_card = banner_mod.RunCard(opt['run_card'])
            else:
                run_card = self.run_card

            # add the conversion from the lhaid to the pdf set names
            if amcatnlo and run_card['pdlabel']=='lhapdf':
                pdfsetsdir=self.get_lhapdf_pdfsetsdir()
                pdfsets=self.get_lhapdf_pdfsets_list(pdfsetsdir)
                lhapdfsetname=[]
                for lhaid in run_card['lhaid']:
                    if lhaid in pdfsets:
                        lhapdfsetname.append(pdfsets[lhaid]['filename'])
                    else:
                        raise MadGraph5Error("lhaid %s is not a valid PDF identification number. This can be due to the use of an outdated version of LHAPDF, or %s is not a LHAGlue number corresponding to a central PDF set (but rather one of the error sets)." % (lhaid,lhaid))
                run_card['lhapdfsetname']=lhapdfsetname
            run_card.write_include_file(opt['output_dir'])

        if mode in ['MadLoop', 'all']:
            if os.path.exists(pjoin(self.me_dir, 'Cards', 'MadLoopParams.dat')):          
                self.MadLoopparam = banner_mod.MadLoopParam(pjoin(self.me_dir, 
                                                  'Cards', 'MadLoopParams.dat'))
                # write the output file
                self.MadLoopparam.write(pjoin(self.me_dir,"SubProcesses",
                                                           "MadLoopParams.dat"))

        if mode in ['param', 'all']:
            if os.path.exists(pjoin(self.me_dir, 'Source', 'MODEL', 'mp_coupl.inc')):
                param_card = check_param_card.ParamCardMP(opt['param_card'])
            else:
                param_card = check_param_card.ParamCard(opt['param_card'])
            outfile = pjoin(opt['output_dir'], 'param_card.inc')
            ident_card = pjoin(self.me_dir,'Cards','ident_card.dat')
            if os.path.isfile(pjoin(self.me_dir,'bin','internal','ufomodel','restrict_default.dat')):
                default = pjoin(self.me_dir,'bin','internal','ufomodel','restrict_default.dat')
            elif os.path.isfile(pjoin(self.me_dir,'bin','internal','ufomodel','param_card.dat')):
                default = pjoin(self.me_dir,'bin','internal','ufomodel','param_card.dat')
            elif not os.path.exists(pjoin(self.me_dir,'bin','internal','ufomodel')):
                fsock = open(pjoin(self.me_dir,'Source','param_card.inc'),'w')
                fsock.write(' ')
                fsock.close()
                return
            else:
                subprocess.call(['python', 'write_param_card.py'],
                             cwd=pjoin(self.me_dir,'bin','internal','ufomodel'))
                default = pjoin(self.me_dir,'bin','internal','ufomodel','param_card.dat')


            if amcatnlo and not keepwidth:
                # force particle in final states to have zero width
                pids = self.get_pid_final_initial_states()
                # check those which are charged under qcd
                if not MADEVENT and pjoin(self.me_dir,'bin','internal') not in sys.path:
                        sys.path.insert(0,pjoin(self.me_dir,'bin','internal'))

                #Ensure that the model that we are going to load is the current
                #one.
                to_del = [name  for name in sys.modules.keys()
                                                if name.startswith('internal.ufomodel')
                                                or name.startswith('ufomodel')]
                for name in to_del:
                    del(sys.modules[name])

                import ufomodel as ufomodel
                zero = ufomodel.parameters.ZERO
                no_width = [p for p in ufomodel.all_particles
                        if (str(p.pdg_code) in pids or str(-p.pdg_code) in pids)
                           and p.color != 1 and p.width != zero]
                done = []
                for part in no_width:
                    if abs(part.pdg_code) in done:
                        continue
                    done.append(abs(part.pdg_code))
                    param = param_card['decay'].get((part.pdg_code,))

                    if  param.value != 0:
                        logger.info('''For gauge cancellation, the width of \'%s\' has been set to zero.'''\
                                    % part.name,'$MG:color:BLACK')
                        param.value = 0

            param_card.write_inc_file(outfile, ident_card, default)

    def get_model(self):
        """return the model related to this process"""

        if self.options['mg5_path']:
            sys.path.append(self.options['mg5_path'])
            import models.import_ufo as import_ufo
            with misc.MuteLogger(['madgraph.model'],[50]):
                out= import_ufo.import_model(pjoin(self.me_dir,'bin','internal','ufomodel'))
            return out
        elif self.mother:
            return self.mother._curr_model
        else:
            return None

    def ask_edit_cards(self, cards, mode='fixed', plot=True, first_cmd=None):
        """ """
        if not self.options['madanalysis_path']:
            plot = False

        self.ask_edit_card_static(cards, mode, plot, self.options['timeout'],
                                  self.ask, first_cmd=first_cmd)

    @staticmethod
    def ask_edit_card_static(cards, mode='fixed', plot=True,
                             timeout=0, ask=None, **opt):
        if not ask:
            ask = CommonRunCmd.ask

        def path2name(path):
            if '_card' in path:
                return path.split('_card')[0]
            elif path == 'delphes_trigger.dat':
                return 'trigger'
            elif path == 'input.lhco':
                return 'lhco'
            elif path == 'MadLoopParams.dat':
                return 'MadLoopParams'
            else:
                raise Exception, 'Unknow cards name %s' % path

        # Ask the user if he wants to edit any of the files
        #First create the asking text
        question = """Do you want to edit a card (press enter to bypass editing)?\n"""
        possible_answer = ['0', 'done']
        card = {0:'done'}

        for i, card_name in enumerate(cards):
            imode = path2name(card_name)
            possible_answer.append(i+1)
            possible_answer.append(imode)
            question += '  %s / %-10s : %s\n' % (i+1, imode, card_name)
            card[i+1] = imode
        if plot:
            question += '  9 / %-10s : plot_card.dat\n' % 'plot'
            possible_answer.append(9)
            possible_answer.append('plot')
            card[9] = 'plot'

        if 'param_card.dat' in cards:
            # Add the path options
            question += ' you can also\n'
            question += '   - enter the path to a valid card or banner.\n'
            question += '   - use the \'set\' command to modify a parameter directly.\n'
            question += '     The set option works only for param_card and run_card.\n'
            question += '     Type \'help set\' for more information on this command.\n'
            question += '   - call an external program (ASperGE/MadWidth/...).\n'
            question += '     Type \'help\' for the list of available command\n'
        else:
            question += ' you can also\n'
            question += '   - enter the path to a valid card.\n'
        if 'transfer_card.dat' in cards:
            question += '   - use the \'change_tf\' command to set a transfer functions.\n'

        out = 'to_run'
        while out not in ['0', 'done']:
            out = ask(question, '0', possible_answer, timeout=int(1.5*timeout),
                              path_msg='enter path', ask_class = AskforEditCard,
                              cards=cards, mode=mode, **opt)

    @staticmethod
    def detect_card_type(path):
        """detect the type of the card. Return value are
           banner
           param_card.dat
           run_card.dat
           pythia_card.dat
           pythia8_card.dat
           plot_card.dat
           pgs_card.dat
           delphes_card.dat
           delphes_trigger.dat
           shower_card.dat [aMCatNLO]
           FO_analyse_card.dat [aMCatNLO]
           madspin_card.dat [MS]
           transfer_card.dat [MW]
           madweight_card.dat [MW]
           
           Please update the unit-test: test_card_type_recognition when adding
           cards.
        """

        fulltext = open(path).read(50000)
        if fulltext == '':
            logger.warning('File %s is empty' % path)
            return 'unknown'
        
        to_search = ['<MGVersion>',           # banner
                     '<mg5proccard>' 
                     'ParticlePropagator',    # Delphes
                     'ExecutionPath', 
                     'Treewriter', 
                     'CEN_max_tracker',
                     '#TRIGGER CARD',         # delphes_trigger.dat
                     'parameter set name',    # pgs_card
                     'muon eta coverage',
                    'req_acc_FO',
                    'MSTP',
                    'b_stable',
                    'FO_ANALYSIS_FORMAT',
                    'MSTU',
                    'Begin Minpts',
                    'gridpack',
                    'ebeam1',
                    'block\s+mw_run',
                    'BLOCK',
                    'DECAY',
                    'launch',
                    'madspin',
                    'transfer_card\.dat',
                    'set',
                    'main:numberofevents'   # pythia8  
                    ]
        
        
        text = re.findall('(%s)' % '|'.join(to_search), fulltext, re.I)
        text = [t.lower() for t in text]
        if '<mgversion>' in text or '<mg5proccard>' in text:
            return 'banner'
        elif 'particlepropagator' in text or 'executionpath' in text or 'treewriter' in text:
            return 'delphes_card.dat'
        elif 'cen_max_tracker' in text:
            return 'delphes_card.dat'
        elif '#trigger card' in text:
            return 'delphes_trigger.dat'
        elif 'parameter set name' in text:
            return 'pgs_card.dat'
        elif 'muon eta coverage' in text:
            return 'pgs_card.dat'
        elif 'mstp' in text and not 'b_stable' in text:
            return 'pythia_card.dat'
        elif 'begin minpts' in text:
            return 'plot_card.dat'
        elif ('gridpack' in text and 'ebeam1' in text) or \
                ('req_acc_fo' in text and 'ebeam1' in text):
            return 'run_card.dat'
        elif any(t.endswith('mw_run') for t in text):
            return 'madweight_card.dat'
        elif 'transfer_card.dat' in text:
            return 'transfer_card.dat'
        elif 'block' in text and 'decay' in text: 
            return 'param_card.dat'
        elif 'b_stable' in text:
            return 'shower_card.dat'
        elif 'fo_analysis_format' in text:
            return 'FO_analyse_card.dat'
        elif 'main:numberofevents' in text:
            return 'pythia8_card.dat'            
        elif 'launch' in text:
            # need to separate madspin/reweight.
            # decay/set can be in both...
            if 'madspin' in text:
                return 'madspin_card.dat'
            if 'decay' in text:
                # need to check if this a line like "decay w+" or "set decay"
                if re.search("(^|;)\s*decay", fulltext):
                    return 'madspin_card.dat'
                else:
                    return 'reweight_card.dat'
            else:
                return 'reweight_card.dat'
        else:
            return 'unknown'


    ############################################################################
    def get_available_tag(self):
        """create automatically a tag"""

        used_tags = [r['tag'] for r in self.results[self.run_name]]
        i=0
        while 1:
            i+=1
            if 'tag_%s' %i not in used_tags:
                return 'tag_%s' % i


    ############################################################################
    @misc.mute_logger(names=['madgraph.various.histograms',
                                          'internal.histograms'],levels=[20,20])
    def generate_Pythia8_HwU_plots(self, plot_root_path,
                                   merging_scale_name, observable_name, 
                                   data_path):
        """Generated the HwU plots from Pythia8 driver output for a specific
        observable."""
        
        # Make sure that the file is present
        if not os.path.isfile(data_path):
            return False

        # Load the HwU file.
        histos = histograms.HwUList(data_path, consider_reweights='ALL',run_id=0)
        if len(histos)==0:
            return False

        # Now also plot the max vs min merging scale
        merging_scales_available = [label[1] for label in \
                  histos[0].bins.weight_labels if 
                  histograms.HwU.get_HwU_wgt_label_type(label)=='merging_scale']
        if len(merging_scales_available)>=2:
            min_merging_scale = min(merging_scales_available)
            max_merging_scale = max(merging_scales_available)
        else:
            min_merging_scale = None
            max_merging_scale = None

        # jet_samples_to_keep = None means that all jet_samples are kept
        histo_output_options = {
          'format':'gnuplot', 
          'uncertainties':['scale','pdf','statistical',
                           'merging_scale','alpsfact'], 
          'ratio_correlations':True,
          'arg_string':'Automatic plotting from MG5aMC', 
          'jet_samples_to_keep':None,
          'use_band':['merging_scale','alpsfact'],
          'auto_open':False
        }
        # alpsfact variation only applies to MLM
        if not (int(self.run_card['ickkw'])==1):
            histo_output_options['uncertainties'].pop(
                histo_output_options['uncertainties'].index('alpsfact'))
            histo_output_options['use_band'].pop(
                     histo_output_options['use_band'].index('alpsfact'))

        histos.output(pjoin(plot_root_path,
            'central_%s_%s_plots'%(merging_scale_name,observable_name)),
            **histo_output_options)
        
        for scale in merging_scales_available:
            that_scale_histos = histograms.HwUList(
                               data_path,  run_id=0, merging_scale=scale)
            that_scale_histos.output(pjoin(plot_root_path,
                '%s_%.3g_%s_plots'%(merging_scale_name,scale,observable_name)),
                **histo_output_options)

        # If several merging scales were specified, then it is interesting
        # to compare the summed jet samples for the maximum and minimum
        # merging scale available.
        if not min_merging_scale is None:
            min_scale_histos = histograms.HwUList(data_path, 
                               consider_reweights=[], run_id=0, 
                                        merging_scale=min_merging_scale)
            max_scale_histos = histograms.HwUList(data_path, 
                               consider_reweights=[], run_id=0, 
                                        merging_scale=max_merging_scale)

            # Give the histos types so that the plot labels look good
            for histo in min_scale_histos:
                if histo.type is None:
                    histo.type = '%s=%.4g'%(merging_scale_name, min_merging_scale)
                else:
                    histo.type += '|%s=%.4g'%(merging_scale_name, min_merging_scale)
            for histo in max_scale_histos:
                if histo.type is None:
                    histo.type = '%s=%.4g'%(merging_scale_name, max_merging_scale)
                else:
                    histo.type += '|%s=%.4g'%(merging_scale_name, max_merging_scale)
            
            # Now plot and compare against oneanother the shape for the the two scales
            histograms.HwUList(min_scale_histos+max_scale_histos).output(
                pjoin(plot_root_path,'min_max_%s_%s_comparison'
                                         %(merging_scale_name,observable_name)),
                format='gnuplot', 
                uncertainties=[], 
                ratio_correlations=True,
                arg_string='Automatic plotting from MG5aMC', 
                jet_samples_to_keep=[],
                use_band=[],
                auto_open=False)
        return True
    
    def create_plot(self, mode='parton', event_path=None, output=None, tag=None):
        """create the plot"""

        if not tag:
            tag = self.run_card['run_tag']

        if mode != 'Pythia8':
            madir = self.options['madanalysis_path']
            td = self.options['td_path']
    
            if not madir or not td or \
                not os.path.exists(pjoin(self.me_dir, 'Cards', 'plot_card.dat')):
                return False
        else:
            PY8_plots_root_path = pjoin(self.me_dir,'HTML',
                                               self.run_name,'%s_PY8_plots'%tag)
            
        if 'ickkw' in self.run_card:
            if int(self.run_card['ickkw']) and mode == 'Pythia':
                self.update_status('Create matching plots for Pythia', level='pythia')
                # recover old data if none newly created
                if not os.path.exists(pjoin(self.me_dir,'Events','events.tree')):
                    misc.gunzip(pjoin(self.me_dir,'Events',
                          self.run_name, '%s_pythia_events.tree.gz' % tag), keep=True,
                               stdout=pjoin(self.me_dir,'Events','events.tree'))
                    files.mv(pjoin(self.me_dir,'Events',self.run_name, tag+'_pythia_xsecs.tree'),
                         pjoin(self.me_dir,'Events','xsecs.tree'))
    
                # Generate the matching plots
                misc.call([self.dirbin+'/create_matching_plots.sh',
                           self.run_name, tag, madir],
                                stdout = os.open(os.devnull, os.O_RDWR),
                                cwd=pjoin(self.me_dir,'Events'))
    
                #Clean output
                misc.gzip(pjoin(self.me_dir,"Events","events.tree"),
                          stdout=pjoin(self.me_dir,'Events',self.run_name, tag + '_pythia_events.tree.gz'))
                files.mv(pjoin(self.me_dir,'Events','xsecs.tree'),
                         pjoin(self.me_dir,'Events',self.run_name, tag+'_pythia_xsecs.tree'))
            
            elif mode == 'Pythia8' and (int(self.run_card['ickkw'])==1  or \
                  self.run_card['ktdurham']>0.0 or self.run_card['ptlund']>0.0):
                
                self.update_status('Create matching plots for Pythia8',
                                                                level='pythia8')

                # Create the directory if not existing at this stage
                if not os.path.isdir(PY8_plots_root_path):
                    os.makedirs(PY8_plots_root_path)

                merging_scale_name = 'qCut' if int(self.run_card['ickkw'])==1 \
                                                                      else 'TMS'

                djr_path = pjoin(self.me_dir,'Events',
                                             self.run_name, '%s_djrs.dat' % tag)
                pt_path = pjoin(self.me_dir,'Events',
                                             self.run_name, '%s_pts.dat' % tag)
                for observable_name, data_path in [('djr',djr_path),
                                                   ('pt',pt_path)]:
                    if not self.generate_Pythia8_HwU_plots(
                                    PY8_plots_root_path, merging_scale_name,
                                                     observable_name,data_path):
                        return False

        if mode == 'Pythia8':
            plot_files = glob.glob(pjoin(PY8_plots_root_path,'*.gnuplot'))
            if not misc.which('gnuplot'):
                logger.warning("Install gnuplot to be able to view the plots"+\
                               " generated at :\n   "+\
                               '\n   '.join('%s.gnuplot'%p for p in plot_files))
                return True
            for plot in plot_files:
                command = ['gnuplot',plot]
                try:
                    subprocess.call(command,cwd=PY8_plots_root_path,stderr=subprocess.PIPE)
                except Exception as e:
                    logger.warning("Automatic processing of the Pythia8 "+\
                            "merging plots with gnuplot failed. Try the"+\
                            " following command by hand:\n   %s"%(' '.join(command))+\
                            "\nException was: %s"%str(e))
                    return False

            plot_files = glob.glob(pjoin(PY8_plots_root_path,'*.pdf'))
            if len(plot_files)>0:
                # Add an html page
                html = "<html>\n<head>\n<TITLE>PLOT FOR PYTHIA8</TITLE>"
                html+= '<link rel=stylesheet href="../../mgstyle.css" type="text/css">\n</head>\n<body>\n'
                html += "<h2> Plot for Pythia8 </h2>\n"
                html += '<a href=../../../crossx.html>return to summary</a><br>'
                html += "<table>\n<tr> <td> <b>Obs.</b> </td> <td> <b>Type of plot</b> </td> <td><b> PDF</b> </td> <td><b> input file</b> </td> </tr>\n"
                def sorted_plots(elem):
                    name = os.path.basename(elem[1])
                    if 'central' in name:
                        return -100
                    if 'min_max' in name:
                        return -10
                    merging_re = re.match(r'^.*_(\d+)_.*$',name)
                    if not merging_re is None:
                        return int(merging_re.group(1))
                    else:
                        return 1e10
                djr_plot_files = sorted(
                            (('DJR',p) for p in plot_files if '_djr_' in p),
                            key = sorted_plots)
                pt_plot_files = sorted(
                            (('Pt',p) for p in plot_files if '_pt_' in p),
                            key = sorted_plots)
                last_obs = None            
                for obs, one_plot in djr_plot_files+pt_plot_files:
                    if obs!=last_obs:
                        # Add a line between observables
                        html += "<tr><td></td></tr>"
                        last_obs = obs
                    name = os.path.basename(one_plot).replace('.pdf','')
                    short_name = name
                    for dummy in ['_plots','_djr','_pt']:
                        short_name = short_name.replace(dummy,'')
                    short_name = short_name.replace('_',' ')                        
                    if 'min max' in short_name:
                        short_name = "%s comparison with min/max merging scale"%obs
                    if 'central' in short_name:
                        short_name = "Merging uncertainty band around central scale"
                    html += "<tr><td>%(obs)s</td><td>%(sn)s</td><td> <a href=./%(n)s.pdf>PDF</a> </td><td> <a href=./%(n)s.HwU>HwU</a> <a href=./%(n)s.gnuplot>GNUPLOT</a> </td></tr>\n" %\
                                        {'obs':obs, 'sn': short_name, 'n': name}
                html += '</table>\n'
                html += '<a href=../../../bin/internal/plot_djrs.py> Example of code to plot the above with matplotlib </a><br><br>'
                html+='</body>\n</html>'
                ff=open(pjoin(PY8_plots_root_path, 'index.html'),'w')
                ff.write(html)
            return True

        if not event_path:
            if mode == 'parton':
                possibilities=[
                    pjoin(self.me_dir, 'Events', 'unweighted_events.lhe'),
                    pjoin(self.me_dir, 'Events', 'unweighted_events.lhe.gz'),
                    pjoin(self.me_dir, 'Events', self.run_name, 'unweighted_events.lhe'),
                    pjoin(self.me_dir, 'Events', self.run_name, 'unweighted_events.lhe.gz')]
                for event_path in possibilities:
                    if os.path.exists(event_path):
                        break
                output = pjoin(self.me_dir, 'HTML',self.run_name, 'plots_parton.html')

            elif mode == 'Pythia':
                event_path = pjoin(self.me_dir, 'Events','pythia_events.lhe')
                output = pjoin(self.me_dir, 'HTML',self.run_name,
                              'plots_pythia_%s.html' % tag)
            elif mode == 'PGS':
                event_path = pjoin(self.me_dir, 'Events', self.run_name,
                                   '%s_pgs_events.lhco' % tag)
                output = pjoin(self.me_dir, 'HTML',self.run_name,
                              'plots_pgs_%s.html' % tag)
            elif mode == 'Delphes':
                event_path = pjoin(self.me_dir, 'Events', self.run_name,'%s_delphes_events.lhco' % tag)
                output = pjoin(self.me_dir, 'HTML',self.run_name,
                              'plots_delphes_%s.html' % tag)
            elif mode == "shower":
                event_path = pjoin(self.me_dir, 'Events','pythia_events.lhe')
                output = pjoin(self.me_dir, 'HTML',self.run_name,
                              'plots_shower_%s.html' % tag)
                if not self.options['pythia-pgs_path']:
                    return
            else:
                raise self.InvalidCmd, 'Invalid mode %s' % mode
        elif mode == 'reweight' and not output:
                output = pjoin(self.me_dir, 'HTML',self.run_name,
                              'plots_%s.html' % tag)

        if not os.path.exists(event_path):
            if os.path.exists(event_path+'.gz'):
                misc.gunzip('%s.gz' % event_path)
            else:
                raise self.InvalidCmd, 'Events file %s does not exist' % event_path
        elif event_path.endswith(".gz"):
             misc.gunzip(event_path)
             event_path = event_path[:-3]

             
        self.update_status('Creating Plots for %s level' % mode, level = mode.lower())

        mode = mode.lower()
        if mode not in ['parton', 'reweight']:
            plot_dir = pjoin(self.me_dir, 'HTML', self.run_name,'plots_%s_%s' % (mode.lower(),tag))
        elif mode == 'parton':
            plot_dir = pjoin(self.me_dir, 'HTML', self.run_name,'plots_parton')
        else:
            plot_dir =pjoin(self.me_dir, 'HTML', self.run_name,'plots_%s' % (tag))

        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        files.ln(pjoin(self.me_dir, 'Cards','plot_card.dat'), plot_dir, 'ma_card.dat')

        try:
            proc = misc.Popen([os.path.join(madir, 'plot_events')],
                            stdout = open(pjoin(plot_dir, 'plot.log'),'w'),
                            stderr = subprocess.STDOUT,
                            stdin=subprocess.PIPE,
                            cwd=plot_dir)
            proc.communicate('%s\n' % event_path)
            del proc
            #proc.wait()
            misc.call(['%s/plot' % self.dirbin, madir, td],
                            stdout = open(pjoin(plot_dir, 'plot.log'),'a'),
                            stderr = subprocess.STDOUT,
                            cwd=plot_dir)

            misc.call(['%s/plot_page-pl' % self.dirbin,
                                os.path.basename(plot_dir),
                                mode],
                            stdout = open(pjoin(plot_dir, 'plot.log'),'a'),
                            stderr = subprocess.STDOUT,
                            cwd=pjoin(self.me_dir, 'HTML', self.run_name))

            shutil.move(pjoin(self.me_dir, 'HTML',self.run_name ,'plots.html'),
                                                                         output)

            logger.info("Plots for %s level generated, see %s" % \
                         (mode, output))
        except OSError, error:
            logger.error('fail to create plot: %s. Please check that MadAnalysis is correctly installed.' % error)

        self.update_status('End Plots for %s level' % mode, level = mode.lower(),
                                                                 makehtml=False)

        return True

    def run_hep2lhe(self, banner_path = None):
        """Run hep2lhe on the file Events/pythia_events.hep"""

        if not self.options['pythia-pgs_path']:
            raise self.InvalidCmd, 'No pythia-pgs path defined'

        pydir = pjoin(self.options['pythia-pgs_path'], 'src')
        eradir = self.options['exrootanalysis_path']

        # Creating LHE file
        if misc.is_executable(pjoin(pydir, 'hep2lhe')):
            self.update_status('Creating shower LHE File (for plot)', level='pythia')
            # Write the banner to the LHE file
            out = open(pjoin(self.me_dir,'Events','pythia_events.lhe'), 'w')
            #out.writelines('<LesHouchesEvents version=\"1.0\">\n')
            out.writelines('<!--\n')
            out.writelines('# Warning! Never use this file for detector studies!\n')
            out.writelines('-->\n<!--\n')
            if banner_path:
                out.writelines(open(banner_path).read().replace('<LesHouchesEvents version="1.0">',''))
            out.writelines('\n-->\n')
            out.close()

            self.cluster.launch_and_wait(self.dirbin+'/run_hep2lhe',
                                         argument= [pydir],
                                        cwd=pjoin(self.me_dir,'Events'),
                                        stdout=os.devnull)

            logger.info('Warning! Never use this lhe file for detector studies!')
            # Creating ROOT file
            if eradir and misc.is_executable(pjoin(eradir, 'ExRootLHEFConverter')):
                self.update_status('Creating Pythia LHE Root File', level='pythia')
                try:
                    misc.call([eradir+'/ExRootLHEFConverter',
                             'pythia_events.lhe',
                             pjoin(self.run_name, '%s_pythia_lhe_events.root' % self.run_tag)],
                            cwd=pjoin(self.me_dir,'Events'))
                except Exception, error:
                    misc.sprint('ExRootLHEFConverter fails', str(error),
                                                                     log=logger)
                    pass

    def store_result(self):
        """Dummy routine, to be overwritten by daughter classes"""

        pass
    
    ############################################################################
    def do_reweight(self, line):
        """ syntax is "reweight RUN_NAME"
            Allow to reweight the events generated with a new choices of model
            parameter. Description of the methods are available here
            cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Reweight
        """
        
        if '-from_cards' in line and not os.path.exists(pjoin(self.me_dir, 'Cards', 'reweight_card.dat')):
            return
        
        # Check that MG5 directory is present .
        if MADEVENT and not self.options['mg5_path']:
            raise self.InvalidCmd, '''The module reweight requires that MG5 is installed on the system.
            You can install it and set its path in ./Cards/me5_configuration.txt'''
        elif MADEVENT:
            sys.path.append(self.options['mg5_path'])
        try:
            import madgraph.interface.reweight_interface as reweight_interface
        except ImportError:
            raise self.ConfigurationError, '''Can\'t load Reweight module.
            The variable mg5_path might not be correctly configured.'''
        

                        
        if not '-from_cards' in line:
            self.keep_cards(['reweight_card.dat'], ignore=['*'])
            self.ask_edit_cards(['reweight_card.dat'], 'fixed', plot=False)        

        # load the name of the event file
        args = self.split_arg(line) 

        if not self.force_run:
            # forbid this function to create an empty item in results.
            if self.run_name and self.results.current and  self.results.current['cross'] == 0:
                self.results.delete_run(self.run_name, self.run_tag)
            self.results.save()
            # we want to run this in a separate shell to avoid hard f2py crash
            command =  [sys.executable]
            if os.path.exists(pjoin(self.me_dir, 'bin', 'madevent')):
                command.append(pjoin(self.me_dir, 'bin', 'internal','madevent_interface.py'))
            else:
                command.append(pjoin(self.me_dir, 'bin', 'internal', 'amcatnlo_run_interface.py'))
            if not isinstance(self, cmd.CmdShell):
                command.append('--web')
            command.append('reweight')
            if self.run_name:
                command.append(self.run_name)
            else:
                command += args
            if '-from_cards' not in command:
                command.append('-from_cards')
            p = misc.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, cwd=self.me_dir)
            while p.poll() is None:
                line = p.stdout.readline()
                if any(t in line for t in ['INFO:', 'WARNING:', 'CRITICAL:', 'ERROR:', 'root:','KEEP:']) and \
                   not '***********' in line:
                        print line[:-1].replace('INFO', 'REWEIGHT').replace('KEEP:','')
                elif __debug__ and line:
                    logger.debug(line[:-1])
            if p.returncode !=0:
                logger.error("Reweighting failed")
                return
            self.results = self.load_results_db()
            # forbid this function to create an empty item in results.
            try:
                if self.results[self.run_name][-2]['cross']==0:
                    self.results.delete_run(self.run_name,self.results[self.run_name][-2]['tag'])
            except:
                pass
            try:
                if self.results.current['cross'] == 0 and self.run_name:
                    self.results.delete_run(self.run_name, self.run_tag)
            except:
                pass                    
            # re-define current run     
            try:
                self.results.def_current(self.run_name, self.run_tag)
            except Exception:
                pass
            return              

        self.to_store.append('event')
        # forbid this function to create an empty item in results.
        if self.results.current['cross'] == 0 and self.run_name:
            self.results.delete_run(self.run_name, self.run_tag)

        self.check_decay_events(args) 
        # args now alway content the path to the valid files
        reweight_cmd = reweight_interface.ReweightInterface(args[0], mother=self)
        #reweight_cmd.use_rawinput = False
        #reweight_cmd.mother = self
        wgt_names = reweight_cmd.get_weight_names()
        if wgt_names == [''] and reweight_cmd.has_nlo:
            self.update_status('Running Reweighting (LO approximate)', level='madspin')
        else:
            self.update_status('Running Reweighting', level='madspin')
        
        path = pjoin(self.me_dir, 'Cards', 'reweight_card.dat')
        reweight_cmd.raw_input=False
        reweight_cmd.me_dir = self.me_dir
        reweight_cmd.import_command_file(path)
        reweight_cmd.do_quit('')
            
        logger.info("quit rwgt")
        
        
        
        # re-define current run
        try:
            self.results.def_current(self.run_name, self.run_tag)
        except Exception:
            pass

    ############################################################################
    def do_pgs(self, line):
        """launch pgs"""
        
        args = self.split_arg(line)
        # Check argument's validity
        if '--no_default' in args:
            no_default = True
            args.remove('--no_default')
        else:
            no_default = False

        if no_default and not os.path.exists(pjoin(self.me_dir, 'Cards', 'pgs_card.dat')):
            logger.info('No pgs_card detected, so not run pgs')
            return

        # Check all arguments
        # This might launch a gunzip in another thread. After the question
        # This thread need to be wait for completion. (This allow to have the
        # question right away and have the computer working in the same time)
        # if lock is define this a locker for the completion of the thread
        lock = self.check_pgs(args,  no_default=no_default)

        # Check that the pgs_card exists. If not copy the default
        if not os.path.exists(pjoin(self.me_dir, 'Cards', 'pgs_card.dat')):
            files.cp(pjoin(self.me_dir, 'Cards', 'pgs_card_default.dat'),
                     pjoin(self.me_dir, 'Cards', 'pgs_card.dat'))
            logger.info('No pgs card found. Take the default one.')

        if not (no_default or self.force):
            self.ask_edit_cards(['pgs_card.dat'])

        self.update_status('prepare PGS run', level=None)

        pgsdir = pjoin(self.options['pythia-pgs_path'], 'src')
        eradir = self.options['exrootanalysis_path']
        madir = self.options['madanalysis_path']
        td = self.options['td_path']

        # Compile pgs if not there
        if not misc.is_executable(pjoin(pgsdir, 'pgs')):
            logger.info('No PGS executable -- running make')
            misc.compile(cwd=pgsdir)

        self.update_status('Running PGS', level='pgs')

        tag = self.run_tag
        # Update the banner with the pgs card
        banner_path = pjoin(self.me_dir, 'Events', self.run_name, '%s_%s_banner.txt' % (self.run_name, self.run_tag))
        if os.path.exists(pjoin(self.me_dir, 'Source', 'banner_header.txt')):
            self.banner.add(pjoin(self.me_dir, 'Cards','pgs_card.dat'))
            self.banner.write(banner_path)
        else:
            open(banner_path, 'w').close()

        ########################################################################
        # now pass the event to a detector simulator and reconstruct objects
        ########################################################################
        if lock:
            lock.wait()
        # Prepare the output file with the banner
        ff = open(pjoin(self.me_dir, 'Events', 'pgs_events.lhco'), 'w')
        if os.path.exists(pjoin(self.me_dir, 'Source', 'banner_header.txt')):
            text = open(banner_path).read()
            text = '#%s' % text.replace('\n','\n#')
            dico = self.results[self.run_name].get_current_info()
            text +='\n##  Integrated weight (pb)  : %.4g' % dico['cross']
            text +='\n##  Number of Event         : %s\n' % dico['nb_event']
            ff.writelines(text)
        ff.close()

        try:
            os.remove(pjoin(self.me_dir, 'Events', 'pgs.done'))
        except Exception:
            pass

        pgs_log = pjoin(self.me_dir, 'Events', self.run_name, "%s_pgs.log" % tag)
        self.cluster.launch_and_wait('../bin/internal/run_pgs',
                            argument=[pgsdir], cwd=pjoin(self.me_dir,'Events'),
                            stdout=pgs_log, stderr=subprocess.STDOUT)

        if not os.path.exists(pjoin(self.me_dir, 'Events', 'pgs.done')):
            logger.error('Fail to create LHCO events')
            return
        else:
            os.remove(pjoin(self.me_dir, 'Events', 'pgs.done'))

        if os.path.getsize(banner_path) == os.path.getsize(pjoin(self.me_dir, 'Events','pgs_events.lhco')):
            misc.call(['cat pgs_uncleaned_events.lhco >>  pgs_events.lhco'],
                            cwd=pjoin(self.me_dir, 'Events'))
            os.remove(pjoin(self.me_dir, 'Events', 'pgs_uncleaned_events.lhco '))

        # Creating Root file
        if eradir and misc.is_executable(pjoin(eradir, 'ExRootLHCOlympicsConverter')):
            self.update_status('Creating PGS Root File', level='pgs')
            try:
                misc.call([eradir+'/ExRootLHCOlympicsConverter',
                             'pgs_events.lhco',pjoin('%s/%s_pgs_events.root' % (self.run_name, tag))],
                            cwd=pjoin(self.me_dir, 'Events'))
            except Exception:
                logger.warning('fail to produce Root output [problem with ExRootAnalysis')
        if os.path.exists(pjoin(self.me_dir, 'Events', 'pgs_events.lhco')):
            # Creating plots
            files.mv(pjoin(self.me_dir, 'Events', 'pgs_events.lhco'),
                    pjoin(self.me_dir, 'Events', self.run_name, '%s_pgs_events.lhco' % tag))
            self.create_plot('PGS')
            misc.gzip(pjoin(self.me_dir, 'Events', self.run_name, '%s_pgs_events.lhco' % tag))

        self.update_status('finish', level='pgs', makehtml=False)

    ############################################################################                                                                                                           
    def do_compute_widths(self, line):
        """Require MG5 directory: Compute automatically the widths of a set 
        of particles"""



        args = self.split_arg(line)
        opts = self.check_compute_widths(args)

        from madgraph.interface.master_interface import MasterCmd
        cmd = MasterCmd()
        self.define_child_cmd_interface(cmd, interface=False)
        cmd.exec_cmd('set automatic_html_opening False --no_save')
        if not opts['path']:
            opts['path'] = pjoin(self.me_dir, 'Cards', 'param_card.dat')
            if not opts['force'] :
                self.ask_edit_cards(['param_card'],[], plot=False)
        
        
        line = 'compute_widths %s %s' % \
                (' '.join([str(i) for i in opts['particles']]),
                 ' '.join('--%s=%s' % (key,value) for (key,value) in opts.items()
                        if key not in ['model', 'force', 'particles'] and value))
        cmd.exec_cmd(line, model=opts['model'])
        self.child = None
        del cmd

    ############################################################################ 
    def do_print_results(self, line):
        """Not in help:Print the cross-section/ number of events for a given run"""
        
        args = self.split_arg(line)
        options={'path':None, 'mode':'w', 'format':'full'}
        for arg in list(args):
            if arg.startswith('--') and '=' in arg:
                name,value=arg.split('=',1)
                name = name [2:]
                options[name] = value
                args.remove(arg)
        
        
        if len(args) > 0:
            run_name = args[0]
        else:
            for i, run_name in enumerate(self.results.order):
                for j, one_result in enumerate(self.results[run_name]):
                    if i or j:
                        options['mode'] = "a"
                    if options['path']:
                        self.print_results_in_file(one_result, options['path'], options['mode'], options['format'])
                    else:
                        self.print_results_in_shell(one_result)
            return

        if run_name not in self.results:
            raise self.InvalidCmd('%s is not a valid run_name or it doesn\'t have any information' \
                                  % run_name)

            
        if len(args) == 2:
            tag = args[1]
            if tag.isdigit():
                tag = int(tag) - 1
                if len(self.results[run_name]) < tag:
                    raise self.InvalidCmd('Only %s different tag available' % \
                                                    len(self.results[run_name]))
                data = self.results[run_name][tag]
            else:
                data = self.results[run_name].return_tag(tag)
        else:
            data = self.results[run_name].return_tag(None) # return the last
        
        if options['path']:
            self.print_results_in_file(data, options['path'], options['mode'], options['format'])
        else:
            self.print_results_in_shell(data)

    def configure_directory(self, *args, **opts):
        """ All action require before any type of run. Typically overloaded by
        daughters if need be."""
        pass

    ############################################################################
    # Start of MadAnalysis5 related function
    ############################################################################
    
    def check_madanalysis5(self, args, mode='parton'):
        """Check the argument for the madanalysis5 command
        syntax: madanalysis5_parton [NAME]
        """

        MA5_options = {'MA5_stdout_lvl':'default'}
        
        stdout_level_tags = [a for a in args if a.startswith('--MA5_stdout_lvl=')]
        for slt in stdout_level_tags:
            lvl = slt.split('=')[1].strip()
            try:
                # It is likely an int
                MA5_options['MA5_stdout_lvl']=int(lvl)
            except ValueError:
                try:
                    # Maybe the user used something like 'logging.INFO'
                    MA5_options['MA5_stdout_lvl']=eval(lvl)
                except:
                    try:
                        MA5_options['MA5_stdout_lvl']=eval('logging.%s'%lvl)
                    except:
                        raise InvalidCmd("MA5 output level specification"+\
                                                 " '%s' is incorrect."%str(lvl))
            args.remove(slt)

        if mode=='parton':
            # We will attempt to run MA5 on the parton level output
            # found in the last run if not specified.
            MA5_options['inputs'] = '*.lhe'
        elif mode=='hadron':
            # We will run MA5 on all sources of post-partonic output we
            # can find if not specified. PY8 is a keyword indicating shower
            # piped to MA5.
            MA5_options['inputs'] = ['fromCard']
        else:
            raise MadGraph5Error('Mode %s not reckognized'%mode+
                                             ' in function check_madanalysis5.')
        # If not madanalysis5 path
        if not self.options['madanalysis5_path']:
            logger.info('Now trying to read the configuration file again'+
                                                   ' to find MadAnalysis5 path')
            self.set_configuration()
            
        if not self.options['madanalysis5_path'] or not \
            os.path.exists(pjoin(self.options['madanalysis5_path'],'bin','ma5')):
            error_msg = 'No valid MadAnalysis5 path set.\n'
            error_msg += 'Please use the set command to define the path and retry.\n'
            error_msg += 'You can also define it in the configuration file.\n'
            error_msg += 'Finally, it can be installed automatically using the'
            error_msg += ' install command.\n'
            raise self.InvalidCmd(error_msg)

        # Now make sure that the corresponding default card exists
        if not os.path.isfile(pjoin(self.me_dir,
                               'Cards','madanalysis5_%s_card.dat'%mode)):
            raise self.InvalidCmd('Your installed version of MadAnalysis5 and/or'+\
                    ' MadGraph5_aMCatNLO does not seem to support analysis at'+
                                                            '%s level.'%mode)
        
        tag = [a for a in args if a.startswith('--tag=')]
        if tag: 
            args.remove(tag[0])
            tag = tag[0][6:]

        if len(args) == 0 and not self.run_name:
            if self.results.lastrun:
                args.insert(0, self.results.lastrun)
            else:
                raise self.InvalidCmd('No run name currently defined. '+
                                                 'Please add this information.')
        
        if len(args) >= 1:
            if mode=='parton' and args[0] != self.run_name and \
             not os.path.exists(pjoin(self.me_dir,'Events',args[0], 
             'unweighted_events.lhe.gz')) and not os.path.exists(
                                           pjoin(self.me_dir,'Events',args[0])):
                raise self.InvalidCmd('No events file in the %s run.'%args[0])
            self.set_run_name(args[0], tag, level='madanalysis5_%s'%mode)            
        else:
            if tag:
                self.run_card['run_tag'] = args[0]
            self.set_run_name(self.run_name, tag, level='madanalysis5_%s'%mode)  
        
        if mode=='parton':
            if any(t for t in args if t.startswith('--input=')):
                raise InvalidCmd('The option --input=<input_file> is not'+
                  ' available when running partonic MadAnalysis5 analysis. The'+
                      ' .lhe output of the selected run is used automatically.')
            input_file = pjoin(self.me_dir,'Events',self.run_name, 'unweighted_events.lhe')
            MA5_options['inputs'] = '%s.gz'%input_file
            if not os.path.exists('%s.gz'%input_file):
                if os.path.exists(input_file):
                    misc.gzip(input_file, keep=True, stdout=output_file)
                else:
                    logger.warning("LHE event file not found in \n%s\ns"%input_file+
                                       "AboringParton level MA5 analysis is disabled.")                    
    
        if mode=='hadron':
            # Make sure to store current results (like Pythia8 hep files)
            # so that can be found here
            self.store_result()
            
            hadron_tag = [t for t in args if t.startswith('--input=')]
            if hadron_tag and hadron_tag[0][8:]:
                hadron_inputs = hadron_tag[0][8:].split(',')
            
            # If not set above, then we must read it from the card
            elif MA5_options['inputs'] == ['fromCard']:
                hadron_inputs = banner_mod.MadAnalysis5Card(pjoin(self.me_dir,
                'Cards','madanalysis5_hadron_card.dat'),mode='hadron')['inputs']

            # Make sure the corresponding input files are present and unfold
            # potential wildcard while making their path absolute as well.
            MA5_options['inputs'] = []
            special_source_tags = []
            for htag in hadron_inputs:
                # Possible pecial tag for MA5 run inputs
                if htag in special_source_tags:
                    # Special check/actions
                    continue
                # Check if the specified file exists and is not a wildcard
                if os.path.isfile(htag) or (os.path.exists(htag) and 
                                          stat.S_ISFIFO(os.stat(htag).st_mode)):
                    MA5_options['inputs'].append(htag)
                    continue

                # Now select one source per tag, giving priority to unzipped 
                # files with 'events' in their name (case-insensitive).
                file_candidates = glob.glob(pjoin(self.me_dir,'Events',self.run_name,htag))+\
                            glob.glob(pjoin(self.me_dir,'Events',self.run_name,'%s.gz'%htag))
                priority_files = [f for f in file_candidates if 
                                self.run_card['run_tag'] in os.path.basename(f)]
                priority_files = [f for f in priority_files if
                                        'EVENTS' in os.path.basename(f).upper()]
                # Make sure to always prefer the original partonic event file
                for f in file_candidates:
                    if os.path.basename(f).startswith('unweighted_events.lhe'):
                        priority_files.append(f)
                if priority_files:
                    MA5_options['inputs'].append(priority_files[-1])
                    continue
                if file_candidates:
                    MA5_options['inputs'].append(file_candidates[-1])
                    continue

        return MA5_options
    
    def ask_madanalysis5_run_configuration(self, runtype='parton',mode=None):
        """Ask the question when launching madanalysis5.
        In the future we can ask here further question about the MA5 run, but
        for now we just edit the cards"""

        cards = ['madanalysis5_%s_card.dat'%runtype]
        self.keep_cards(cards)
        
        if self.force:
            return runtype
        
        # This heavy-looking structure of auto is just to mimick what is done
        # for ask_pythia_configuration
        auto=False
        if mode=='auto':
            auto=True
        if auto:
            self.ask_edit_cards(cards, mode='auto', plot=False)
        else:
            self.ask_edit_cards(cards, plot=False)

        # For now, we don't pass any further information and simply return the
        # input mode asked for
        mode = runtype 
        return mode

    def complete_madanalysis5_hadron(self,text, line, begidx, endidx):
        "Complete the madanalysis5 command"
        args = self.split_arg(line[0:begidx], error=False)
        if len(args) == 1:
            #return valid run_name
            data = []
            for name in banner_mod.MadAnalysis5Card._default_hadron_inputs:
                data += glob.glob(pjoin(self.me_dir, 'Events', '*','%s'%name))
                data += glob.glob(pjoin(self.me_dir, 'Events', '*','%s.gz'%name))
            data = [n.rsplit('/',2)[1] for n in data]
            tmp1 =  self.list_completion(text, data)
            if not self.run_name:
                return tmp1
            else:
                tmp2 = self.list_completion(text, ['-f',
                '--MA5_stdout_lvl=','--input=','--no_default', '--tag='], line)
                return tmp1 + tmp2
            
        elif '--MA5_stdout_lvl=' in line and not any(arg.startswith(
                                          '--MA5_stdout_lvl=') for arg in args):
            return self.list_completion(text, 
                ['--MA5_stdout_lvl=%s'%opt for opt in 
                ['logging.INFO','logging.DEBUG','logging.WARNING',
                                                'logging.CRITICAL','90']], line)
        elif '--input=' in line and not any(arg.startswith(
                                                  '--input=') for arg in args):
            return self.list_completion(text, ['--input=%s'%opt for opt in
         (banner_mod.MadAnalysis5Card._default_hadron_inputs +['path'])], line)
        else:
            return self.list_completion(text, ['-f', 
                '--MA5_stdout_lvl=','--input=','--no_default', '--tag='], line)

    def do_madanalysis5_hadron(self, line):
        """launch MadAnalysis5 at the hadron level."""
        return self.run_madanalysis5(line,mode='hadron')

    def run_madanalysis5(self, line, mode='parton'):
        """launch MadAnalysis5 at the parton level or at the hadron level with
        a specific command line."""
        
        # Check argument's validity
        args = self.split_arg(line)
        
        if '--no_default' in args:
            no_default = True
            args.remove('--no_default')
        else:
            no_default = False

        if no_default:
            # Called issued by MG5aMC itself during a generate_event action
            if mode=='parton' and not os.path.exists(pjoin(self.me_dir, 'Cards',
                                               'madanalysis5_parton_card.dat')):
                return
            if mode=='hadron' and not os.path.exists(pjoin(self.me_dir, 'Cards',
                                               'madanalysis5_hadron_card.dat')):
                return
        else:
            # Called issued by the user itself and only MA5 will be run.
            # we must therefore ask wheter the user wants to edit the card
            self.ask_madanalysis5_run_configuration(runtype=mode)        
        
        if not self.options['madanalysis5_path'] or \
            all(not os.path.exists(pjoin(self.me_dir, 'Cards',card)) for card in
               ['madanalysis5_parton_card.dat','madanalysis5_hadron_card.dat']):
            if no_default:
                return
            else:
                raise InvalidCmd('You must have MadAnalysis5 available to run'+
           " this command. Consider installing it with the 'install' function.")

        if not self.run_name:
            MA5_opts = self.check_madanalysis5(args, mode=mode)
            self.configure_directory(html_opening =False)
        else:
            # initialize / remove lhapdf mode        
            self.configure_directory(html_opening =False)
            MA5_opts = self.check_madanalysis5(args, mode=mode)

        # Now check that there is at least one input to run
        if MA5_opts['inputs']==[]:
            if no_default:
                logger.warning('No hadron level input found to run MadAnalysis5 on.'+
                                         ' Skipping its hadron-level analysis.')
                return
            else:
                raise self.InvalidCmd('\nNo input files specified or availabled for'+
        ' this MadAnalysis5 hadron-level run.\nPlease double-check the options of this'+
        ' MA5 command (or card) and which output files\nare currently in the chosen'+
        " run directory '%s'."%self.run_name)

        MA5_card = banner_mod.MadAnalysis5Card(pjoin(self.me_dir, 'Cards',
                                    'madanalysis5_%s_card.dat'%mode), mode=mode)
        
        MA5_cmds_list = MA5_card.get_MA5_cmds(MA5_opts['inputs'],
                pjoin(self.me_dir,'MA5_%s_ANALYSIS'%mode.upper()),
                run_dir_path = pjoin(self.me_dir,'Events', self.run_name),
                UFO_model_path=pjoin(self.me_dir,'bin','internal','ufomodel'),
                run_tag = self.run_tag)

############################# TEMPORARY ########################################
#       Here's how to print the MA5 commands generated by MG5aMC
#        for MA5_runtag, MA5_cmds in MA5_cmds_list:
#            misc.sprint('****************************************')
#            misc.sprint('* Commands for MA5 runtag %s:'%MA5_runtag)
#            misc.sprint('\n'+('\n'.join('* %s'%cmd for cmd in MA5_cmds)))
#            misc.sprint('****************************************')
############################# TEMPORARY ########################################
        
        logger.info('-------------------------------------------','$MG:color:GREEN')
        logger.info('Now starting MadAnalysis5 [arXiv:1206.1599]','$MG:color:GREEN')
        logger.info('-------------------------------------------','$MG:color:GREEN')
        if mode=='hadron':
            logger.info('Hadron input files considered:')
            for input in MA5_opts['inputs']:
                logger.info('  --> %s'%input)
        elif mode=='parton':
            logger.info('Parton input file considered:')
            logger.info('  --> %s'%MA5_opts['inputs'])
        logger.info('-------------------------------------------','$MG:color:GREEN')

        # Obtain a main MA5 interpreter
        # Ideally we would like to do it all with a single interpreter
        # but we'd need a way to reset it for this.
        if MA5_opts['MA5_stdout_lvl']=='default':
            if MA5_card['stdout_lvl'] is None:
                MA5_lvl = self.options['stdout_level']
            else:
                MA5_lvl = MA5_card['stdout_lvl']                
        else:
            MA5_lvl = MA5_opts['MA5_stdout_lvl']

        MA5_interpreter = misc.get_MadAnalysis5_interpreter(
                self.options['mg5_path'], 
                self.options['madanalysis5_path'],
                logstream=sys.stdout,
                loglevel=MA5_lvl,
                forced=True)

        # Now loop over the different MA5_runs
        for MA5_runtag, MA5_cmds in MA5_cmds_list:
            
            # Make sure to properly initialize MA5 interpreter
            if mode=='hadron':
                MA5_interpreter.init_reco()
            else:
                MA5_interpreter.init_parton()
            
            if MA5_runtag!='default':
                if MA5_runtag.startswith('_reco_'):
                    logger.info("MadAnalysis5 now running the reconstruction '%s'..."%
                                              MA5_runtag[6:],'$MG:color:GREEN')
                elif MA5_runtag=='Recasting':
                    logger.info("MadAnalysis5 now running the recasting...",
                                                              '$MG:color:GREEN') 
                else:
                    logger.info("MadAnalysis5 now running the '%s' analysis..."%
                                                   MA5_runtag,'$MG:color:GREEN')                    
            # Now the magic, finally call MA5.
            MA5_interpreter.load(MA5_cmds)
                     
            if MA5_runtag.startswith('_reco_'):
                # When doing a reconstruction we must first copy the event file
                # created and then directly proceed to the next batch of instructions
                # There can be several output directory if there were several
                # input files.
                links_created=[]
                for i, input in enumerate(MA5_opts['inputs']):
                    # Make sure it is not an lhco or root input, which would not
                    # undergo any reconstruction of course.
                    if not banner_mod.MadAnalysis5Card.events_can_be_reconstructed(input):
                        continue
                    reco_output = pjoin(self.me_dir,
                           'MA5_%s_ANALYSIS%s_%d'%(mode.upper(),MA5_runtag,i+1))
                    # Look for either a root or .lhe.gz output
                    reco_event_file = glob.glob(pjoin(reco_output,'Output',
                      '_reco_events','*.lhe.gz')) + glob.glob(pjoin(
                                  reco_output,'Output','_reco_events','*.root'))
                    if len(reco_event_file)==0:
                        raise MadGraph5Error, "MadAnalysis5 failed to produce the "+\
                  "reconstructed event file for reconstruction '%s'."%MA5_runtag[6:]
                    reco_event_file = reco_event_file[0]
                    # move the reconstruction output to the HTML directory
                    shutil.move(reco_output,pjoin(self.me_dir,'HTML',
                                 self.run_name,'%s_MA5_%s_ANALYSIS%s_%d'%
                                    (self.run_tag,mode.upper(),MA5_runtag,i+1)))
                    # link the reconstructed event file to the run directory
                    links_created.append(os.path.basename(reco_event_file))
                    files.ln(pjoin(self.me_dir,'HTML',self.run_name,
                      '%s_MA5_%s_ANALYSIS%s_%d'%(self.run_tag,mode.upper(),
                      MA5_runtag,i+1),'Output','_reco_events',links_created[-1]),
                                      pjoin(self.me_dir,'Events',self.run_name))
                    
                logger.info("MadAnalysis5 successfully completed the reconstruction "+
                  "'%s'. Links to the reconstructed event files are:"%MA5_runtag[6:])
                for link in links_created:
                    logger.info('  --> %s'%pjoin(self.me_dir,'Events',self.run_name,link))
                continue

            if MA5_runtag.upper()=='RECASTING':
                target = pjoin(self.me_dir,'MA5_%s_ANALYSIS_%s'\
              %(mode.upper(),MA5_runtag),'Output','CLs_output_summary.dat')
            else:
                target = pjoin(self.me_dir,'MA5_%s_ANALYSIS_%s'\
                                  %(mode.upper(),MA5_runtag),'PDF','main.pdf')
            if not os.path.isfile(target):
                raise MadGraph5Error, "MadAnalysis5 failed to produced "+\
          "an output for the parton analysis '%s' in\n   %s"%(MA5_runtag,target)

            # Copy the PDF report or CLs in the Events/run directory.
            if MA5_runtag.upper()=='RECASTING':
                carboncopy_name = '%s_MA5_CLs.dat'%(self.run_tag)
            else:
                carboncopy_name = '%s_MA5_%s_analysis_%s.pdf'%(
                                                   self.run_tag,mode,MA5_runtag)
            shutil.copy(target, pjoin(self.me_dir,'Events',self.run_name,carboncopy_name))
            if MA5_runtag!='default':
                logger.info("MadAnalysis5 successfully completed the "+
                  "%s. Reported results are placed in:"%("analysis '%s'"%MA5_runtag 
                           if MA5_runtag.upper()!='RECASTING' else "recasting"))
            else:
                logger.info("MadAnalysis5 successfully completed the analysis."+
                                            " Reported results are placed in:")
            logger.info('  --> %s'%pjoin(self.me_dir,'Events',self.run_name,carboncopy_name))

            # Copy the entire analysis in the HTML directory
            shutil.move(pjoin(self.me_dir,'MA5_%s_ANALYSIS_%s'\
              %(mode.upper(),MA5_runtag)), pjoin(self.me_dir,'HTML',self.run_name,
                '%s_MA5_%s_ANALYSIS_%s'%(self.run_tag,mode.upper(),MA5_runtag)))

        # Set the number of events and cross-section to the last one 
        # (maybe do something smarter later)
        new_details={}
        for detail in ['nb_event','cross','error']:
            new_details[detail] = \
                      self.results[self.run_name].get_current_info()[detail]
        for detail in new_details:
            self.results.add_detail(detail,new_details[detail])

        self.update_status('finish', level='madanalysis5_%s'%mode,
                                                                 makehtml=False)
 
        #Update the banner
        self.banner.add(pjoin(self.me_dir, 'Cards',
                                               'madanalysis5_%s_card.dat'%mode))
        banner_path = pjoin(self.me_dir,'Events', self.run_name,
                               '%s_%s_banner.txt'%(self.run_name, self.run_tag))
        self.banner.write(banner_path)
 
        if not no_default:
            logger.info('Find more information about this run on the HTML local page')
            logger.info('  --> %s'%pjoin(self.me_dir,'index.html'))
    
    ############################################################################
    # End of MadAnalysis5 related function
    ############################################################################
    
    def do_delphes(self, line):
        """ run delphes and make associate root file/plot """

        args = self.split_arg(line)
        # Check argument's validity
        if '--no_default' in args:
            no_default = True
            args.remove('--no_default')
        else:
            no_default = False
            
        if no_default and  not os.path.exists(pjoin(self.me_dir, 'Cards', 'delphes_card.dat')):
            logger.info('No delphes_card detected, so not run Delphes')
            return
            
        # Check all arguments
        filepath = self.check_delphes(args, nodefault=no_default)
        if no_default and not filepath:
            return # no output file but nothing to do either.
        
        self.update_status('prepare delphes run', level=None)

        if os.path.exists(pjoin(self.options['delphes_path'], 'data')):
            delphes3 = False
            prog = '../bin/internal/run_delphes'
            if filepath and '.hepmc' in filepath[:-10]:
                raise self.InvalidCmd, 'delphes2 do not support hepmc'
        else:
            delphes3 = True
            prog =  '../bin/internal/run_delphes3'

        # Check that the delphes_card exists. If not copy the default and
        # ask for edition of the card.
        if not os.path.exists(pjoin(self.me_dir, 'Cards', 'delphes_card.dat')):
            if no_default:
                logger.info('No delphes_card detected, so not running Delphes')
                return
            files.cp(pjoin(self.me_dir, 'Cards', 'delphes_card_default.dat'),
                     pjoin(self.me_dir, 'Cards', 'delphes_card.dat'))
            logger.info('No delphes card found. Take the default one.')
        if not delphes3 and not os.path.exists(pjoin(self.me_dir, 'Cards', 'delphes_trigger.dat')):
            files.cp(pjoin(self.me_dir, 'Cards', 'delphes_trigger_default.dat'),
                     pjoin(self.me_dir, 'Cards', 'delphes_trigger.dat'))
        if not (no_default or self.force):
            if delphes3:
                self.ask_edit_cards(['delphes_card.dat'], args)
            else:
                self.ask_edit_cards(['delphes_card.dat', 'delphes_trigger.dat'], args)

        self.update_status('Running Delphes', level=None)

        delphes_dir = self.options['delphes_path']
        tag = self.run_tag
        if os.path.exists(pjoin(self.me_dir, 'Source', 'banner_header.txt')):
            self.banner.add(pjoin(self.me_dir, 'Cards','delphes_card.dat'))
            if not delphes3:
                self.banner.add(pjoin(self.me_dir, 'Cards','delphes_trigger.dat'))
            self.banner.write(pjoin(self.me_dir, 'Events', self.run_name, '%s_%s_banner.txt' % (self.run_name, tag)))

        cross = self.results[self.run_name].get_current_info()['cross']

        delphes_log = pjoin(self.me_dir, 'Events', self.run_name, "%s_delphes.log" % tag)
        self.cluster.launch_and_wait(prog,
                        argument= [delphes_dir, self.run_name, tag, str(cross), filepath],
                        stdout=delphes_log, stderr=subprocess.STDOUT,
                        cwd=pjoin(self.me_dir,'Events'))

        if not os.path.exists(pjoin(self.me_dir, 'Events',
                                self.run_name, '%s_delphes_events.lhco' % tag)):
            logger.warning('We do not create lhco by default anymore')


        #eradir = self.options['exrootanalysis_path']
        madir = self.options['madanalysis_path']
        td = self.options['td_path']

        if os.path.exists(pjoin(self.me_dir, 'Events',
                                self.run_name, '%s_delphes_events.lhco' % tag)):
            # Creating plots
            self.create_plot('Delphes')

        if os.path.exists(pjoin(self.me_dir, 'Events', self.run_name,  '%s_delphes_events.lhco' % tag)):
            misc.gzip(pjoin(self.me_dir, 'Events', self.run_name, '%s_delphes_events.lhco' % tag))

        self.update_status('delphes done', level='delphes', makehtml=False)

    ############################################################################
    def get_pid_final_initial_states(self):
        """Find the pid of all particles in the final and initial states"""
        pids = set()
        subproc = [l.strip() for l in open(pjoin(self.me_dir,'SubProcesses',
                                                                 'subproc.mg'))]
        nb_init = self.ninitial
        pat = re.compile(r'''DATA \(IDUP\(I,\d+\),I=1,\d+\)/([\+\-\d,\s]*)/''', re.I)
        for Pdir in subproc:
            text = open(pjoin(self.me_dir, 'SubProcesses', Pdir, 'born_leshouche.inc')).read()
            group = pat.findall(text)
            for particles in group:
                particles = particles.split(',')
                pids.update(set(particles))

        return pids

    ############################################################################
    def get_pdf_input_filename(self):
        """return the name of the file which is used by the pdfset"""

        if self.options["cluster_local_path"] and \
               os.path.exists(self.options["cluster_local_path"]) and \
               self.options['run_mode'] ==1:
            # no need to transfer the pdf.
            return ''
        
        def check_cluster(path):
            if not self.options["cluster_local_path"] or \
                        os.path.exists(self.options["cluster_local_path"]) or\
                        self.options['run_mode'] !=1:
                return path
            main = self.options["cluster_local_path"]
            if os.path.isfile(path):
                filename = os.path.basename(path)
            possible_path = [pjoin(main, filename),
                             pjoin(main, "lhadpf", filename),
                             pjoin(main, "Pdfdata", filename)]
            if any(os.path.exists(p) for p in possible_path):
                return " "
            else:
                return path
                             

        if hasattr(self, 'pdffile') and self.pdffile:
            return self.pdffile
        else:
            for line in open(pjoin(self.me_dir,'Source','PDF','pdf_list.txt')):
                data = line.split()
                if len(data) < 4:
                    continue
                if data[1].lower() == self.run_card['pdlabel'].lower():
                    self.pdffile = check_cluster(pjoin(self.me_dir, 'lib', 'Pdfdata', data[2]))
                    return self.pdffile
            else:
                # possible when using lhapdf
                path = pjoin(self.me_dir, 'lib', 'PDFsets')
                if os.path.exists(path):
                    self.pdffile = path
                else:
                    self.pdffile = " "
                return self.pdffile
                      
    ############################################################################
    def do_open(self, line):
        """Open a text file/ eps file / html file"""

        args = self.split_arg(line)
        # Check Argument validity and modify argument to be the real path
        self.check_open(args)
        file_path = args[0]

        misc.open_file(file_path)

    ############################################################################
    def do_set(self, line, log=True):
        """Set an option, which will be default for coming generations/outputs
        """
        # cmd calls automaticaly post_set after this command.


        args = self.split_arg(line)
        # Check the validity of the arguments
        self.check_set(args)
        # Check if we need to save this in the option file
        if args[0] in self.options_configuration and '--no_save' not in args:
            self.do_save('options --auto')

        if args[0] == "stdout_level":
            if args[1].isdigit():
                logging.root.setLevel(int(args[1]))
                logging.getLogger('madgraph').setLevel(int(args[1]))
            else:
                logging.root.setLevel(eval('logging.' + args[1]))
                logging.getLogger('madgraph').setLevel(eval('logging.' + args[1]))
            if log: logger.info('set output information to level: %s' % args[1])
        elif args[0] == "fortran_compiler":
            if args[1] == 'None':
                args[1] = None
            self.options['fortran_compiler'] = args[1]
            current = misc.detect_current_compiler(pjoin(self.me_dir,'Source','make_opts'), 'fortran')
            if current != args[1] and args[1] != None:
                misc.mod_compilator(self.me_dir, args[1], current, 'gfortran')
        elif args[0] == "cpp_compiler":
            if args[1] == 'None':
                args[1] = None
            self.options['cpp_compiler'] = args[1]
            current = misc.detect_current_compiler(pjoin(self.me_dir,'Source','make_opts'), 'cpp')
            if current != args[1] and args[1] != None:
                misc.mod_compilator(self.me_dir, args[1], current, 'cpp')
        elif args[0] == "run_mode":
            if not args[1] in [0,1,2,'0','1','2']:
                raise self.InvalidCmd, 'run_mode should be 0, 1 or 2.'
            self.cluster_mode = int(args[1])
            self.options['run_mode'] =  self.cluster_mode
        elif args[0] in  ['cluster_type', 'cluster_queue', 'cluster_temp_path']:
            if args[1] == 'None':
                args[1] = None
            self.options[args[0]] = args[1]
            # cluster (re)-initialization done later
            # self.cluster update at the end of the routine
        elif args[0] in ['cluster_nb_retry', 'cluster_retry_wait', 'cluster_size']:
            self.options[args[0]] = int(args[1])
            # self.cluster update at the end of the routine
        elif args[0] == 'nb_core':
            if args[1] == 'None':
                import multiprocessing
                self.nb_core = multiprocessing.cpu_count()
                self.options['nb_core'] = self.nb_core
                return
            if not args[1].isdigit():
                raise self.InvalidCmd('nb_core should be a positive number')
            self.nb_core = int(args[1])
            self.options['nb_core'] = self.nb_core
        elif args[0] == 'timeout':
            self.options[args[0]] = int(args[1])
        elif args[0] == 'cluster_status_update':
            if '(' in args[1]:
                data = ' '.join([a for a in args[1:] if not a.startswith('-')])
                data = data.replace('(','').replace(')','').replace(',',' ').split()
                first, second = data[:2]
            else:
                first, second = args[1:3]

            self.options[args[0]] = (int(first), int(second))
        elif args[0] == 'notification_center':
            if args[1] in ['None','True','False']:
                self.allow_notification_center = eval(args[1])
                self.options[args[0]] = eval(args[1])
            else:
                raise self.InvalidCmd('Not a valid value for notification_center')
        elif args[0] in self.options:
            if args[1] in ['None','True','False']:
                self.options[args[0]] = ast.literal_eval(args[1])
            elif args[0].endswith('path'):
                if os.path.exists(args[1]):
                    self.options[args[0]] = args[1]
                elif os.path.exists(pjoin(self.me_dir, args[1])):
                    self.options[args[0]] = pjoin(self.me_dir, args[1])
                else:
                    raise self.InvalidCmd('Not a valid path: keep previous value: \'%s\'' % self.options[args[0]])
            else:
                self.options[args[0]] = args[1]

    def post_set(self, stop, line):
        """Check if we need to save this in the option file"""
        try:
            args = self.split_arg(line)
            if 'cluster' in args[0] or args[0] == 'run_mode':
                self.configure_run_mode(self.options['run_mode'])


            # Check the validity of the arguments
            self.check_set(args)

            if args[0] in self.options_configuration and '--no_save' not in args:
                self.exec_cmd('save options --auto')
            elif args[0] in self.options_madevent:
                logger.info('This option will be the default in any output that you are going to create in this session.')
                logger.info('In order to keep this changes permanent please run \'save options\'')
            return stop
        except self.InvalidCmd:
            return stop

    def configure_run_mode(self, run_mode):
        """change the way to submit job 0: single core, 1: cluster, 2: multicore"""

        self.cluster_mode = run_mode
        self.options['run_mode'] = run_mode

        if run_mode == 2:
            if not self.options['nb_core']:
                import multiprocessing
                self.options['nb_core'] = multiprocessing.cpu_count()
            nb_core = self.options['nb_core']
        elif run_mode == 0:
            nb_core = 1



        if run_mode in [0, 2]:
            self.cluster = cluster.MultiCore(
                             **self.options)
            self.cluster.nb_core = nb_core
        #cluster_temp_path=self.options['cluster_temp_path'],

        if self.cluster_mode == 1:
            opt = self.options
            cluster_name = opt['cluster_type']
            if cluster_name in cluster.from_name:
                self.cluster = cluster.from_name[cluster_name](**opt)
            else:
                if MADEVENT and ('mg5_path' not in self.options or not self.options['mg5_path']):
                    raise self.InvalidCmd('%s not native cluster type and not MG5aMC found to check for plugin')
                elif MADEVENT:
                    mg5dir = self.options['mg5_path']
                    if mg5dir not in sys.path:
                        sys.path.append(mg5dir)
                else:
                    mg5dir = MG5DIR
                # Check if a plugin define this type of cluster
                # check for PLUGIN format
                for plug in os.listdir(pjoin(mg5dir, 'PLUGIN')):
                    if os.path.exists(pjoin(mg5dir, 'PLUGIN', plug, '__init__.py')):
                        try:
                            __import__('PLUGIN.%s' % plug)
                        except Exception, error:
                            logger.critical('plugin directory %s fail to be loaded. Please check it', plug)
                            continue
                        plugin = sys.modules['PLUGIN.%s' % plug]  
                        if not hasattr(plugin, 'new_cluster'):
                            continue
                        if not misc.is_plugin_supported(plugin):
                            continue              
                        if cluster_name in plugin.new_cluster:
                            logger.info("cluster handling will be done with PLUGIN: %s" % plug,'$MG:color:BLACK')
                            self.cluster = plugin.new_cluster[cluster_name](**opt)
                            break
                else:
                    raise self.InvalidCmd, "%s is not recognized as a supported cluster format." % cluster_name
                
                
    def check_param_card(self, path, run=True, dependent=False):
        """
        1) Check that no scan parameter are present
        2) Check that all the width are define in the param_card.
        - If a scan parameter is define. create the iterator and recall this fonction 
          on the first element.
        - If some width are set on 'Auto', call the computation tools.
        3) if dependent is on True check for dependent parameter (automatic for scan)"""
        
        pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M)  
        pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
        text = open(path).read()
               
        if pattern_scan.search(text):
            if not isinstance(self, cmd.CmdShell):
                # we are in web mode => forbid scan due to security risk
                raise Exception, "Scan are not allowed in web mode"
            # at least one scan parameter found. create an iterator to go trough the cards
            main_card = check_param_card.ParamCardIterator(text)
            self.param_card_iterator = main_card
            first_card = main_card.next(autostart=True)
            first_card.write(path)
            return self.check_param_card(path, run, dependent=True)
        
        pdg_info = pattern_width.findall(text)
        if pdg_info:
            if run:
                logger.info('Computing the width set on auto in the param_card.dat')
                has_nlo = any(nlo.lower()=="@nlo" for _,nlo in pdg_info)
                pdg = [pdg for pdg,nlo in pdg_info]
                if not has_nlo:
                    self.do_compute_widths('%s %s' % (' '.join(pdg), path))
                else:
                    self.do_compute_widths('%s %s --nlo' % (' '.join(pdg), path)) 
            else:
                logger.info('''Some width are on Auto in the card. 
    Those will be computed as soon as you have finish the edition of the cards.
    If you want to force the computation right now and being able to re-edit
    the cards afterwards, you can type \"compute_wdiths\".''')
        
        if dependent:   
            card = check_param_card.ParamCard(path)
            AskforEditCard.update_dependent(self, self.me_dir, card, path, timer=20)
        
        return

    def add_error_log_in_html(self, errortype=None):
        """If a ME run is currently running add a link in the html output"""

        # Be very carefull to not raise any error here (the traceback
        #will be modify in that case.)
        if hasattr(self, 'results') and hasattr(self.results, 'current') and\
                self.results.current and 'run_name' in self.results.current and \
                hasattr(self, 'me_dir'):
            name = self.results.current['run_name']
            tag = self.results.current['tag']
            self.debug_output = pjoin(self.me_dir, '%s_%s_debug.log' % (name,tag))
            if errortype:
                self.results.current.debug = errortype
            else:
                self.results.current.debug = self.debug_output

        else:
            #Force class default
            self.debug_output = CommonRunCmd.debug_output
        if os.path.exists('ME5_debug') and not 'ME5_debug' in self.debug_output:
            os.remove('ME5_debug')
        if not 'ME5_debug' in self.debug_output:
            os.system('ln -s %s ME5_debug &> /dev/null' % self.debug_output)


    def do_quit(self, line):
        """Not in help: exit """

        if not self.force_run:
            try:
                os.remove(pjoin(self.me_dir,'RunWeb'))
            except Exception:
                pass
        try:
            self.store_result()
        except Exception:
            # If nothing runs they they are no result to update
            pass

        try:
            self.update_status('', level=None)
        except Exception, error:
            pass
        devnull = open(os.devnull, 'w')
        try:
            misc.call(['./bin/internal/gen_cardhtml-pl'], cwd=self.me_dir,
                        stdout=devnull, stderr=devnull)
        except Exception:
            pass
        devnull.close()

        return super(CommonRunCmd, self).do_quit(line)

    # Aliases
    do_EOF = do_quit
    do_exit = do_quit


    def update_status(self, status, level, makehtml=True, force=True,
                      error=False, starttime = None, update_results=True,
                      print_log=True):
        """ update the index status """

        if makehtml and not force:
            if hasattr(self, 'next_update') and time.time() < self.next_update:
                return
            else:
                self.next_update = time.time() + 3

        if print_log:
            if isinstance(status, str):
                if '<br>' not  in status:
                    logger.info(status)
            elif starttime:
                running_time = misc.format_timer(time.time()-starttime)
                logger.info(' Idle: %s,  Running: %s,  Completed: %s [ %s ]' % \
                           (status[0], status[1], status[2], running_time))
            else:
                logger.info(' Idle: %s,  Running: %s,  Completed: %s' % status[:3])

        if 'arXiv' in status:
            if '[' in status:
                status = status.split('[',1)[0]
            else:
                status = status.split('arXiv',1)[0]

        if update_results:
            self.results.update(status, level, makehtml=makehtml, error=error)

    ############################################################################
    def keep_cards(self, need_card=[], ignore=[]):
        """Ask the question when launching generate_events/multi_run"""

        check_card = ['pythia_card.dat', 'pgs_card.dat','delphes_card.dat',
                      'delphes_trigger.dat', 'madspin_card.dat', 'shower_card.dat',
                      'reweight_card.dat','pythia8_card.dat',
                      'madanalysis5_parton_card.dat','madanalysis5_hadron_card.dat',]

        cards_path = pjoin(self.me_dir,'Cards')
        for card in check_card:
            if card in ignore or (ignore == ['*'] and card not in need_card):
                continue
            if card not in need_card:
                if os.path.exists(pjoin(cards_path, card)):
                    files.mv(pjoin(cards_path, card), pjoin(cards_path, '.%s' % card))
            else:
                if not os.path.exists(pjoin(cards_path, card)):
                    if os.path.exists(pjoin(cards_path, '.%s' % card)):
                        files.mv(pjoin(cards_path, '.%s' % card), pjoin(cards_path, card))
                    else:
                        default = card.replace('.dat', '_default.dat')
                        files.cp(pjoin(cards_path, default),pjoin(cards_path, card))

    ############################################################################
    def set_configuration(self, config_path=None, final=True, initdir=None, amcatnlo=False):
        """ assign all configuration variable from file
            ./Cards/mg5_configuration.txt. assign to default if not define """

        if not hasattr(self, 'options') or not self.options:
            self.options = dict(self.options_configuration)
            self.options.update(self.options_madgraph)
            self.options.update(self.options_madevent)

        if not config_path:
            if os.environ.has_key('MADGRAPH_BASE'):
                config_path = pjoin(os.environ['MADGRAPH_BASE'],'mg5_configuration.txt')
                self.set_configuration(config_path=config_path, final=False)
            if 'HOME' in os.environ:
                config_path = pjoin(os.environ['HOME'],'.mg5',
                                                        'mg5_configuration.txt')
                if os.path.exists(config_path):
                    self.set_configuration(config_path=config_path,  final=False)
            if amcatnlo:
                me5_config = pjoin(self.me_dir, 'Cards', 'amcatnlo_configuration.txt')
            else:
                me5_config = pjoin(self.me_dir, 'Cards', 'me5_configuration.txt')
            self.set_configuration(config_path=me5_config, final=False, initdir=self.me_dir)

            if self.options.has_key('mg5_path') and self.options['mg5_path']:
                MG5DIR = self.options['mg5_path']
                config_file = pjoin(MG5DIR, 'input', 'mg5_configuration.txt')
                self.set_configuration(config_path=config_file, final=False,initdir=MG5DIR)
            else:
                self.options['mg5_path'] = None
            return self.set_configuration(config_path=me5_config, final=final,initdir=self.me_dir)

        config_file = open(config_path)

        # read the file and extract information
        logger.info('load configuration from %s ' % config_file.name)
        for line in config_file:
            
            if '#' in line:
                line = line.split('#',1)[0]
            line = line.replace('\n','').replace('\r\n','')
            try:
                name, value = line.split('=')
            except ValueError:
                pass
            else:
                name = name.strip()
                value = value.strip()
                if name.endswith('_path') and not name.startswith('cluster'):
                    path = value
                    if os.path.isdir(path):
                        self.options[name] = os.path.realpath(path)
                        continue
                    if not initdir:
                        continue
                    path = pjoin(initdir, value)
                    if os.path.isdir(path):
                        self.options[name] = os.path.realpath(path)
                        continue
                else:
                    self.options[name] = value
                    if value.lower() == "none":
                        self.options[name] = None

        if not final:
            return self.options # the return is usefull for unittest

        # Treat each expected input
        # delphes/pythia/... path
        for key in self.options:
            # Final cross check for the path
            if key.endswith('path') and not key.startswith("cluster"):
                path = self.options[key]
                if path is None:
                    continue
                if os.path.isdir(path):
                    self.options[key] = os.path.realpath(path)
                    continue
                path = pjoin(self.me_dir, self.options[key])
                if os.path.isdir(path):
                    self.options[key] = os.path.realpath(path)
                    continue
                elif self.options.has_key('mg5_path') and self.options['mg5_path']:
                    path = pjoin(self.options['mg5_path'], self.options[key])
                    if os.path.isdir(path):
                        self.options[key] = os.path.realpath(path)
                        continue
                self.options[key] = None
            elif key.startswith('cluster') and key != 'cluster_status_update':
                if key in ('cluster_nb_retry','cluster_wait_retry'):
                    self.options[key] = int(self.options[key]) 
                if hasattr(self,'cluster'):
                    del self.cluster
                pass
            elif key == 'automatic_html_opening':
                if self.options[key] in ['False', 'True']:
                    self.options[key] =ast.literal_eval(self.options[key])
            elif key == "notification_center":
                if self.options[key] in ['False', 'True']:
                    self.allow_notification_center =ast.literal_eval(self.options[key])
                    self.options[key] =ast.literal_eval(self.options[key])
            elif key not in ['text_editor','eps_viewer','web_browser','stdout_level',
                              'complex_mass_scheme', 'gauge', 'group_subprocesses']:
                # Default: try to set parameter
                try:
                    self.do_set("%s %s --no_save" % (key, self.options[key]), log=False)
                except self.InvalidCmd:
                    logger.warning("Option %s from config file not understood" \
                                   % key)

        # Configure the way to open a file:
        misc.open_file.configure(self.options)
        self.configure_run_mode(self.options['run_mode'])
        return self.options

    @staticmethod
    def find_available_run_name(me_dir):
        """ find a valid run_name for the current job """

        name = 'run_%02d'
        data = [int(s[4:j]) for s in os.listdir(pjoin(me_dir,'Events')) for 
                j in range(4,len(s)+1) if \
                s.startswith('run_') and s[4:j].isdigit()]
        return name % (max(data+[0])+1)


    ############################################################################
    def do_decay_events(self,line):
        """Require MG5 directory: decay events with spin correlations
        """

        if '-from_cards' in line and not os.path.exists(pjoin(self.me_dir, 'Cards', 'madspin_card.dat')):
            return

        # First need to load MadSpin

        # Check that MG5 directory is present .
        if MADEVENT and not self.options['mg5_path']:
            raise self.InvalidCmd, '''The module decay_events requires that MG5 is installed on the system.
            You can install it and set its path in ./Cards/me5_configuration.txt'''
        elif MADEVENT:
            sys.path.append(self.options['mg5_path'])
        try:
            import MadSpin.decay as decay
            import MadSpin.interface_madspin as interface_madspin
        except ImportError:
            if __debug__:
                raise
            else:
                raise self.ConfigurationError, '''Can\'t load MadSpin
            The variable mg5_path might not be correctly configured.'''

        self.update_status('Running MadSpin', level='madspin')
        if not '-from_cards' in line:
            self.keep_cards(['madspin_card.dat'], ignore=['*'])
            self.ask_edit_cards(['madspin_card.dat'], 'fixed', plot=False)
        self.help_decay_events(skip_syntax=True)

        # load the name of the event file
        args = self.split_arg(line)
        self.check_decay_events(args)
        # args now alway content the path to the valid files
        madspin_cmd = interface_madspin.MadSpinInterface(args[0])
        madspin_cmd.update_status = lambda *x,**opt: self.update_status(*x, level='madspin',**opt)

        path = pjoin(self.me_dir, 'Cards', 'madspin_card.dat')

        madspin_cmd.import_command_file(path)

        # create a new run_name directory for this output
        i = 1
        while os.path.exists(pjoin(self.me_dir,'Events', '%s_decayed_%i' % (self.run_name,i))):
            i+=1
        new_run = '%s_decayed_%i' % (self.run_name,i)
        evt_dir = pjoin(self.me_dir, 'Events')

        os.mkdir(pjoin(evt_dir, new_run))
        current_file = args[0].replace('.lhe', '_decayed.lhe')
        new_file = pjoin(evt_dir, new_run, os.path.basename(args[0]))
        if not os.path.exists(current_file):
            if os.path.exists(current_file+'.gz'):
                current_file += '.gz'
                new_file += '.gz'
            else:
                logger.error('MadSpin fails to create any decayed file.')
                return

        files.mv(current_file, new_file)
        logger.info("The decayed event file has been moved to the following location: ")
        logger.info(new_file)

        if hasattr(self, 'results'):
            current = self.results.current
            nb_event = self.results.current['nb_event']
            if not nb_event:
                current = self.results[self.run_name][0]
                nb_event = current['nb_event']

            cross = current['cross']
            error = current['error']
            self.results.add_run( new_run, self.run_card)
            self.results.add_detail('nb_event', int(nb_event*madspin_cmd.efficiency))
            self.results.add_detail('cross', madspin_cmd.cross)#cross * madspin_cmd.branching_ratio)
            self.results.add_detail('error', madspin_cmd.error+ cross * madspin_cmd.err_branching_ratio)
            self.results.add_detail('run_mode', current['run_mode'])

        self.run_name = new_run
        self.banner = madspin_cmd.banner
        self.banner.add(path)
        self.banner.write(pjoin(self.me_dir,'Events',self.run_name, '%s_%s_banner.txt' %
                                (self.run_name, self.run_tag)))
        self.update_status('MadSpin Done', level='parton', makehtml=False)
        if 'unweighted' in os.path.basename(args[0]):
            self.create_plot('parton')

    def complete_decay_events(self, text, line, begidx, endidx):
        args = self.split_arg(line[0:begidx], error=False)
        if len(args) == 1:
            return self.complete_plot(text, line, begidx, endidx)
        else:
            return

    def complete_print_results(self,text, line, begidx, endidx):
        "Complete the print results command"
        args = self.split_arg(line[0:begidx], error=False) 
        if len(args) == 1:
            #return valid run_name
            data = misc.glob(pjoin('*','unweighted_events.lhe.gz'),
                             pjoin(self.me_dir, 'Events')) 

            data = [n.rsplit('/',2)[1] for n in data]
            tmp1 =  self.list_completion(text, data)
            return tmp1        
        else:
            data = misc.glob('*_pythia_events.hep.gz', pjoin(self.me_dir, 'Events', args[0]))
            data = [os.path.basename(p).rsplit('_',1)[0] for p in data]
            data += ["--mode=a", "--mode=w", "--path=", "--format=short"]
            tmp1 =  self.list_completion(text, data)
            return tmp1
            
    def help_print_result(self):
        logger.info("syntax: print_result [RUN] [TAG] [options]")
        logger.info("-- show in text format the status of the run (cross-section/nb-event/...)")
        logger.info("--path= defines the path of the output file.")
        logger.info("--mode=a allow to add the information at the end of the file.")
        logger.info("--format=short (only if --path is define)")
        logger.info("        allows to have a multi-column output easy to parse")

    ############################################################################
    def do_check_events(self, line):
        """ Run some sanity check on the generated events."""
                
        # Check that MG5 directory is present .
        if MADEVENT and not self.options['mg5_path']:
            raise self.InvalidCmd, '''The module reweight requires that MG5 is installed on the system.
            You can install it and set its path in ./Cards/me5_configuration.txt'''
        elif MADEVENT:
            sys.path.append(self.options['mg5_path'])
        try:
            import madgraph.interface.reweight_interface as reweight_interface
        except ImportError:
            raise self.ConfigurationError, '''Can\'t load Reweight module.
            The variable mg5_path might not be correctly configured.'''
              

        # load the name of the event file
        args = self.split_arg(line) 
        self.check_check_events(args) 
        # args now alway content the path to the valid files
        reweight_cmd = reweight_interface.ReweightInterface(args[0], allow_madspin=True)
        reweight_cmd.mother = self
        self.update_status('Running check on events', level='check')
        
        reweight_cmd.check_events()
        
    ############################################################################
    def complete_check_events(self, text, line, begidx, endidx):
        args = self.split_arg(line[0:begidx], error=False)

        if len(args) == 1 and os.path.sep not in text:
            #return valid run_name
            data = misc.glob(pjoin('*','*events.lhe*'), pjoin(self.me_dir, 'Events'))
            data = [n.rsplit('/',2)[1] for n in data]
            return  self.list_completion(text, data, line)
        else:
            return self.path_completion(text,
                                        os.path.join('.',*[a for a in args \
                                                    if a.endswith(os.path.sep)]))

    def complete_reweight(self,text, line, begidx, endidx):
        "Complete the pythia command"
        args = self.split_arg(line[0:begidx], error=False)

        #return valid run_name
        data = misc.glob(pjoin('*','*events.lhe*'), pjoin(self.me_dir, 'Events'))
        data = [n.rsplit('/',2)[1] for n in data]
        if not '-f' in args:
            data.append('-f')
        tmp1 =  self.list_completion(text, data)
        return tmp1 



    def complete_compute_widths(self, text, line, begidx, endidx, formatting=True):
        "Complete the compute_widths command"

        args = self.split_arg(line[0:begidx])
        
        if args[-1] in  ['--path=', '--output=']:
            completion = {'path': self.path_completion(text)}
        elif line[begidx-1] == os.path.sep:
            current_dir = pjoin(*[a for a in args if a.endswith(os.path.sep)])
            if current_dir.startswith('--path='):
                current_dir = current_dir[7:]
            if current_dir.startswith('--output='):
                current_dir = current_dir[9:]                
            completion = {'path': self.path_completion(text, current_dir)}
        else:
            completion = {}            
            completion['options'] = self.list_completion(text, 
                            ['--path=', '--output=', '--min_br=0.\$', '--nlo',
                             '--precision_channel=0.\$', '--body_decay='])            
        
        return self.deal_multiple_categories(completion, formatting)
        

    def update_make_opts(self):
        """update the make_opts file writing the environmental variables
        stored in make_opts_var"""
        make_opts = os.path.join(self.me_dir, 'Source', 'make_opts')
        
        # Set some environment variables common to all interfaces
        if not hasattr(self,'options') or not 'pythia8_path' in self.options or \
           not self.options['pythia8_path'] or \
           not os.path.isfile(pjoin(self.options['pythia8_path'],'bin','pythia8-config')):
            self.make_opts_var['PYTHIA8_PATH']='NotInstalled'
        else:
            self.make_opts_var['PYTHIA8_PATH']=self.options['pythia8_path']

        self.make_opts_var['MG5AMC_VERSION'] = misc.get_pkg_info()['version']

        return self.update_make_opts_full(make_opts, self.make_opts_var)

    @staticmethod
    def update_make_opts_full(path, def_variables, keep_old=True):
        """update the make_opts file writing the environmental variables
        of def_variables.
        if a value of the dictionary is None then it is not written.
        """
        make_opts = path
        pattern = re.compile(r'^(\w+)\s*=\s*(.*)$',re.DOTALL)
        diff = False # set to True if one varible need to be updated 
                     #if on False the file is not modify
        
        tag = '#end_of_make_opts_variables\n'
        make_opts_variable = True # flag to say if we are in edition area or not
        content = []
        variables = dict(def_variables)
        need_keys = variables.keys()
        for line in open(make_opts):
            line = line.strip()
            if make_opts_variable: 
                if line.startswith('#') or not line:
                    if line.startswith('#end_of_make_opts_variables'):
                        make_opts_variable = False
                    continue
                elif pattern.search(line):
                    key, value = pattern.search(line).groups()
                    if key not in variables:
                        variables[key] = value
                    elif value !=  variables[key]:
                        diff=True
                    else:
                        need_keys.remove(key)
                else: 
                    misc.sprint("end on line", line)
                    make_opts_variable = False
                    content.append(line)
            else:                  
                content.append(line)
                     
        if need_keys:
            diff=True #This means that new definition are added to the file. 

        content_variables = '\n'.join('%s=%s' % (k,v) for k, v in variables.items() if v is not None)
        content_variables += '\n%s' % tag

        if diff:
            with open(make_opts, 'w') as fsock: 
                fsock.write(content_variables + '\n'.join(content))
        return       


# lhapdf-related functions
    def link_lhapdf(self, libdir, extra_dirs = []):
        """links lhapdf into libdir"""

        lhapdf_version = self.get_lhapdf_version()
        logger.info('Using LHAPDF v%s interface for PDFs' % lhapdf_version)
        lhalibdir = subprocess.Popen([self.options['lhapdf'], '--libdir'],
                 stdout = subprocess.PIPE).stdout.read().strip()

        if lhapdf_version.startswith('5.'):
            pdfsetsdir = subprocess.Popen([self.options['lhapdf'], '--pdfsets-path'],
                 stdout = subprocess.PIPE).stdout.read().strip()
        else:
            pdfsetsdir = subprocess.Popen([self.options['lhapdf'], '--datadir'],
                 stdout = subprocess.PIPE).stdout.read().strip()
        
        self.lhapdf_pdfsets = self.get_lhapdf_pdfsets_list(pdfsetsdir)
        # link the static library in lib
        lhalib = 'libLHAPDF.a'

        if os.path.exists(pjoin(libdir, lhalib)):
            files.rm(pjoin(libdir, lhalib))
        files.ln(pjoin(lhalibdir, lhalib), libdir)
        # just create the PDFsets dir, the needed PDF set will be copied at run time
        if not os.path.isdir(pjoin(libdir, 'PDFsets')):
            os.mkdir(pjoin(libdir, 'PDFsets'))
        self.make_opts_var['lhapdf'] = self.options['lhapdf']
        self.make_opts_var['lhapdfversion'] = lhapdf_version[0]
        self.make_opts_var['lhapdf_config'] = self.options['lhapdf']


    def get_characteristics(self, path=None):
        """reads the proc_characteristics file and initialises the correspondant
        dictionary"""
        
        if not path:
            path = os.path.join(self.me_dir, 'SubProcesses', 'proc_characteristics')
        
        self.proc_characteristics = banner_mod.ProcCharacteristic(path)
        return self.proc_characteristics


    def copy_lhapdf_set(self, lhaid_list, pdfsets_dir):
        """copy (if needed) the lhapdf set corresponding to the lhaid in lhaid_list 
        into lib/PDFsets"""

        if not hasattr(self, 'lhapdf_pdfsets'):
            self.lhapdf_pdfsets = self.get_lhapdf_pdfsets_list(pdfsets_dir)

        pdfsetname=set()
        for lhaid in lhaid_list:
            if isinstance(lhaid, (int,float)):
                try:
                    if lhaid in self.lhapdf_pdfsets:
                        pdfsetname.add(self.lhapdf_pdfsets[lhaid]['filename'])
                    else:
                        raise MadGraph5Error('lhaid %s not valid input number for the current lhapdf' % lhaid )
                except KeyError:
                    if self.lhapdf_version.startswith('5'):
                        raise MadGraph5Error(\
                            ('invalid lhaid set in th run_card: %d .\nPlease note that some sets' % lhaid) + \
                             '(eg MSTW 90%CL error sets) \nare not available in aMC@NLO + LHAPDF 5.x.x')
                    else:
                        logger.debug('%d not found in pdfsets.index' % lhaid)
            else:
                pdfsetname.add(lhaid)

        # check if the file exists, otherwise install it:
        # also check that the PDFsets dir exists, otherwise create it.
        # if fails, install the lhapdfset into lib/PDFsets
        if not os.path.isdir(pdfsets_dir):
            try:
                os.mkdir(pdfsets_dir)
            except OSError:
                pdfsets_dir = pjoin(self.me_dir, 'lib', 'PDFsets')
        elif os.path.exists(pjoin(self.me_dir, 'lib', 'PDFsets')):
            #clean previous set of pdf used
            for name in os.listdir(pjoin(self.me_dir, 'lib', 'PDFsets')):
                if name not in pdfsetname:
                    try:
                        if os.path.isdir(pjoin(self.me_dir, 'lib', 'PDFsets', name)):
                            shutil.rmtree(pjoin(self.me_dir, 'lib', 'PDFsets', name))
                        else:
                            os.remove(pjoin(self.me_dir, 'lib', 'PDFsets', name))
                    except Exception, error:
                        logger.debug('%s', error)
        
        if self.options["cluster_local_path"]:
            lhapdf_cluster_possibilities = [self.options["cluster_local_path"],
                                      pjoin(self.options["cluster_local_path"], "lhapdf"),
                                      pjoin(self.options["cluster_local_path"], "lhapdf", "pdfsets"),
                                      pjoin(self.options["cluster_local_path"], "..", "lhapdf"),
                                      pjoin(self.options["cluster_local_path"], "..", "lhapdf", "pdfsets"),
                                      pjoin(self.options["cluster_local_path"], "..", "lhapdf","pdfsets", "6.1")
                                      ]
        else:
            lhapdf_cluster_possibilities = []

        for pdfset in pdfsetname:
        # Check if we need to copy the pdf
            if self.options["cluster_local_path"] and self.options["run_mode"] == 1 and \
                any((os.path.exists(pjoin(d, pdfset)) for d in lhapdf_cluster_possibilities)):
    
                os.environ["LHAPATH"] = [d for d in lhapdf_cluster_possibilities if os.path.exists(pjoin(d, pdfset))][0]
                os.environ["CLUSTER_LHAPATH"] = os.environ["LHAPATH"]
                # no need to copy it
                if os.path.exists(pjoin(pdfsets_dir, pdfset)):
                    try:
                        if os.path.isdir(pjoin(pdfsets_dir, name)):
                            shutil.rmtree(pjoin(pdfsets_dir, name))
                        else:
                            os.remove(pjoin(pdfsets_dir, name))
                    except Exception, error:
                        logger.debug('%s', error)
        
            #check that the pdfset is not already there
            elif not os.path.exists(pjoin(self.me_dir, 'lib', 'PDFsets', pdfset)) and \
               not os.path.isdir(pjoin(self.me_dir, 'lib', 'PDFsets', pdfset)):
    
                if pdfset and not os.path.exists(pjoin(pdfsets_dir, pdfset)):
                    self.install_lhapdf_pdfset(pdfsets_dir, pdfset)
    
                if os.path.exists(pjoin(pdfsets_dir, pdfset)):
                    files.cp(pjoin(pdfsets_dir, pdfset), pjoin(self.me_dir, 'lib', 'PDFsets'))
                elif os.path.exists(pjoin(os.path.dirname(pdfsets_dir), pdfset)):
                    files.cp(pjoin(os.path.dirname(pdfsets_dir), pdfset), pjoin(self.me_dir, 'lib', 'PDFsets'))
            
    def install_lhapdf_pdfset(self, pdfsets_dir, filename):
        """idownloads and install the pdfset filename in the pdfsets_dir"""
        lhapdf_version = self.get_lhapdf_version()
        local_path = pjoin(self.me_dir, 'lib', 'PDFsets')
        return self.install_lhapdf_pdfset_static(self.options['lhapdf'],
                                             pdfsets_dir, filename,
                                             lhapdf_version=lhapdf_version,
                                             alternate_path=local_path)
                                                 

    @staticmethod
    def install_lhapdf_pdfset_static(lhapdf_config, pdfsets_dir, filename, 
                                        lhapdf_version=None, alternate_path=None):
        """idownloads and install the pdfset filename in the pdfsets_dir.
        Version which can be used independently of the class.
        local path is used if the global installation fails.
        """
           
        if not lhapdf_version:
            lhapdf_version = subprocess.Popen([lhapdf_config, '--version'], 
                        stdout = subprocess.PIPE).stdout.read().strip()
        if not pdfsets_dir:
            pdfsets_dir = subprocess.Popen([lhapdf_config, '--datadir'], 
                        stdout = subprocess.PIPE).stdout.read().strip()
                                
        if isinstance(filename, int):
            pdf_info = CommonRunCmd.get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version)
            filename = pdf_info[filename]['filename']
        
        if os.path.exists(pjoin(pdfsets_dir, filename)):
            logger.debug('%s is already present in %s', filename, pdfsets_dir)
            return
             
        logger.info('Trying to download %s' % filename)

        if lhapdf_version.startswith('5.'):

            # use the lhapdf-getdata command, which is in the same path as
            # lhapdf-config
            getdata = lhapdf_config.replace('lhapdf-config', ('lhapdf-getdata'))
            misc.call([getdata, filename], cwd = pdfsets_dir)

        elif lhapdf_version.startswith('6.'):
            # use the "lhapdf install xxx" command, which is in the same path as
            # lhapdf-config
            getdata = lhapdf_config.replace('lhapdf-config', ('lhapdf'))

            misc.call([getdata, 'install', filename], cwd = pdfsets_dir)

        else:
            raise MadGraph5Error('Not valid LHAPDF version: %s' % lhapdf_version)
        
        # check taht the file has been installed in the global dir
        if os.path.exists(pjoin(pdfsets_dir, filename)) or \
           os.path.isdir(pjoin(pdfsets_dir, filename)):
            logger.info('%s successfully downloaded and stored in %s' \
                    % (filename, pdfsets_dir))
        #otherwise (if v5) save it locally
        elif lhapdf_version.startswith('5.'):
            logger.warning('Could not download %s into %s. Trying to save it locally' \
                    % (filename, pdfsets_dir))
            CommonRunCmd.install_lhapdf_pdfset_static(lhapdf_config, alternate_path, filename,
                                                      lhapdf_version=lhapdf_version)
        elif lhapdf_version.startswith('6.') and '.LHgrid' in filename:
            logger.info('Could not download %s: Try %s', filename, filename.replace('.LHgrid',''))
            return CommonRunCmd.install_lhapdf_pdfset_static(lhapdf_config, pdfsets_dir, 
                                                              filename.replace('.LHgrid',''), 
                                        lhapdf_version, alternate_path)
            
        else:
            raise MadGraph5Error, \
                'Could not download %s into %s. Please try to install it manually.' \
                    % (filename, pdfsets_dir)



    def get_lhapdf_pdfsets_list(self, pdfsets_dir):
        """read the PDFsets.index file, which should be located in the same
        place as pdfsets_dir, and return a list of dictionaries with the information
        about each pdf set"""
        lhapdf_version = self.get_lhapdf_version()
        return self.get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version)

    @staticmethod
    def get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version):

        if lhapdf_version.startswith('5.'):
            if os.path.exists('%s.index' % pdfsets_dir):
                indexfile = '%s.index' % pdfsets_dir
            else:
                raise MadGraph5Error, 'index of lhapdf file not found'
            pdfsets_lines = \
                    [l for l in open(indexfile).read().split('\n') if l.strip() and \
                        not '90cl' in l]
            lhapdf_pdfsets = dict( (int(l.split()[0]), {'lhaid': int(l.split()[0]),
                          'pdflib_ntype': int(l.split()[1]),
                          'pdflib_ngroup': int(l.split()[2]),
                          'pdflib_nset': int(l.split()[3]),
                          'filename': l.split()[4],
                          'lhapdf_nmem': int(l.split()[5]),
                          'q2min': float(l.split()[6]),
                          'q2max': float(l.split()[7]),
                          'xmin': float(l.split()[8]),
                          'xmax': float(l.split()[9]),
                          'description': l.split()[10]}) \
                         for l in pdfsets_lines)

        elif lhapdf_version.startswith('6.'):
            pdfsets_lines = \
                    [l for l in open(pjoin(pdfsets_dir, 'pdfsets.index')).read().split('\n') if l.strip()]
            lhapdf_pdfsets = dict( (int(l.split()[0]), 
                        {'lhaid': int(l.split()[0]),
                          'filename': l.split()[1]}) \
                         for l in pdfsets_lines)

        else:
            raise MadGraph5Error('Not valid LHAPDF version: %s' % lhapdf_version)

        return lhapdf_pdfsets


    def get_lhapdf_version(self):
        """returns the lhapdf version number"""
        if not hasattr(self, 'lhapdfversion'):
            try:
                self.lhapdf_version = \
                    subprocess.Popen([self.options['lhapdf'], '--version'], 
                        stdout = subprocess.PIPE).stdout.read().strip()
            except OSError, error:
                if error.errno == 2:
                    raise Exception, 'lhapdf executable (%s) is not found on your system. Please install it and/or indicate the path to the correct executable in input/mg5_configuration.txt' % self.options['lhapdf']
                else:
                    raise
                
        # this will be removed once some issues in lhapdf6 will be fixed
        if self.lhapdf_version.startswith('6.0'):
            raise MadGraph5Error('LHAPDF 6.0.x not supported. Please use v6.1 or later')

        return self.lhapdf_version


    def get_lhapdf_pdfsetsdir(self):
        lhapdf_version = self.get_lhapdf_version()

        # check if the LHAPDF_DATA_PATH variable is defined
        if 'LHAPDF_DATA_PATH' in os.environ.keys() and os.environ['LHAPDF_DATA_PATH']:
            datadir = os.environ['LHAPDF_DATA_PATH']

        elif lhapdf_version.startswith('5.'):
            datadir = subprocess.Popen([self.options['lhapdf'], '--pdfsets-path'],
                         stdout = subprocess.PIPE).stdout.read().strip()

        elif lhapdf_version.startswith('6.'):
            datadir = subprocess.Popen([self.options['lhapdf'], '--datadir'],
                         stdout = subprocess.PIPE).stdout.read().strip()

        return datadir

    def get_lhapdf_libdir(self):
        lhapdf_version = self.get_lhapdf_version()

        if lhapdf_version.startswith('5.'):
            libdir = subprocess.Popen([self.options['lhapdf-config'], '--libdir'],
                         stdout = subprocess.PIPE).stdout.read().strip()

        elif lhapdf_version.startswith('6.'):
            libdir = subprocess.Popen([self.options['lhapdf'], '--libs'],
                         stdout = subprocess.PIPE).stdout.read().strip()

        return libdir

class AskforEditCard(cmd.OneLinePathCompletion):
    """A class for asking a question where in addition you can have the
    set command define and modifying the param_card/run_card correctly"""

    all_card_name = ['param_card', 'run_card', 'pythia_card', 'pythia8_card', 
                     'madweight_card', 'MadLoopParams', 'shower_card']

    special_shortcut = {'ebeam':['run_card ebeam1 %(0)s', 'run_card ebeam2 %(0)s'],
                        'lpp': ['run_card lpp1 %(0)s', 'run_card lpp2 %(0)s' ],
                        'lhc': ['run_card lpp1 1', 'run_card lpp2 1', 'run_card ebeam1 %(0)s*1000/2', 'run_card ebeam2 %(0)s*1000/2'],
                        'lep': ['run_card lpp1 0', 'run_card lpp2 0', 'run_card ebeam1 %(0)s/2', 'run_card ebeam2 %(0)s/2'],
                        'ilc': ['run_card lpp1 0', 'run_card lpp2 0', 'run_card ebeam1 %(0)s/2', 'run_card ebeam2 %(0)s/2'],
                        'lcc':['run_card lpp1 1', 'run_card lpp2 1', 'run_card ebeam1 %(0)s*1000/2', 'run_card ebeam2 %(0)s*1000/2'],
                        'fixed_scale': ['run_card fixed_fac_scale T', 'run_card fixed_ren_scale T', 'run_card scale %(0)s', 'run_card dsqrt_q2fact1 %(0)s' ,'run_card dsqrt_q2fact2 %(0)s'],
                        }

    special_shortcut_help = {              
    'ebeam' : 'syntax: set ebeam VALUE:\n      This parameter sets the energy to both beam to the value in GeV',
    'lpp'   : 'syntax: set ebeam  VALUE:\n'+\
              '   Set the type of beam to a given value for both beam\n'+\
              '   0 : means no PDF\n'+\
              '   1 : means proton PDF\n'+\
              '  -1 : means proton PDF\n'+\
              '   2 : means PDF for elastic photon emited from a proton\n'+\
              '   3 : means PDF for elastic photon emited from an electron',
    'lhc'   : 'syntax: set lhc VALUE:\n      Set for a proton-proton collision with that given center of mass energy (in TeV)',
    'lep'   : 'syntax: set lep VALUE:\n      Set for a electron-positron collision with that given center of mass energy (in GeV)',
    'fixed_scale' : 'syntax: set fixed_scale VALUE:\n      Set all scales to the give value (in GeV)',              
    }
    
    def load_default(self):
        """ define all default variable. No load of card here.
            This allow to subclass this class and just change init and still have
            all variables defined."""
    
        self.me_dir = None
        self.param_card = None
        self.run_card = {}
        self.pname2block = {}
        self.conflict = []
        self.restricted_value = {}
        self.mode = ''
        self.cards = [] 
        self.run_set = []
        self.has_mw = False
        self.has_ml = False   
        self.has_shower = False
        self.has_PY8 = False
        self.paths = {}

    
    def define_paths(self, **opt):
        # Initiation
        if 'pwd' in opt:
            self.me_dir = opt['pwd']
        elif 'mother_interface' in opt:
            self.mother_interface = opt['mother_interface']     
        if not hasattr(self, 'me_dir') or not self.me_dir:
            self.me_dir = self.mother_interface.me_dir
        
        #define paths
        self.paths['param'] = pjoin(self.me_dir,'Cards','param_card.dat')
        self.paths['param_default'] = pjoin(self.me_dir,'Cards','param_card_default.dat')
        self.paths['run'] = pjoin(self.me_dir,'Cards','run_card.dat')
        self.paths['run_default'] = pjoin(self.me_dir,'Cards','run_card_default.dat')
        self.paths['transfer'] =pjoin(self.me_dir,'Cards','transfer_card.dat')
        self.paths['MadWeight'] =pjoin(self.me_dir,'Cards','MadWeight_card.dat')
        self.paths['MadWeight_default'] =pjoin(self.me_dir,'Cards','MadWeight_card_default.dat')
        self.paths['ML'] =pjoin(self.me_dir,'Cards','MadLoopParams.dat')
        self.paths['shower'] = pjoin(self.me_dir,'Cards','shower_card.dat')
        self.paths['shower_default'] = pjoin(self.me_dir,'Cards','shower_card_default.dat')
        self.paths['pythia'] =pjoin(self.me_dir, 'Cards','pythia_card.dat')
        self.paths['PY8'] = pjoin(self.me_dir, 'Cards','pythia8_card.dat')
        self.paths['PY8_default'] = pjoin(self.me_dir, 'Cards','pythia8_card_default.dat')
        self.paths['madspin_default'] = pjoin(self.me_dir,'Cards/madspin_card_default.dat')
        self.paths['madspin'] = pjoin(self.me_dir,'Cards/madspin_card.dat')
        self.paths['reweight'] = pjoin(self.me_dir,'Cards','reweight_card.dat')
        self.paths['delphes'] = pjoin(self.me_dir,'Cards','delphes_card.dat')
        

    def __init__(self, question, cards=[], mode='auto', *args, **opt):

        self.load_default()        
        self.define_paths(**opt)
        cmd.OneLinePathCompletion.__init__(self, question, *args, **opt)


        try:
            self.param_card = check_param_card.ParamCard(self.paths['param'])
        except (check_param_card.InvalidParamCard, ValueError) as e:
            logger.error('Current param_card is not valid. We are going to use the default one.')
            logger.error('problem detected: %s' % e)
            files.cp(pjoin(self.me_dir,'Cards','param_card_default.dat'),
                     pjoin(self.me_dir,'Cards','param_card.dat'))
            self.param_card = check_param_card.ParamCard(pjoin(self.me_dir,'Cards','param_card.dat'))
        default_param = check_param_card.ParamCard(pjoin(self.me_dir,'Cards','param_card_default.dat'))
        self.param_card_default = default_param
        
        try:
            self.run_card = banner_mod.RunCard(self.paths['run'])
        except IOError:
            self.run_card = {}
        try:
            run_card_def = banner_mod.RunCard(self.paths['run_default'])
        except IOError:
            run_card_def = {}

        self.pname2block = {}
        self.conflict = []
        self.restricted_value = {}
        self.mode = mode
        self.cards = cards

        # Read the comment of the param_card_default to find name variable for
        # the param_card also check which value seems to be constrained in the
        # model.
        self.pname2block, self.restricted_value = \
                                              default_param.analyze_param_card()

        if run_card_def:
            self.run_set = run_card_def.keys() + self.run_card.hidden_param
        elif self.run_card:
            self.run_set = self.run_card.keys()
        else:
            self.run_set = []
        # check for conflict with run_card
        for var in self.pname2block:
            if var in self.run_set:
                self.conflict.append(var)        
                

        self.has_delphes = False        
        if 'delphes_card.dat' in cards:
            self.has_delphes = True

        #check if Madweight_card is present:
        self.has_mw = False
        if 'madweight_card.dat' in cards:
            
            self.do_change_tf = self.mother_interface.do_define_transfer_fct
            self.complete_change_tf = self.mother_interface.complete_define_transfer_fct
            self.help_change_tf = self.mother_interface.help_define_transfer_fct
            if not os.path.exists(self.paths['transfer']):
                logger.warning('No transfer function currently define. Please use the change_tf command to define one.')
            
            
            self.has_mw = True
            try:
                import madgraph.madweight.Cards as mwcards
            except:
                import internal.madweight.Cards as mwcards
            self.mw_card = mwcards.Card(self.paths['MadWeight'])
            self.mw_card = self.mw_card.info
            self.mw_vars = []
            for key in self.mw_card:
                if key == 'comment': 
                    continue
                for key2 in self.mw_card.info[key]:
                    if isinstance(key2, str) and not key2.isdigit():
                        self.mw_vars.append(key2)
            
            # check for conflict with run_card/param_card
            for var in self.pname2block:                
                if var in self.mw_vars:
                    self.conflict.append(var)           
            for var in self.mw_vars:
                if var in self.run_card:
                    self.conflict.append(var)
                    
        #check if MadLoopParams.dat is present:
        self.has_ml = False
        if os.path.isfile(self.paths['ML']):
            self.has_ml = True
            self.MLcard = banner_mod.MadLoopParam(self.paths['ML'])
            self.MLcardDefault = banner_mod.MadLoopParam()
            
            self.ml_vars = [k.lower() for k in self.MLcard.keys()]
            # check for conflict
            for var in self.ml_vars:
                if var in self.run_card:
                    self.conflict.append(var)
                if var in self.pname2block:
                    self.conflict.append(var)
                if self.has_mw and var in self.mw_vars:
                    self.conflict.append(var)

        #check if shower_card is present:
        self.has_shower = False
        if 'shower_card.dat' in cards:
            self.has_shower = True
            try:
                import madgraph.various.shower_card as showercards
            except:
                import internal.shower_card as showercards
            self.shower_card = showercards.ShowerCard(self.paths['shower'])
            self.shower_vars = self.shower_card.keys()
            
            # check for conflict with run_card/param_card
            for var in self.pname2block:                
                if var in self.shower_vars:
                    self.conflict.append(var)           
            for var in self.shower_vars:
                if var in self.run_card:
                    self.conflict.append(var)

        #check if pythia8_card.dat is present:
        self.has_PY8 = False
        if 'pythia8_card.dat' in cards:
            self.has_PY8 = True
            self.PY8Card = banner_mod.PY8Card(self.paths['PY8'])
            self.PY8CardDefault = banner_mod.PY8Card()
            
            self.py8_vars = [k.lower() for k in self.PY8Card.keys() if 
                                     k.lower() not in self.PY8Card.hidden_param]
            # check for conflict
            for var in self.py8_vars:
                if var in self.run_card:
                    self.conflict.append(var)
                if var in self.pname2block:
                    self.conflict.append(var)
                if self.has_mw and var in self.mw_vars:
                    self.conflict.append(var)
                if self.has_ml and var in self.ml_vars:
                    self.conflict.append(var)

    def do_help(self, line, conflict_raise=False, banner=True):    
     try:       
         
        if banner:                      
            logger.info('*** HELP MESSAGE ***', '$MG:color:BLACK')
         
        args = self.split_arg(line)
        # handle comand related help
        if len(args) == 1 and hasattr(self, 'do_%s' % args[0]):
            out = cmd.BasicCmd.do_help(self, line)
            if banner:
                logger.info('*** END HELP ***', '$MG:color:BLACK')  
            return out      
        # check for special shortcut.
        # special shortcut:
        if args[0] in self.special_shortcut:    
            if args[0] in self.special_shortcut_help:
                print self.special_shortcut_help[args[0]]
            if banner:
                logger.info('*** END HELP ***', '$MG:color:BLACK')  
            return       
        
        start = 0
        card = ''
        if args[0] in self.all_card_name:
            start += 1
            card = args[0]
            if len(args) == 1:
                if args[0] == 'param_card':
                    logger.info("Param_card information: ", '$MG:color:BLUE')
                    print "File to define the various model parameter"
                    logger.info("List of the Block defined:",'$MG:color:BLUE')
                    print "\t".join(self.param_card.keys())
                elif hasattr(self, args[0]):
                    logger.info("%s information: " % args[0], '$MG:color:BLUE')
                    print(eval('self.%s' % args[0]).__doc__)
                    logger.info("List of parameter associated", '$MG:color:BLUE')
                    print "\t".join(eval('self.%s' % args[0]).keys())
            if banner:
                logger.info('*** END HELP ***', '$MG:color:BLACK')  
            return 
                    
                
        #### RUN CARD
        if args[start] in [l.lower() for l in self.run_card.keys()] and card in ['', 'run_card']:
            if args[start] not in self.run_set:
                args[start] = [l for l in self.run_set if l.lower() == args[start]][0]

            if args[start] in self.conflict and not conflict_raise:
                conflict_raise = True
                logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:color:BLACK')
                logger.info('**   If not explicitely speficy this parameter  will modif the run_card file', '$MG:color:BLACK')

            self.run_card.do_help(args[start])
        ### PARAM_CARD WITH BLOCK NAME -----------------------------------------
        elif (args[start] in self.param_card or args[start] == 'width') \
                                                  and card in ['','param_card']:
            if args[start] in self.conflict and not conflict_raise:
                conflict_raise = True
                logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:color:BLACK')
                logger.info('**   If not explicitely speficy this parameter  will modif the param_card file', '$MG:color:BLACK')
                 
            if args[start] == 'width':
                args[start] = 'decay'
                
            if len(args) == start+1:
                self.param_card.do_help(args[start], tuple())
                key = None
            elif args[start+1] in self.pname2block:
                all_var = self.pname2block[args[start+1]]
                key = None
                for bname, lhaid in all_var:
                    if bname == args[start]:
                        key = lhaid
                        break
                else:
                    logger.warning('%s is not part of block "%s" but "%s". please correct.' %
                                    (args[start+1], args[start], bname))
            else:
                try:
                    key = tuple([int(i) for i in args[start+1:]])
                except ValueError:
                    logger.warning('Failed to identify LHA information')
                    return            
            
            if key in self.param_card[args[start]].param_dict:
                self.param_card.do_help(args[start], key, default=self.param_card_default)
            elif key:
                logger.warning('invalid information: %s not defined in the param_card' % (key,))
        # PARAM_CARD NO BLOCK NAME ---------------------------------------------
        elif args[start] in self.pname2block and card in ['','param_card']:                
            if args[start] in self.conflict and not conflict_raise:
                conflict_raise = True
                logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:color:BLACK')
                logger.info('**   If not explicitely speficy this parameter  will modif the param_card file', '$MG:color:BLACK')
                 
            all_var = self.pname2block[args[start]]
            for bname, lhaid in all_var:
                new_line = 'param_card %s %s %s' % (bname,
                   ' '.join([ str(i) for i in lhaid]), ' '.join(args[start+1:]))
                self.do_help(new_line, conflict_raise=True, banner=False) 
        else:
            print "no help available" 
          
        if banner:                      
            logger.info('*** END HELP ***', '$MG:color:BLACK')    
        #raw_input('press enter to quit the help')
        return        
     except Exception, error:
         import traceback
         traceback.print_exc()
         print error    

    def complete_help(self, text, line, begidx, endidx):
     try:
        possibilities = self.complete_set(text, line, begidx, endidx,formatting=False)
        if line == '':
            possibilities['Defined command'] = cmd.BasicCmd.completenames(self, text, line)#, begidx, endidx)
        return self.deal_multiple_categories(possibilities)
     except Exception, error:
         import traceback
         traceback.print_exc()
         print error


    def complete_set(self, text, line, begidx, endidx, formatting=True):
        """ Complete the set command"""

        prev_timer = signal.alarm(0) # avoid timer if any
        if prev_timer:
            nb_back = len(line)
            self.stdout.write('\b'*nb_back + '[timer stopped]\n')
            self.stdout.write(line)
            self.stdout.flush()

        possibilities = {}
        allowed = {}
        args = self.split_arg(line[0:begidx])
        if args[-1] in ['Auto', 'default']:
            return
        if len(args) == 1:
            allowed = {'category':'', 'run_card':'', 'block':'all', 'param_card':'','shortcut':''}
            if self.has_mw:
                allowed['madweight_card'] = ''
                allowed['mw_block'] = 'all'
            if self.has_shower:
                allowed['shower_card'] = ''
            if self.has_ml:
                allowed['madloop_card'] = ''
            if self.has_PY8:
                allowed['pythia8_card'] = ''
            if self.has_delphes:
                allowed['delphes_card'] = ''
                
        elif len(args) == 2:
            if args[1] == 'run_card':
                allowed = {'run_card':'default'}
            elif args[1] == 'param_card':
                allowed = {'block':'all', 'param_card':'default'}
            elif args[1] in self.param_card.keys():
                allowed = {'block':args[1]}
            elif args[1] == 'width':
                allowed = {'block': 'decay'}
            elif args[1] == 'MadWeight_card':
                allowed = {'madweight_card':'default', 'mw_block': 'all'}
            elif args[1] == 'MadLoop_card':
                allowed = {'madloop_card':'default'}
            elif args[1] == 'pythia8_card':
                allowed = {'pythia8_card':'default'}                
            elif self.has_mw and args[1] in self.mw_card.keys():
                allowed = {'mw_block':args[1]}
            elif args[1] == 'shower_card':
                allowed = {'shower_card':'default'}
            elif args[1] == 'delphes_card':
                allowed = {'delphes_card':'default'}
            else:
                allowed = {'value':''}
        else:
            start = 1
            if args[1] in  ['run_card', 'param_card', 'MadWeight_card', 'shower_card', 
                            'MadLoop_card','pythia8_card','delphes_card']:
                start = 2
            if args[-1] in self.pname2block.keys():
                allowed['value'] = 'default'
            elif args[start] in self.param_card.keys() or args[start] == 'width':
                if args[start] == 'width':
                    args[start] = 'decay'
                    
                if args[start+1:]:
                    allowed = {'block':(args[start], args[start+1:])}
                else:
                    allowed = {'block':args[start]}
            elif self.has_mw and args[start] in self.mw_card.keys():
                if args[start+1:]:
                    allowed = {'mw_block':(args[start], args[start+1:])}
                else:
                    allowed = {'mw_block':args[start]}     
            #elif len(args) == start +1:
            #        allowed['value'] = ''
            else: 
                allowed['value'] = ''

        if 'category' in allowed.keys():
            categories = ['run_card', 'param_card']
            if self.has_mw:
                categories.append('MadWeight_card')
            if self.has_shower:
                categories.append('shower_card')
            if self.has_ml:
                categories.append('MadLoop_card')
            if self.has_PY8:
                categories.append('pythia8_card')
            if self.has_delphes:
                categories.append('delphes_card')
            
            possibilities['category of parameter (optional)'] = \
                          self.list_completion(text, categories)
        
        if 'shortcut' in allowed.keys():
            possibilities['special values'] = self.list_completion(text, self.special_shortcut.keys()+['qcut', 'showerkt'])

        if 'run_card' in allowed.keys():
            opts = self.run_set
            if allowed['run_card'] == 'default':
                opts.append('default')

            possibilities['Run Card'] = self.list_completion(text, opts)

        if 'param_card' in allowed.keys():
            opts = self.pname2block.keys()
            if allowed['param_card'] == 'default':
                opts.append('default')
            possibilities['Param Card'] = self.list_completion(text, opts)
            
        if 'madweight_card' in allowed.keys():
            opts = self.mw_vars + [k for k in self.mw_card.keys() if k !='comment']
            if allowed['madweight_card'] == 'default':
                opts.append('default')
            possibilities['MadWeight Card'] = self.list_completion(text, opts)            

        if 'madloop_card' in allowed.keys():
            opts = self.ml_vars
            if allowed['madloop_card'] == 'default':
                opts.append('default')
            possibilities['MadLoop Parameter'] = self.list_completion(text, opts)
                                
        if 'pythia8_card' in allowed.keys():
            opts = self.py8_vars
            if allowed['pythia8_card'] == 'default':
                opts.append('default')
            possibilities['Pythia8 Parameter'] = self.list_completion(text, opts)
                                
        if 'shower_card' in allowed.keys():
            opts = self.shower_vars + [k for k in self.shower_card.keys() if k !='comment']
            if allowed['shower_card'] == 'default':
                opts.append('default')
            possibilities['Shower Card'] = self.list_completion(text, opts)            

        if 'delphes_card' in allowed:
            if allowed['delphes_card'] == 'default':
                opts = ['default', 'atlas', 'cms']
            possibilities['Delphes Card'] = self.list_completion(text, opts)              

        if 'value' in allowed.keys():
            opts = ['default']
            if 'decay' in args:
                opts.append('Auto')
                opts.append('Auto@NLO')
            elif args[-1] in self.pname2block and self.pname2block[args[-1]][0][0] == 'decay':
                opts.append('Auto')
                opts.append('Auto@NLO')
            possibilities['Special Value'] = self.list_completion(text, opts)

        if 'block' in allowed.keys():
            if allowed['block'] == 'all':
                allowed_block = [i for i in self.param_card.keys() if 'qnumbers' not in i]
                allowed_block.append('width')
                possibilities['Param Card Block' ] = \
                                       self.list_completion(text, allowed_block)
            elif isinstance(allowed['block'], basestring):
                block = self.param_card[allowed['block']].param_dict
                ids = [str(i[0]) for i in block
                          if (allowed['block'], i) not in self.restricted_value]
                possibilities['Param Card id' ] = self.list_completion(text, ids)
                varname = [name for name, all_var in self.pname2block.items()
                                               if any((bname == allowed['block']
                                                   for bname,lhaid in all_var))]
                possibilities['Param card variable'] = self.list_completion(text,
                                                                        varname)
            else:
                block = self.param_card[allowed['block'][0]].param_dict
                nb = len(allowed['block'][1])
                ids = [str(i[nb]) for i in block if len(i) > nb and \
                            [str(a) for a in i[:nb]] == allowed['block'][1]]

                if not ids:
                    if tuple([int(i) for i in allowed['block'][1]]) in block:
                        opts = ['default']
                        if allowed['block'][0] == 'decay':
                            opts.append('Auto')
                            opts.append('Auto@NLO')
                        possibilities['Special value'] = self.list_completion(text, opts)
                possibilities['Param Card id' ] = self.list_completion(text, ids)

        if 'mw_block' in allowed.keys():
            if allowed['mw_block'] == 'all':
                allowed_block = [i for i in self.mw_card.keys() if 'comment' not in i]
                possibilities['MadWeight Block' ] = \
                                       self.list_completion(text, allowed_block)
            elif isinstance(allowed['mw_block'], basestring):
                block = self.mw_card[allowed['mw_block']]
                ids = [str(i[0]) if isinstance(i, tuple) else str(i) for i in block]
                possibilities['MadWeight Card id' ] = self.list_completion(text, ids)
            else:
                block = self.mw_card[allowed['mw_block'][0]]
                nb = len(allowed['mw_block'][1])
                ids = [str(i[nb]) for i in block if isinstance(i, tuple) and\
                           len(i) > nb and \
                           [str(a) for a in i[:nb]] == allowed['mw_block'][1]]
                
                if not ids:
                    if tuple([i for i in allowed['mw_block'][1]]) in block or \
                                      allowed['mw_block'][1][0] in block.keys():
                        opts = ['default']
                        possibilities['Special value'] = self.list_completion(text, opts)
                possibilities['MadWeight Card id' ] = self.list_completion(text, ids) 

        return self.deal_multiple_categories(possibilities, formatting)

    def do_set(self, line):
        """ edit the value of one parameter in the card"""
        
        args = self.split_arg(line)
        # fix some formatting problem
        if len(args)==1 and '=' in args[-1]:
            arg1, arg2 = args.pop(-1).split('=',1)
            args += [arg1, arg2]
        if '=' in args:
            args.remove('=')

        # special shortcut:
        if args[0] in self.special_shortcut:
            if len(args) == 1:
                values = {}
            elif len(args) == 2: 
                targettype = float
                if args[1].strip().isdigit():
                    targettype = int
                 
                try:  
                    values = {'0': targettype(args[1])}
                except ValueError as e:
                    logger.warning("Wrong argument: The last entry should be a number.")
                    return
            else:
                logger.warning("too many argument for this command")
                return
            
            for arg in self.special_shortcut[args[0]]:
                try:
                    text = arg % values
                except KeyError:
                    logger.warning("This command requires one argument")
                    return
                except Exception as e:
                    logger.warning(str(e))
                    return
                else:
                    self.do_set(arg % values)
            return

        
        start = 0
        if len(args) < 2:
            logger.warning('Invalid set command %s (need two arguments)' % line)
            return            

        # Special case for the qcut value
        if args[0].lower() == 'qcut':
            pythia_path = self.paths['pythia']
            if os.path.exists(pythia_path):
                logger.info('add line QCUT = %s in pythia_card.dat' % args[1])
                p_card = open(pythia_path,'r').read()
                p_card, n = re.subn('''^\s*QCUT\s*=\s*[\de\+\-\.]*\s*$''',
                                    ''' QCUT = %s ''' % args[1], \
                                    p_card, flags=(re.M+re.I))
                if n==0:
                    p_card = '%s \n QCUT= %s' % (p_card, args[1])
                with open(pythia_path, 'w') as fsock: 
                    fsock.write(p_card)
                return
        # Special case for the showerkt value
        if args[0].lower() == 'showerkt':
            pythia_path = self.paths['pythia']
            if os.path.exists(pythia_path):
                logger.info('add line SHOWERKT = %s in pythia_card.dat' % args[1].upper())
                p_card = open(pythia_path,'r').read()
                p_card, n = re.subn('''^\s*SHOWERKT\s*=\s*[default\de\+\-\.]*\s*$''',
                                    ''' SHOWERKT = %s ''' % args[1].upper(), \
                                    p_card, flags=(re.M+re.I))
                if n==0:
                    p_card = '%s \n SHOWERKT= %s' % (p_card, args[1].upper())
                with open(pythia_path, 'w') as fsock:
                    fsock.write(p_card)
                return
            

        card = '' #store which card need to be modify (for name conflict)
        if args[0] == 'madweight_card':
            if not self.mw_card:
                logger.warning('Invalid Command: No MadWeight card defined.')
                return
            args[0] = 'MadWeight_card'
        
        if args[0] == 'shower_card':
            if not self.shower_card:
                logger.warning('Invalid Command: No Shower card defined.')
                return
            args[0] = 'shower_card'

        if args[0] == "madloop_card":
            if not self.has_ml:
                logger.warning('Invalid Command: No MadLoopParam card defined.')
                return
            args[0] = 'MadLoop_card'

        if args[0] == "pythia8_card":
            if not self.has_PY8:
                logger.warning('Invalid Command: No Pythia8 card defined.')
                return
            args[0] = 'pythia8_card'
            
        if args[0] == 'delphes_card':
            if not self.has_delphes:
                logger.warning('Invalid Command: No Delphes card defined.')
                return
            if args[1] == 'atlas':
                logger.info("set default ATLAS configuratin for Delphes")
                files.cp(pjoin(self.me_dir,'Cards', 'delphes_card_ATLAS.dat'),
                         pjoin(self.me_dir,'Cards', 'delphes_card.dat'))
                return
            elif args[1] == 'cms':
                logger.info("set default CMS configuratin for Delphes")
                files.cp(pjoin(self.me_dir,'Cards', 'delphes_card_CMS.dat'),
                         pjoin(self.me_dir,'Cards', 'delphes_card.dat'))
                return
            
            
        if args[0] in ['run_card', 'param_card', 'MadWeight_card', 'shower_card',
                       'delphes_card']:
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                files.cp(self.paths['%s_default' %args[0][:-5]], self.paths[args[0][:-5]])
                if args[0] == 'param_card':
                    self.param_card = check_param_card.ParamCard(self.paths['param'])
                elif args[0] == 'run_card':
                    self.run_card = banner_mod.RunCard(self.paths['run'])
                elif args[0] == 'shower_card':
                    self.shower_card = shower_card_mod.ShowerCard(self.paths['shower'])
                return
            else:
                card = args[0]
            start=1
            if len(args) < 3:
                logger.warning('Invalid set command: %s (not enough arguments)' % line)
                return
            
        elif args[0] in ['MadLoop_card']:
            if args[1] == 'default':
                logging.info('replace MadLoopParams.dat by the default card')
                self.MLcard = banner_mod.MadLoopParam(self.MLcardDefault)
                self.MLcard.write(self.paths['ML'],
                                  commentdefault=True)
                return
            else:
                card = args[0]
            start=1
            if len(args) < 3:
                logger.warning('Invalid set command: %s (not enough arguments)' % line)
                return
        elif args[0] in ['pythia8_card']:
            if args[1] == 'default':
                logging.info('replace pythia8_card.dat by the default card')
                self.PY8Card = banner_mod.PY8Card(self.PY8CardDefault)
                self.PY8Card.write(pjoin(self.me_dir,'Cards','pythia8_card.dat'),
                          pjoin(self.me_dir,'Cards','pythia8_card_default.dat'),
                          print_only_visible=True)
                return
            else:
                card = args[0]
            start=1
            if len(args) < 3:
                logger.warning('Invalid set command: %s (not enough arguments)' % line)
                return
        elif args[0] in ['madspin_card']:
            if args[1] == 'default':
                logging.info('replace madspin_card.dat by the default card')
                files.cp(self.paths['MS_default'], self.paths['madspin'])
                return
            else:
                logger.warning("""Command set not allowed for modifying the madspin_card. 
                    Check the command \"decay\" instead.""")
                return

        #### RUN CARD
        if args[start] in [l.lower() for l in self.run_card.keys()] and card in ['', 'run_card']:
            if args[start] not in self.run_set:
                args[start] = [l for l in self.run_set if l.lower() == args[start]][0]

            if args[start] in self.conflict and card == '':
                text  = 'Ambiguous name (present in more than one card). Will assume it to be referred to run_card.\n'
                text += 'If this is not intended, please reset it in the run_card and specify the relevant card to \n'
                text += 'edit, in the format < set card parameter value >'
                logger.warning(text)

            if args[start+1] == 'default':
                default = banner_mod.RunCard(self.paths['run_default'])
                if args[start] in default.keys():
                    self.setR(args[start],default[args[start]])
                else:
                    logger.info('remove information %s from the run_card' % args[start])
                    del self.run_card[args[start]]
            else:
                if args[0].startswith('sys_') or \
                   args[0] in self.run_card.list_parameter or \
                   args[0] in self.run_card.dict_parameter:
                    val = ' '.join(args[start+1:])
                    val = val.split('#')[0]
                else:
                    val = args[start+1]
                self.setR(args[start], val)
            self.run_card.write(self.paths['run'], self.paths['run_default'])
            
        ### PARAM_CARD WITH BLOCK NAME -----------------------------------------
        elif (args[start] in self.param_card or args[start] == 'width') \
                                                  and card in ['','param_card']:
            #special treatment for scan
            if any(t.startswith('scan') for t in args):
                index = [i for i,t in enumerate(args) if t.startswith('scan')][0]
                args = args[:index] + [' '.join(args[index:])]
                
            if args[start] in self.conflict and card == '':
                text  = 'ambiguous name (present in more than one card). Please specify which card to edit'
                text += ' in the format < set card parameter value>'
                logger.warning(text)
                return
            
            if args[start] == 'width':
                args[start] = 'decay'

            if args[start+1] in self.pname2block:
                all_var = self.pname2block[args[start+1]]
                key = None
                for bname, lhaid in all_var:
                    if bname == args[start]:
                        key = lhaid
                        break
                else:
                    logger.warning('%s is not part of block "%s" but "%s". please correct.' %
                                    (args[start+1], args[start], bname))
                    return
            else:
                try:
                    key = tuple([int(i) for i in args[start+1:-1]])
                except ValueError:
                    if args[start] == 'decay' and args[start+1:-1] == ['all']:
                        for key in self.param_card[args[start]].param_dict:
                            if (args[start], key) in self.restricted_value:
                                continue
                            else:
                                self.setP(args[start], key, args[-1])
                        self.param_card.write(self.paths['param'])
                        return
                    logger.warning('invalid set command %s (failed to identify LHA information)' % line)
                    return

            if key in self.param_card[args[start]].param_dict:
                if (args[start], key) in self.restricted_value:
                    text = "Note that this parameter seems to be ignore by MG.\n"
                    text += "MG will use instead the expression: %s\n" % \
                                      self.restricted_value[(args[start], key)]
                    text += "You need to match this expression for external program (such pythia)."
                    logger.warning(text)

                if args[-1].lower() in ['default', 'auto', 'auto@nlo'] or args[-1].startswith('scan'):
                    self.setP(args[start], key, args[-1])
                else:
                    try:
                        value = float(args[-1])
                    except Exception:
                        logger.warning('Invalid input: Expected number and not \'%s\'' \
                                                                     % args[-1])
                        return
                    self.setP(args[start], key, value)
            else:
                logger.warning('invalid set command %s' % line)
                return
            self.param_card.write(self.paths['param'])
        
        # PARAM_CARD NO BLOCK NAME ---------------------------------------------
        elif args[start] in self.pname2block and card in ['','param_card']:
            if args[start] in self.conflict and card == '':
                text  = 'ambiguous name (present in more than one card). Please specify which card to edit'
                text += ' in the format < set card parameter value>'
                logger.warning(text)
                return
            
            all_var = self.pname2block[args[start]]
            for bname, lhaid in all_var:
                new_line = 'param_card %s %s %s' % (bname,
                   ' '.join([ str(i) for i in lhaid]), ' '.join(args[start+1:]))
                self.do_set(new_line)
            if len(all_var) > 1:
                logger.warning('This variable correspond to more than one parameter in the param_card.')
                for bname, lhaid in all_var:
                    logger.warning('   %s %s' % (bname, ' '.join([str(i) for i in lhaid])))
                logger.warning('all listed variables have been modified')
                
        # MadWeight_card with block name ---------------------------------------
        elif self.has_mw and (args[start] in self.mw_card and args[start] != 'comment') \
                                              and card in ['','MadWeight_card']:
            
            if args[start] in self.conflict and card == '':
                text  = 'ambiguous name (present in more than one card). Please specify which card to edit'
                text += ' in the format < set card parameter value>'
                logger.warning(text)
                return
                       
            block = args[start]
            name = args[start+1]
            value = args[start+2:]
            self.setM(block, name, value)
            self.mw_card.write(self.paths['MadWeight'])        
        
        # MadWeight_card NO Block name -----------------------------------------
        elif self.has_mw and args[start] in self.mw_vars \
                                             and card in ['', 'MadWeight_card']:
            
            if args[start] in self.conflict and card == '':
                text  = 'ambiguous name (present in more than one card). Please specify which card to edit'
                text += ' in the format < set card parameter value>'
                logger.warning(text)
                return

            block = [b for b, data in self.mw_card.items() if args[start] in data]
            if len(block) > 1:
                logger.warning('%s is define in more than one block: %s.Please specify.'
                               % (args[start], ','.join(block)))
                return
           
            block = block[0]
            name = args[start]
            value = args[start+1:]
            self.setM(block, name, value)
            self.mw_card.write(self.paths['MadWeight'])
             
        # MadWeight_card New Block ---------------------------------------------
        elif self.has_mw and args[start].startswith('mw_') and len(args[start:]) == 3\
                                                    and card == 'MadWeight_card':
            block = args[start]
            name = args[start+1]
            value = args[start+2]
            self.setM(block, name, value)
            self.mw_card.write(self.paths['MadWeight'])    

        #### SHOWER CARD
        elif self.has_shower and args[start].lower() in [l.lower() for l in \
                       self.shower_card.keys()] and card in ['', 'shower_card']:
            if args[start] not in self.shower_card:
                args[start] = [l for l in self.shower_card if l.lower() == args[start].lower()][0]

            if args[start] in self.conflict and card == '':
                text  = 'ambiguous name (present in more than one card). Please specify which card to edit'
                text += ' in the format < set card parameter value>'
                logger.warning(text)
                return

            if args[start+1].lower() == 'default':
                default = shower_card_mod.ShowerCard(self.paths['shower_default'])
                if args[start] in default.keys():
                    self.shower_card.set_param(args[start],default[args[start]], self.paths['shower'])
                else:
                    logger.info('remove information %s from the shower_card' % args[start])
                    del self.shower_card[args[start]]
            elif args[start+1].lower() in ['t','.true.','true']:
                self.shower_card.set_param(args[start],'.true.',self.paths['shower'])
            elif args[start+1].lower() in ['f','.false.','false']:
                self.shower_card.set_param(args[start],'.false.',self.paths['shower'])
            elif args[start] in ['analyse', 'extralibs', 'extrapaths', 'includepaths'] or\
                                                  args[start].startswith('dm_'):
                #case sensitive parameters
                args = line.split()
                args_str = ' '.join(str(a) for a in args[start+1:len(args)])
                self.shower_card.set_param(args[start],args_str,pjoin(self.me_dir,'Cards','shower_card.dat'))
            else:
                args_str = ' '.join(str(a) for a in args[start+1:len(args)])
                self.shower_card.set_param(args[start],args_str,self.paths['shower'])
     
        # MadLoop Parameter  ---------------------------------------------------
        elif self.has_ml and args[start] in self.ml_vars \
                                               and card in ['', 'MadLoop_card']:
        
            if args[start] in self.conflict and card == '':
                text = 'ambiguous name (present in more than one card). Please specify which card to edit'
                logger.warning(text)
                return

            if args[start+1] == 'default':
                value = self.MLcardDefault[args[start]]
                default = True
            else:
                value = args[start+1]
                default = False
            self.setML(args[start], value, default=default)
            self.MLcard.write(self.paths['ML'],
                              commentdefault=True)

        # Pythia8 Parameter  ---------------------------------------------------
        elif self.has_PY8 and (card == 'pythia8_card' or (card == '' and \
             args[start] in self.PY8Card)):

            if args[start] in self.conflict and card == '':
                text = 'ambiguous name (present in more than one card). Please specify which card to edit'
                logger.warning(text)
                return

            if args[start+1] == 'default':
                value = self.PY8CardDefault[args[start]]
                default = True
            else:
                value = ' '.join(args[start+1:])
                default = False
            self.setPY8(args[start], value, default=default)
            self.PY8Card.write(pjoin(self.me_dir,'Cards','pythia8_card.dat'),
                          pjoin(self.me_dir,'Cards','pythia8_card_default.dat'),
                          print_only_visible=True)
                
        #INVALID --------------------------------------------------------------
        else:      
            logger.warning('invalid set command %s ' % line)
            return

    def setM(self, block, name, value):
        
        if isinstance(value, list) and len(value) == 1:
            value = value[0]
            
        if block not in self.mw_card:
            logger.warning('block %s was not present in the current MadWeight card. We are adding it' % block)
            self.mw_card[block] = {}
        elif name not in self.mw_card[block]:
            logger.info('name %s was not present in the block %s for the current MadWeight card. We are adding it' % (name,block),'$MG:color:BLACK')
        if value == 'default':
            import madgraph.madweight.Cards as mwcards
            mw_default = mwcards.Card(self.paths['MadWeight_default'])
            try:
                value = mw_default[block][name]
            except KeyError:
                logger.info('removing id "%s" from Block "%s" '% (name, block))
                if name in self.mw_card[block]:
                    del self.mw_card[block][name]
                return
        if value:
            logger.info('modify madweight_card information BLOCK "%s" with id "%s" set to %s' %\
                    (block, name, value), '$MG:color:BLACK')
        else:
            logger.value("Invalid command: No value. To set default value. Use \"default\" as value")
            return
        
        self.mw_card[block][name] = value
    
    def setR(self, name, value):
        logger.info('modify parameter %s of the run_card.dat to %s' % (name, value),'$MG:color:BLACK')
        self.run_card.set(name, value, user=True)

    def setML(self, name, value, default=False):
        
        try:
            self.MLcard.set(name, value, user=True)
        except Exception, error:
            logger.warning("Fail to change parameter. Please Retry. Reason: %s." % error)
            return
        logger.info('modify parameter %s of the MadLoopParam.dat to %s' % (name, value),'$MG:color:BLACK')
        if default and name.lower() in self.MLcard.user_set:
            self.MLcard.user_set.remove(name.lower())

    def setPY8(self, name, value, default=False):
        try:
            self.PY8Card.userSet(name, value)
        except Exception, error:
            logger.warning("Fail to change parameter. Please Retry. Reason: %s." % error)
            return
        logger.info('modify parameter %s of the pythia8_card.dat to %s' % (name, value), '$MG:color:BLACK')
        if default and name.lower() in self.PY8Card.user_set:
            self.PY8Card.user_set.remove(name.lower())

    def setP(self, block, lhaid, value):
        if isinstance(value, str):
            value = value.lower()
            if value == 'default':
                default = check_param_card.ParamCard(self.paths['param_default'])
                value = default[block].param_dict[lhaid].value

            elif value in ['auto', 'auto@nlo']:
                if 'nlo' in value:
                    value = 'Auto@NLO'
                else:
                    value = 'Auto'
                if block != 'decay':
                    logger.warning('Invalid input: \'Auto\' value only valid for DECAY')
                    return
            elif value.startswith('scan'):
                if ':' not in value:
                    logger.warning('Invalid input: \'scan\' mode requires a \':\' before the definition.')
                    return
                tag = value.split(':')[0]
                tag = tag[4:].strip()
                if tag and not tag.isdigit():
                    logger.warning('Invalid input: scan tag need to be integer and not "%s"' % tag)
                    return
                
                
                pass
            else:
                try:
                    value = float(value)
                except ValueError:
                    logger.warning('Invalid input: \'%s\' not valid intput.'% value)

        logger.info('modify param_card information BLOCK %s with id %s set to %s' %\
                    (block, lhaid, value), '$MG:color:BLACK')
        self.param_card[block].param_dict[lhaid].value = value
    
    def check_card_consistency(self):
        """This is run on quitting the class. Apply here all the self-consistency
        rule that you want. Do the modification via the set command."""

        # if NLO reweighting is ON: ensure that we keep the rwgt information
        if 'reweight' in self.allow_arg and 'run' in self.allow_arg and \
            isinstance(self.run_card,banner_mod.RunCardNLO) and \
            not self.run_card['store_rwgt_info']:
            #check if a NLO reweighting is required
                re_pattern = re.compile(r'''^\s*change\s*mode\s* (LO\+NLO|LO|NLO)\s*(?:#|$)''', re.M+re.I)
                text = open(self.paths['reweight']).read()
                options = re_pattern.findall(text)
                if any(o in ['NLO', 'LO+NLO'] for o in options):
                    logger.info('NLO reweighting is on ON. Automatically set store_rwgt_info to True', '$MG:color:BLACK' )
                    self.do_set('run_card store_rwgt_info True')
                    
        # Check the extralibs flag.
        if self.has_shower and isinstance(self.run_card, banner_mod.RunCardNLO):
            modify_extralibs, modify_extrapaths = False,False
            extralibs = self.shower_card['extralibs'].split()
            extrapaths = self.shower_card['extrapaths'].split()
            # remove default stdhep/Fmcfio for recent shower
            if self.run_card['parton_shower'] in ['PYTHIA8', 'HERWIGPP', 'HW7']:
                if 'stdhep' in self.shower_card['extralibs']:
                    extralibs.remove('stdhep')
                    modify_extralibs = True
                if 'Fmcfio' in self.shower_card['extralibs']:
                    extralibs.remove('Fmcfio')     
                    modify_extralibs = True               
            if self.run_card['parton_shower'] == 'PYTHIA8':
                # First check sanity of PY8
                if not self.mother_interface.options['pythia8_path']:
                    raise self.mother.InvalidCmd, 'Pythia8 is not correctly specified  to MadGraph5_aMC@NLO'
                executable = pjoin(self.mother_interface.options['pythia8_path'], 'bin', 'pythia8-config')
                if not os.path.exists(executable):
                    raise self.mother.InvalidCmd, 'Pythia8 is not correctly specified to MadGraph5_aMC@NLO'                
                
                # 2. take the compilation flag of PY8 from pythia8-config
                libs , paths = [], []
                p = misc.subprocess.Popen([executable, '--libs'], stdout=subprocess.PIPE)
                stdout, _ = p. communicate()
                libs = [x[2:] for x in stdout.split() if x.startswith('-l') or paths.append(x[2:])]
                
                # Add additional user-defined compilation flags
                p = misc.subprocess.Popen([executable, '--config'], stdout=subprocess.PIPE)
                stdout, _ = p. communicate()
                for lib in ['-ldl','-lstdc++','-lc++']:
                    if lib in stdout:
                        libs.append(lib[2:])                    

                # This precompiler flag is in principle useful for the analysis if it writes HEPMC
                # events, but there is unfortunately no way for now to specify it in the shower_card.
                supports_HEPMCHACK = '-DHEPMC2HACK' in stdout
                
                #3. ensure that those flag are in the shower card
                for l in libs:
                    if l not in extralibs:
                        modify_extralibs = True
                        extralibs.append(l)
                for L in paths:
                    if L not in extrapaths:
                        modify_extrapaths = True
                        extrapaths.append(L)
                        
            # Apply the required modification
            if modify_extralibs:
                if extralibs:
                    self.do_set('shower_card extralibs %s ' % ' '.join(extralibs))
                else:
                    self.do_set('shower_card extralibs \'\' ')
            if modify_extrapaths:
                if extrapaths:
                    self.do_set('shower_card extrapaths %s ' % ' '.join(extrapaths))
                else:
                    self.do_set('shower_card extrapaths \'\' ')    
    
    def reask(self, *args, **opt):
        
        cmd.OneLinePathCompletion.reask(self,*args, **opt)
        if self.has_mw and not os.path.exists(pjoin(self.me_dir,'Cards','transfer_card.dat')):
            logger.warning('No transfer function currently define. Please use the change_tf command to define one.')
    
    def postcmd(self, stop, line):
        ending_question = cmd.OneLinePathCompletion.postcmd(self,stop,line)

        if ending_question:
            self.check_card_consistency()
            self.do_update_dependent('', timer=20)
            return ending_question
    
    def do_update_dependent(self, line, timer=0):
        """Change the mass/width of particles which are not free parameter for the 
        model."""
        if not self.mother_interface:
            logger.warning('Failed to update dependent parameter. This might create trouble for external program (like MadSpin/shower/...)')
        
        pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
        pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M)
        param_text= open(self.paths['param']).read()
        
        if pattern_scan.search(param_text):
            #for block, key in self.restricted_value:
            #    self.param_card[block].get(key).value = -9.999e-99
            #    self.param_card.write(self.paths['param'])
            return
        elif pattern_width.search(param_text):
            self.do_compute_widths('')
            self.param_card = check_param_card.ParamCard(self.paths['param'])
    
        # calling the routine doing the work    
        self.update_dependent(self.mother_interface, self.me_dir, self.param_card,
                               self.paths['param'], timer)
    
    @staticmethod
    def update_dependent(mecmd, me_dir, param_card, path ,timer=0):
        """static method which can also be called from outside the class
           usefull in presence of scan.
           return if the param_card was updated or not
        """
        modify = True
        class TimeOutError(Exception): 
            pass
        def handle_alarm(signum, frame): 
            raise TimeOutError
        signal.signal(signal.SIGALRM, handle_alarm)
        if timer:
            signal.alarm(timer)
            log_level=30
        else:
            log_level=20
        # Try to load the model in the limited amount of time allowed
        try:
            model = mecmd.get_model()
            signal.alarm(0)
        except TimeOutError:
            logger.warning('The model takes too long to load so we bypass the updating of dependent parameter.\n'+\
                           'This might create trouble for external program (like MadSpin/shower/...)\n'+\
                           'The update can be forced without timer by typing \'update_dependent\' at the time of the card edition')
            modify =False
        except Exception:
            logger.warning('Failed to update dependent parameter. This might create trouble for external program (like MadSpin/shower/...)')
            signal.alarm(0)
        else:
            restrict_card = pjoin(me_dir,'Source','MODEL','param_card_rule.dat')
            if not os.path.exists(restrict_card):
                restrict_card = None
            #restrict_card = None
            if model:
                modify = param_card.update_dependent(model, restrict_card, log_level)
                if modify and path:
                    param_card.write(path)
            else:
                logger.warning('missing MG5aMC code. Fail to update dependent parameter. This might create trouble for program like MadSpin/shower/...')
            
        if log_level==20:
            logger.info('param_card up to date.')
            
        return modify
    
    
    
    
    def check_answer_consistency(self):
        """function called if the code reads a file"""
        self.check_card_consistency()
        self.do_update_dependent('', timer=20) 
      
    def help_set(self):
        '''help message for set'''

        logger.info('********************* HELP SET ***************************')
        logger.info("syntax: set [run_card|param_card|...] NAME [VALUE|default]")
        logger.info("syntax: set [param_card] BLOCK ID(s) [VALUE|default]")
        logger.info('')
        logger.info('-- Edit the param_card/run_card/... and replace the value of the')
        logger.info('    parameter by the value VALUE.')
        logger.info('   ')
        logger.info('-- Example:')
        logger.info('     set run_card ebeam1 4000')
        logger.info('     set ebeam2 4000')
        logger.info('     set lpp1 0')
        logger.info('     set ptj default')
        logger.info('')
        logger.info('     set param_card mass 6 175')
        logger.info('     set mass 25 125.3')
        logger.info('     set mass mh 125')
        logger.info('     set mh 125')
        logger.info('     set decay 25 0.004')
        logger.info('     set decay wh 0.004')
        logger.info('     set vmix 2 1 2.326612e-01')
        logger.info('')
        logger.info('     set param_card default #return all parameter to default')
        logger.info('     set run_card default')
        logger.info('********************* HELP SET ***************************')


    def default(self, line):
        """Default action if line is not recognized"""

        line = line.strip()
        args = line.split()
        if line == '' and self.default_value is not None:
            self.value = self.default_value
        # check if input is a file
        elif hasattr(self, 'do_%s' % args[0]):
            self.do_set(' '.join(args[1:]))
        elif os.path.isfile(line):
            self.copy_file(line)
            self.value = 'repeat'
        elif self.me_dir and os.path.exists(pjoin(self.me_dir, line)):
            self.copy_file(pjoin(self.me_dir,line))
            self.value = 'repeat'            
        elif line.strip() != '0' and line.strip() != 'done' and \
            str(line) != 'EOF' and line.strip() in self.allow_arg:            
            self.open_file(line)
            self.value = 'repeat'
        else:
            self.value = line

        return line

    def do_decay(self, line):
        """edit the madspin_card to define the decay of the associate particle"""
        signal.alarm(0) # avoid timer if any
        path = self.paths['madspin']
        
        if 'madspin_card.dat' not in self.cards or not os.path.exists(path):
            logger.warning("Command decay not valid. Since MadSpin is not available.")
            return
        
        if ">" not in line:
            logger.warning("invalid command for decay. Line ignored")
            return
        
        if "-add" in line:
            # just to have to add the line to the end of the file
            particle = line.split('>')[0].strip()
            text = open(path).read()
            line = line.replace('--add', '').replace('-add','')
            logger.info("change madspin_card to add one decay to %s: %s" %(particle, line.strip()), '$MG:color:BLACK')
            
            text = text.replace('launch', "\ndecay %s\nlaunch\n" % line,1)
        else:
            # Here we have to remove all the previous definition of the decay
            #first find the particle
            particle = line.split('>')[0].strip()
            logger.info("change madspin_card to define the decay of %s: %s" %(particle, line.strip()), '$MG:color:BLACK')
            particle = particle.replace('+','\+').replace('-','\-')
            decay_pattern = re.compile(r"^\s*decay\s+%s\s*>[\s\w+-~]*?$" % particle, re.I+re.M)
            text= open(path).read()
            text = decay_pattern.sub('', text)
            text = text.replace('launch', "\ndecay %s\nlaunch\n" % line,1)

        with open(path,'w') as fsock:
            fsock.write(text) 

        

    def do_compute_widths(self, line):
        signal.alarm(0) # avoid timer if any
        path = self.paths['param']
        pattern = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
        text = open(path).read()
        pdg_info = pattern.findall(text)
        has_nlo = any("@nlo"==nlo.lower() for _, nlo in pdg_info)
        pdg = [p for p,_ in pdg_info]
        
        
        line = '%s %s' % (line, ' '.join(pdg))
        if not '--path' in line:
            line += ' --path=%s' % path
        if has_nlo:
            line += ' --nlo'

        try:
            return self.mother_interface.do_compute_widths(line)
        except InvalidCmd, error:
            logger.error("Invalid command: %s " % error)

    def help_compute_widths(self):
        signal.alarm(0) # avoid timer if any
        return self.mother_interface.help_compute_widths()

    def help_decay(self):
        """help for command decay which modifies MadSpin_card"""
        
        signal.alarm(0) # avoid timer if any
        print '--syntax: decay PROC [--add]'
        print ' '
        print '  modify the madspin_card to modify the decay of the associate particle.'
        print '  and define it to PROC.'
        print '  if --add is present, just add a new decay for the associate particle.'
        
    def complete_compute_widths(self, *args, **opts):
        signal.alarm(0) # avoid timer if any
        return self.mother_interface.complete_compute_widths(*args,**opts)


    def help_add(self):
        """help for add command"""

        logger.info('********************* HELP ADD ***************************')
        logger.info( '-- syntax: add pythia8_card NAME VALUE')
        logger.info( "   add a definition of name in the pythia8_card with the given value")
        logger.info( "   Do not work for the param_card"        )
        logger.info( '-- syntax: add filename line')
        logger.info( '   add the given LINE to the end of the associate file (all file supportedd).')
        logger.info('********************* HELP ADD ***************************') 
    def complete_add(self, text, line, begidx, endidx):
        """ auto-completion for add command"""
        signal.alarm(0) # avoid timer if any
                 
        possibilities = {} 
        cards = [c.rsplit('.',1)[0] for c in self.cards]   
        possibilities['category of parameter (optional)'] = \
                          self.list_completion(text, cards)
                          
        return self.deal_multiple_categories(possibilities)
    
    def do_add(self, line):
        """ syntax: add filename NAME VALUE
            syntax: add filename LINE"""

        args = self.split_arg(line)
        if len(args) == 3 and args[0] in ['pythia8_card', 'pythia8_card.dat'] and self.has_PY8:
            name= args[1]
            value = args[2]
            self.PY8Card.userSet(name, value)
            self.PY8Card.write(pjoin(self.me_dir,'Cards','pythia8_card.dat'),
                          pjoin(self.me_dir,'Cards','pythia8_card_default.dat'),
                          print_only_visible=True)
            logger.info("add in the pythia8_card the parameter \"%s\" with value \"%s\"" % (name, value), '$MG:color:BLACK')
        elif len(args) > 0: 
            if args[0] in self.cards:
                card = args[0]
            elif "%s.dat" % args[0] in self.cards:
                card = "%s.dat" % args[0]
            elif "%s_card.dat" % args[0] in self.cards: 
                card = "%s_card.dat" % args[0]
            elif self.has_ml and args[0].lower() == "madloop":
                card = "MadLoopParams.dat"
            else:
                logger.error("unknow card %s. Please retry." % args[0])
                return
            
            ff = open(pjoin(self.me_dir,'Cards',card),'a')
            ff.write("%s \n" % line.split(None,1)[1])
            ff.close()
            self.reload_card(pjoin(self.me_dir,'Cards',card))
            logger.info("adding at the end of the file %s the line: \"%s\"" %(card, line.split(None,1)[1] ))


    def help_asperge(self):
        """Help associated to the asperge command"""
        signal.alarm(0)

        print '-- syntax: asperge [options]'
        print '   Call ASperGe to diagonalize all mass matrices in the model.'
        print '   This works only if the ASperGE module is part of the UFO model (a subdirectory).'
        print '   If you specify some names after the command (i.e. asperge m1 m2) then ASperGe will only'
        print '   diagonalize the associate mass matrices (here m1 and m2).'

    def complete_asperge(self, text, line, begidx, endidx, formatting=True):
        signal.alarm(0) # avoid timer if any

        blockname = self.pname2block.keys()
        # remove those that we know for sure are not mixing
        wrong = ['decay', 'mass', 'sminput']
        valid = [k for k in blockname if 'mix' in k]
        potential = [k for k in blockname if k not in valid+wrong]
        output = {'Mixing matrices': self.list_completion(text, valid, line),
                  'Other potential valid input': self.list_completion(text, potential, line)}

        return self.deal_multiple_categories(output, formatting)


    def do_asperge(self, line):
        """Running ASperGe"""
        signal.alarm(0) # avoid timer if any

        path = pjoin(self.me_dir,'bin','internal','ufomodel','ASperGE')
        if not os.path.exists(path):
            logger.error('ASperge has not been detected in the current model, therefore it will not be run.')
            return
        elif not os.path.exists(pjoin(path,'ASperGe')):
            logger.info('ASperGe has been detected but is not compiled. Running the compilation now.')
            try:
                misc.compile(cwd=path,shell=True)
            except MadGraph5Error, error:
                logger.error('''ASperGe failed to compile. Note that gsl is needed
     for this compilation to go trough. More information on how to install this package on
     http://www.gnu.org/software/gsl/
     Full compilation log is available at %s''' % pjoin(self.me_dir, 'ASperge_compilation.log'))
                open(pjoin(self.me_dir, 'ASperge_compilation.log'),'w').write(str(error))
                return

        opts = line.split()
        card = self.paths['param']
        logger.info('running ASperGE')
        returncode = misc.call([pjoin(path,'ASperGe'), card, '%s.new' % card] + opts)
        if returncode:
            logger.error('ASperGE fails with status %s' % returncode)
        else:
            logger.info('AsPerGe creates the file succesfully')
        files.mv(card, '%s.beforeasperge' % card)
        files.mv('%s.new' % card, card)



    def copy_file(self, path):
        """detect the type of the file and overwritte the current file"""
        
        if path.endswith('.lhco'):
            #logger.info('copy %s as Events/input.lhco' % (path))
            #files.cp(path, pjoin(self.mother_interface.me_dir, 'Events', 'input.lhco' ))
            self.do_set('mw_run inputfile %s' % os.path.relpath(path, self.mother_interface.me_dir))
            return
        elif path.endswith('.lhco.gz'):
            #logger.info('copy %s as Events/input.lhco.gz' % (path))
            #files.cp(path, pjoin(self.mother_interface.me_dir, 'Events', 'input.lhco.gz' ))
            self.do_set('mw_run inputfile %s' % os.path.relpath(path, self.mother_interface.me_dir))     
            return             
        else:
            card_name = CommonRunCmd.detect_card_type(path)

        if card_name == 'unknown':
            logger.warning('Fail to determine the type of the file. Not copied')
        if card_name != 'banner':
            logger.info('copy %s as %s' % (path, card_name))
            files.cp(path, self.paths[card_name.split('_',1)[0]])
            self.reload_card(self.paths[card_name.split('_',1)[0]])
        elif card_name == 'banner':
            banner_mod.split_banner(path, self.mother_interface.me_dir, proc_card=False)
            logger.info('Splitting the banner in it\'s component')
            if not self.mode == 'auto':
                self.mother_interface.keep_cards(self.cards)
            for card_name in self.cards:
                self.reload_card(pjoin(self.me_dir, 'Cards', card_name))

    def open_file(self, answer):
        """open the file"""
        try:
            me_dir = self.mother_interface.me_dir
        except:
            me_dir = None
            
        if answer.isdigit():
            if answer == '9':
                answer = 'plot'
            else:
                answer = self.cards[int(answer)-1]
        if 'madweight' in answer:
            answer = answer.replace('madweight', 'MadWeight')
        if 'MadLoopParams' in answer:
            answer = self.paths['ML']
        if 'pythia8_card' in answer:
            answer = self.paths['PY8']
        if not '.dat' in answer and not '.lhco' in answer:
            if answer != 'trigger':
                path = self.paths[answer]
            else:
                path = self.paths['delphes']
        elif not '.lhco' in answer:
            if '_' in answer:
                path = self.paths[answer.split('_')[0]]
            else:
                path = pjoin(me_dir, 'Cards', answer)
        else:
            path = pjoin(me_dir, self.mw_card['mw_run']['inputfile'])
            if not os.path.exists(path):
                logger.info('Path in MW_card not existing')
                path = pjoin(me_dir, 'Events', answer)
        #security
        path = path.replace('_card_card','_card')
        try:
            self.mother_interface.exec_cmd('open %s' % path)
        except InvalidCmd, error:
            if str(error) != 'No default path for this file':
                raise
            if answer == 'transfer_card.dat':
                logger.warning('You have to specify a transfer function first!')
            elif answer == 'input.lhco':
                path = pjoin(me_dir,'Events', 'input.lhco')
                ff = open(path,'w')
                ff.write('''No LHCO information imported at current time.
To import a lhco file: Close this file and type the path of your file.
You can also copy/paste, your event file here.''')
                ff.close()
                self.open_file(path)
            else:
                raise
        self.reload_card(path)
        
    def reload_card(self, path): 
        """reload object to have it in sync"""

        if path == self.paths['param']:        
            try:
                self.param_card = check_param_card.ParamCard(path) 
            except (check_param_card.InvalidParamCard, ValueError) as e:
                logger.error('Current param_card is not valid. We are going to use the default one.')
                logger.error('problem detected: %s' % e)
                logger.error('Please re-open the file and fix the problem.')
                logger.warning('using the \'set\' command without opening the file will discard all your manual change')
        elif path == self.paths['run']:
            self.run_card = banner_mod.RunCard(path)
        elif path == self.paths['ML']:
            self.MLcard = banner_mod.MadLoopParam(path)
        elif path == self.paths['PY8']
            # Use the read function so that modified/new parameters are correctly
            # set as 'user_set'
            if not self.PY8Card:
                self.PY8Card = banner_mod.PY8Card(self.paths['PY8_default'])

            self.PY8Card.read(self.paths['PY8'], setter='user')
            self.py8_vars = [k.lower() for k in self.PY8Card.keys() if 
                                     k.lower() not in self.PY8Card.hidden_param]
        elif path == self.paths['MadWeight']:
            try:
                import madgraph.madweight.Cards as mwcards
            except:
                import internal.madweight.Cards as mwcards
            self.mw_card = mwcards.Card(path)
        return path

class EditParamCard(AskforEditCard):
    """a dedicated module for the param"""
    
    special_shortcut ={}
    
    def __init__(self, question, card=[], mode='auto', *args, **opt):
                 
        self.load_default()
        cmd.OneLinePathCompletion.__init__(self, question, *args, **opt)
        if os.path.isfile(card[0]):
            self.param_card = check_param_card.ParamCard(card[0])
            self.paths['param'] = card[0]
            if os.path.isfile(card[0].replace('.dat', '_default.dat')):
                self.paths['param_default'] = card[0].replace('.dat', '_default.dat')
            else:
                self.paths['param_default'] = card[0]
        else:
            raise Exception, 'path %s do not exists' % card[0]
        
        self.pname2block, self.restricted_value = self.param_card.analyze_param_card()
        self.cards=['param']
        
    def do_asperge(self, *args, **opts):
        "Not available"
        logger.warning("asperge not available in this mode")

