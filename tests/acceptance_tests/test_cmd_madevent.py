################################################################################
#
# Copyright (c) 2009 The MadGraph Development team and Contributors
#
# This file is a part of the MadGraph 5 project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph license which should accompany this 
# distribution.
#
# For more information, please visit: http://madgraph.phys.ucl.ac.be
#
################################################################################
from __future__ import division
import subprocess
import unittest
import os
import re
import shutil
import sys
import logging
import time

logger = logging.getLogger('test_cmd')

import tests.unit_tests.iolibs.test_file_writers as test_file_writers

import madgraph.interface.master_interface as MGCmd
import madgraph.interface.madevent_interface as MECmd
import madgraph.interface.launch_ext_program as launch_ext

import madgraph.various.misc as misc
import madgraph.various.lhe_parser as lhe_parser

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
_pickle_path =os.path.join(_file_path, 'input_files')

from madgraph import MG4DIR, MG5DIR, MadGraph5Error, InvalidCmd

pjoin = os.path.join

    
    


#===============================================================================
# TestCmd
#===============================================================================
class TestMECmdShell(unittest.TestCase):
    """this treats all the command not related to MG_ME"""
    
    def generate(self, process, model):
        """Create a process"""

        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass

        interface = MGCmd.MasterCmd()
        interface.onecmd('import model %s' % model)
        if isinstance(process, str):
            interface.onecmd('generate %s' % process)
        else:
            for p in process:
                interface.onecmd('add process %s' % p)
        interface.onecmd('output madevent /tmp/MGPROCESS/ -f')
        pythia_path = pjoin(_file_path, os.path.pardir, 'pythia-pgs')
        if not os.path.exists(pythia_path):
            interface.onecmd('install pythia-pgs')

        misc.compile(cwd=pythia_path)
        if not misc.which('root'):
            raise Exception, 'root is require for this test'
        if not os.path.exists(pjoin(_file_path, os.path.pardir, 'MadAnalysis')):
            interface.onecmd('install MadAnalysis')
        
        self.cmd_line = MECmd.MadEventCmdShell(me_dir= '/tmp/MGPROCESS')
        self.cmd_line.exec_cmd('set automatic_html_opening False')

    @staticmethod
    def join_path(*path):
        """join path and treat spaces"""     
        combine = os.path.join(*path)
        return combine.replace(' ','\ ')        
    
    def do(self, line):
        """ exec a line in the cmd under test """        
        self.cmd_line.exec_cmd(line)
        
  
        
        
            
        
        
    def test_width_computation(self):
        """test the param_card created is correct"""
        
        cmd = os.getcwd()
        self.generate(['Z > l+ l-','Z > j j'], 'sm')
        self.assertEqual(cmd, os.getcwd())
        self.do('calculate_decay_widths -f')        
        
        # test the param_card is correctly written
        self.assertTrue(os.path.exists('/tmp/MGPROCESS/Events/run_01/param_card.dat'))
        
        text = open('/tmp/MGPROCESS/Events/run_01/param_card.dat').read()
        data = text.split('DECAY  23')[1].split('DECAY',1)[0]
        self.assertEqual("""1.492240e+00
#  BR             NDA  ID1    ID2   ...
   2.493165e-01   2    3  -3 # 0.37204
   2.493165e-01   2    1  -1 # 0.37204
   1.944158e-01   2    4  -4 # 0.290115
   1.944158e-01   2    2  -2 # 0.290115
   5.626776e-02   2    -11  11 # 0.083965
   5.626776e-02   2    -13  13 # 0.083965
#
#      PDG        Width""".split('\n'), data.strip().split('\n'))
        
    def test_creating_matched_plot(self):
        """test that the creation of matched plot works"""

        cmd = os.getcwd()
        self.generate('p p > W+', 'sm')
        self.assertEqual(cmd, os.getcwd())        
        shutil.copy(os.path.join(_file_path, 'input_files', 'run_card_matching.dat'),
                    '/tmp/MGPROCESS/Cards/run_card.dat')
        shutil.copy('/tmp/MGPROCESS/Cards/pythia_card_default.dat',
                    '/tmp/MGPROCESS/Cards/pythia_card.dat')
        self.do('generate_events -f')     
        
        
        f1 = self.check_matched_plot(tag='fermi')         
        start = time.time()
        
        #modify the run_card
        run_card = self.cmd_line.run_card
        run_card['nevents'] = 44
        run_card.write('/tmp/MGPROCESS/Cards/run_card.dat',
                                    '/tmp/MGPROCESS/Cards/run_card_default.dat')
            
        self.assertEqual(cmd, os.getcwd())        
        self.do('generate_events -f')
        self.assertEqual(int(self.cmd_line.run_card['nevents']), 44)
        self.do('pythia run_01 -f')
        self.do('quit')
        
        self.assertEqual(int(self.cmd_line.run_card['nevents']), 100)
        
        self.check_parton_output()
        self.check_parton_output('run_02', target_event=44)
        self.check_pythia_output()        
        f2 = self.check_matched_plot(mintime=start, tag='tag_1')        
        
        self.assertNotEqual(f1.split('\n'), f2.split('\n'))
        
        
        self.assertEqual(cmd, os.getcwd())
        
    def test_group_subprocess(self):
        """check that both u u > u u gives the same result"""
        
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass

        mg_cmd = MGCmd.MasterCmd()
        mg_cmd.exec_cmd('set automatic_html_opening False --save')
        mg_cmd.exec_cmd(' generate u u > u u')
        mg_cmd.exec_cmd('output /tmp/MGPROCESS/')
        self.cmd_line = MECmd.MadEventCmdShell(me_dir= '/tmp/MGPROCESS')
        self.cmd_line.exec_cmd('set automatic_html_opening False')
        
        self.do('generate_events -f')
        val1 = self.cmd_line.results.current['cross']
        err1 = self.cmd_line.results.current['error']
        
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass
        
        mg_cmd.exec_cmd('set group_subprocesses False')
        mg_cmd.exec_cmd('generate u u > u u')
        mg_cmd.exec_cmd('output /tmp/MGPROCESS')
        self.cmd_line = MECmd.MadEventCmdShell(me_dir= '/tmp/MGPROCESS')
        self.cmd_line.exec_cmd('set automatic_html_opening False')
        
        self.do('generate_events -f')        
        
        val2 = self.cmd_line.results.current['cross']
        err2 = self.cmd_line.results.current['error']        
        
        self.assertTrue(abs(val2 - val1) / (err1 + err2) < 5)
        target = 1278400
        self.assertTrue(abs(val2 - target) / (err2) < 5)
        #check precision
        self.assertTrue(err2 / val2 < 0.005)
        self.assertTrue(err1 / val1 < 0.005)
        
    def test_e_p_collision(self):
        """check that e p > e j gives the correct result"""
        
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass

        mg_cmd = MGCmd.MasterCmd()
        mg_cmd.exec_cmd('set automatic_html_opening False --save')
        mg_cmd.exec_cmd(' generate e- p  > e- j')
        mg_cmd.exec_cmd('output /tmp/MGPROCESS/')
        self.cmd_line = MECmd.MadEventCmdShell(me_dir= '/tmp/MGPROCESS')
        self.cmd_line.exec_cmd('set automatic_html_opening False')
        shutil.copy(os.path.join(_file_path, 'input_files', 'run_card_ep.dat'),
                    '/tmp/MGPROCESS/Cards/run_card.dat')
        
        self.do('generate_events -f')
        val1 = self.cmd_line.results.current['cross']
        err1 = self.cmd_line.results.current['error']
        
        target = 3864.0
        self.assertTrue(abs(val1 - target) / err1 < 1.)
        
    def test_e_e_collision(self):
        """check that e+ e- > t t~ gives the correct result"""
        
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass

        mg_cmd = MGCmd.MasterCmd()
        mg_cmd.exec_cmd('set automatic_html_opening False --save')
        mg_cmd.exec_cmd(' generate e+ e-  > t t~')
        mg_cmd.exec_cmd('output /tmp/MGPROCESS/')
        self.cmd_line = MECmd.MadEventCmdShell(me_dir= '/tmp/MGPROCESS')
        self.cmd_line.exec_cmd('set automatic_html_opening False')
        shutil.copy(os.path.join(_file_path, 'input_files', 'run_card_ee.dat'),
                    '/tmp/MGPROCESS/Cards/run_card.dat')
        
        self.do('generate_events -f')
        val1 = self.cmd_line.results.current['cross']
        err1 = self.cmd_line.results.current['error']
        
        target = 0.545
        self.assertTrue(abs(val1 - target) / err1 < 1.)
        
    def load_result(self, run_name):
        
        import madgraph.iolibs.save_load_object as save_load_object
        import madgraph.various.gen_crossxhtml as gen_crossxhtml
        
        result = save_load_object.load_from_file('/tmp/MGPROCESS/HTML/results.pkl')
        return result[run_name]

    def check_parton_output(self, run_name='run_01', target_event=100):
        """Check that parton output exists and reach the targert for event"""
                
        # check that the number of event is fine:
        data = self.load_result(run_name)
        self.assertEqual(int(data[0]['nb_event']), target_event)
        self.assertTrue('lhe' in data[0].parton)
                
    def check_pythia_output(self, run_name='run_01'):
        """ """
        # check that the number of event is fine:
        data = self.load_result(run_name)
        self.assertTrue('lhe' in data[0].pythia)
        self.assertTrue('log' in data[0].pythia)

    def check_matched_plot(self, run_name='run_01', mintime=None, tag='fermi'):
        """ """
        path = '/tmp/MGPROCESS/HTML/%(run)s/plots_pythia_%(tag)s/DJR1.ps' % \
                                {'run': run_name, 'tag': tag}
        self.assertTrue(os.path.exists(path))
        
        if mintime:
            self.assertTrue(os.path.getctime(path) > mintime)
        
        return open(path).read()
#===============================================================================
# TestCmd
#===============================================================================
class TestMEfromfile(unittest.TestCase):
    """test that we can launch everything from a single file"""


    def test_add_time_of_flight(self):
        """checking time of flight is working fine"""

        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass
        
        cmd = """import model sm
                set automatic_html_opening False --no-save
                 generate p p > w+ z
                 output /tmp/MGPROCESS -f -nojpeg
                 launch -i 
                 generate_events
                 parton
                 set nevents 100
                 add_time_of_flight --threshold=3e-26
                 pythia
                 """
        open('/tmp/mg5_cmd','w').write(cmd)
        
        devnull =open(os.devnull,'w')
        subprocess.call([pjoin(_file_path, os.path.pardir,'bin','mg5'), 
                         '/tmp/mg5_cmd'],
                         cwd=pjoin(_file_path, os.path.pardir),
                        stdout=devnull, stderr=devnull)

        self.check_parton_output(cross=17.21, error=0.19)
        self.check_pythia_output()
        event = '/tmp/MGPROCESS/Events/run_01/unweighted_events.lhe'
        if not os.path.exists(event):
            os.system('gunzip %s.gz' % event)
        
        has_zero = False
        has_non_zero = False
        for event in lhe_parser.EventFile(event):
            for particle in event:
                if particle.pid in [23,25]:
                    self.assertTrue(particle.vtim ==0 or particle.vtim > 3e-26)
                    if particle.vtim == 0 :
                        has_zero = True
                    else:
                        has_non_zero = True
        self.assertTrue(has_zero)
        self.assertTrue(has_non_zero)
        
        
        


    def test_generation_from_file_1(self):
        """ """
        cwd = os.getcwd()
        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass
        import subprocess
        
        devnull =open(os.devnull,'w')
        pythia_path =  pjoin(_file_path, os.path.pardir, 'pythia-pgs')
        if not os.path.exists(pythia_path):
            p = subprocess.Popen([pjoin(_file_path, os.path.pardir,'bin','mg5')],
                             stdin=subprocess.PIPE,
                             stdout=devnull,stderr=devnull)
            out = p.communicate('install pythia-pgs')
        misc.compile(cwd=pythia_path)
        

        subprocess.call([pjoin(_file_path, os.path.pardir,'bin','mg5'), 
                         pjoin(_file_path, 'input_files','test_mssm_generation')],
                         cwd=pjoin(_file_path, os.path.pardir),
                        stdout=devnull,stderr=devnull)

        
        self.check_parton_output(cross=4.541638, error=0.035)
        self.check_parton_output('run_02', cross=4.541638, error=0.035)
        self.check_pythia_output()
        self.assertEqual(cwd, os.getcwd())
        #

    def load_result(self, run_name):
        
        import madgraph.iolibs.save_load_object as save_load_object
        import madgraph.various.gen_crossxhtml as gen_crossxhtml
        
        result = save_load_object.load_from_file('/tmp/MGPROCESS/HTML/results.pkl')
        return result[run_name]

    def check_parton_output(self, run_name='run_01', target_event=100, cross=0, error=9e99):
        """Check that parton output exists and reach the targert for event"""
                
        # check that the number of event is fine:
        data = self.load_result(run_name)
        self.assertEqual(int(data[0]['nb_event']), target_event)
        self.assertTrue('lhe' in data[0].parton)
        
        if cross:
            self.assertTrue(abs(cross - float(data[0]['cross']))/error < 3)
                
    def check_pythia_output(self, run_name='run_01'):
        """ """
        # check that the number of event is fine:
        data = self.load_result(run_name)
        self.assertTrue('lhe' in data[0].pythia)
        self.assertTrue('log' in data[0].pythia)
        
        

#===============================================================================
# TestCmd
#===============================================================================
class TestMEfromPdirectory(unittest.TestCase):
    """test that we can launch everything from the P directory"""

    

    def generate(self, process, model):
        """Create a process"""

        try:
            shutil.rmtree('/tmp/MGPROCESS/')
        except Exception, error:
            pass

        interface = MGCmd.MasterCmd()
        interface.onecmd('import model %s' % model)
        if isinstance(process, str):
            interface.onecmd('generate %s' % process)
        else:
            for p in process:
                interface.onecmd('add process %s' % p)
        interface.onecmd('set automatic_html_opening False')
        interface.onecmd('output madevent /tmp/MGPROCESS/ -f')

    def load_result(self, run_name):
        
        import madgraph.iolibs.save_load_object as save_load_object
        import madgraph.various.gen_crossxhtml as gen_crossxhtml
        
        result = save_load_object.load_from_file('/tmp/MGPROCESS/HTML/results.pkl')
        return result[run_name]

    def check_parton_output(self, run_name='run_01', target_event=100, cross=0):
        """Check that parton output exists and reach the targert for event"""
                
        # check that the number of event is fine:
        data = self.load_result(run_name)
        self.assertEqual(int(data[0]['nb_event']), target_event)
        self.assertTrue('lhe' in data[0].parton)
        if cross:
            self.assertTrue(abs(cross - float(data[0]['cross']))/float(data[0]['error']) < 3)


        
    def test_run_fromP(self):
        """ """
                
        cmd = os.getcwd()
        self.generate('p p > e+ e-', 'sm')
        self.assertEqual(cmd, os.getcwd())
        shutil.copy(os.path.join(_file_path, 'input_files', 'run_card_matching.dat'),
                    '/tmp/MGPROCESS/Cards/run_card.dat')
        os.chdir('/tmp/MGPROCESS/')
        ff = open('cmd.cmd','w')
        ff.write('set automatic_html_opening False\n')
        ff.write('generate_events -f \n') 
        ff.close()
        if logger.getEffectiveLevel() > 20:
            output = open(os.devnull,'w')
        else:
            output = None
        id = subprocess.call(['./bin/madevent','cmd.cmd'], stdout=output, stderr=output)
        self.assertEqual(id, 0)
        self.check_parton_output(cross=947.9)
        os.chdir(cmd)
        
