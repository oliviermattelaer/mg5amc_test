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

"""Unit test library for the export v4 format routines"""

import StringIO
import copy
import fractions
import os 
import sys

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(os.path.join(root_path, os.path.pardir, os.path.pardir))

import tests.unit_tests as unittest

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.files as files
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.save_load_object as save_load_object
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.core.base_objects as base_objects
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.color_algebra as color
import madgraph.various.diagram_symmetry as diagram_symmetry
import madgraph.various.process_checks as process_checks
import madgraph.various.misc as misc
import madgraph.core.color_amp as color_amp
import madgraph.color_ordering.color_ordered_amplitudes as \
       color_ordered_amplitudes
import madgraph.color_ordering.color_ordered_helas_objects as \
       color_ordered_helas_objects
import madgraph.color_ordering.color_ordered_export_v4 as \
       color_ordered_export_v4
import tests.unit_tests.core.test_helas_objects as test_helas_objects
import tests.unit_tests.iolibs.test_file_writers as test_file_writers
import tests.unit_tests.iolibs.test_helas_call_writers as \
                                            test_helas_call_writers

_file_path = os.path.dirname(os.path.realpath(__file__))
_input_file_path = os.path.join(_file_path, os.path.pardir, os.path.pardir,
                                'input_files')
#===============================================================================
# COExportV4Test
#===============================================================================
class COExportV4Test(unittest.TestCase,
                     test_file_writers.CheckFileCreate):
    """Test class for the color ordered export v4 module"""

    mymodel = base_objects.Model()
    myfortranmodel = color_ordered_export_v4.COFortranUFOHelasCallWriter(mymodel)
    created_files = ['test'
                    ]

    def setUp(self):

        test_file_writers.CheckFileCreate.clean_files
        # Set up model

        self.mypartlist = base_objects.ParticleList()
        self.myinterlist = base_objects.InteractionList()

        # A gluon
        self.mypartlist.append(base_objects.Particle({'name':'g',
                      'antiname':'g',
                      'spin':3,
                      'color':8,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'g',
                      'antitexname':'g',
                      'line':'curly',
                      'charge':0.,
                      'pdg_code':21,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))

        # A quark U and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name':'u',
                      'antiname':'u~',
                      'spin':2,
                      'color':3,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'u',
                      'antitexname':'\bar u',
                      'line':'straight',
                      'charge':2. / 3.,
                      'pdg_code':2,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        antiu = copy.copy(self.mypartlist[1])
        antiu.set('is_part', False)

        # A quark D and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name':'d',
                      'antiname':'d~',
                      'spin':2,
                      'color':3,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'d',
                      'antitexname':'\bar d',
                      'line':'straight',
                      'charge':-1. / 3.,
                      'pdg_code':1,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        antid = copy.copy(self.mypartlist[2])
        antid.set('is_part', False)

        # A quark S and its antiparticle
        self.mypartlist.append(base_objects.Particle({'name':'s',
                      'antiname':'s~',
                      'spin':2,
                      'color':3,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'s',
                      'antitexname':'\bar s',
                      'line':'straight',
                      'charge':-1. / 3.,
                      'pdg_code':3,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        s = self.mypartlist[len(self.mypartlist) - 1]
        antis = copy.copy(s)
        antis.set('is_part', False)

        # A photon
        self.mypartlist.append(base_objects.Particle({'name':'a',
                      'antiname':'a',
                      'spin':3,
                      'color':1,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'\gamma',
                      'antitexname':'\gamma',
                      'line':'wavy',
                      'charge':0.,
                      'pdg_code':22,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))
        gamma = self.mypartlist[len(self.mypartlist) - 1]

        # A electron and positron
        self.mypartlist.append(base_objects.Particle({'name':'e+',
                      'antiname':'e-',
                      'spin':2,
                      'color':1,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'e^+',
                      'antitexname':'e^-',
                      'line':'straight',
                      'charge':-1.,
                      'pdg_code':11,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))
        e = self.mypartlist[len(self.mypartlist) - 1]
        antie = copy.copy(e)
        antie.set('is_part', False)

        # W
        self.mypartlist.append(base_objects.Particle({'name':'w+',
                      'antiname':'w-',
                      'spin':3,
                      'color':0,
                      'mass':'WMASS',
                      'width':'WWIDTH',
                      'texname':'W^+',
                      'antitexname':'W^-',
                      'line':'wavy',
                      'charge':1.,
                      'pdg_code':24,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':False}))

        wplus = self.mypartlist[len(self.mypartlist) - 1]
        wminus = copy.copy(wplus)
        wminus.set('is_part', False)

        # 3 gluon vertex
        self.myinterlist.append(base_objects.Interaction({
                      'id': 1,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[0]] * 3),
                      'color': [color.ColorString([color.f(0,1,2)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'G'},
                      'orders':{'QCD':1}}))

        # 4 gluon vertex
        self.myinterlist.append(base_objects.Interaction({
                      'id': 2,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[0]] * 4),
                      'color': [color.ColorString([color.f(1, 2, -1),
                                                   color.f(-1, 0, 3)]),
                                color.ColorString([color.f(1, 3, -1),
                                                   color.f(-1, 0, 2)]),
                                color.ColorString([color.f(2, 3, -1),
                                                   color.f(-1, 0, 1)])],
                      'lorentz':['VVVV4', 'VVVV3', 'VVVV1'],
            'couplings':{(0, 0):'GG', (1, 1):'GG', (2, 2):'GG'},
                      'orders':{'QCD':2}}))

        # Gluon and photon couplings to quarks
        self.myinterlist.append(base_objects.Interaction({
                      'id': 3,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[1], \
                                             antiu, \
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.T(2,0,1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'QCD':1}}))

        self.myinterlist.append(base_objects.Interaction({
                      'id': 4,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[1], \
                                             antiu, \
                                             gamma]),
                      'color': [color.ColorString([color.T(0,1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        self.myinterlist.append(base_objects.Interaction({
                      'id': 5,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[2], \
                                             antid, \
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.T(2,0,1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'QCD':1}}))

        self.myinterlist.append(base_objects.Interaction({
                      'id': 6,
                      'particles': base_objects.ParticleList(\
                                            [self.mypartlist[2], \
                                             antid, \
                                             gamma]),
                      'color': [color.ColorString([color.T(0,1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        self.myinterlist.append(base_objects.Interaction({
                      'id': 7,
                      'particles': base_objects.ParticleList(\
                                            [s, 
                                             antis,
                                             self.mypartlist[0]]),
                      'color': [color.ColorString([color.T(2,0,1)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'QCD':1}}))

        # Coupling of e to gamma

        self.myinterlist.append(base_objects.Interaction({
                      'id': 9,
                      'particles': base_objects.ParticleList(\
                                            [e, \
                                             antie, \
                                             gamma]),
                      'color': [color.ColorString([color.T(1,0)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        # Coupling of u and d to W

        self.myinterlist.append(base_objects.Interaction({
                      'id': 10,
                      'particles': base_objects.ParticleList(\
                                            [antid, \
                                             self.mypartlist[1], \
                                             wminus]),
                      'color': [color.ColorString([color.T(1,0)])],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        # Coupling of d and u to W

        self.myinterlist.append(base_objects.Interaction({
                      'id': 11,
                      'particles': base_objects.ParticleList(\
                                            [antiu, \
                                             self.mypartlist[2], \
                                             wplus]),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        self.mymodel.set('particles', self.mypartlist)
        self.mymodel.set('interactions', self.myinterlist)

        self.myfortranmodel.downcase = False

    tearDown = test_file_writers.CheckFileCreate.clean_files

    def test_export_co_matrix_element_v4_standalone_uux_uuxgg(self):
        """Test the result of exporting u u~ > u u~ g g to standalone"""

        # Create the amplitude
        myleglist = base_objects.LegList([\
            base_objects.Leg({'id':2, 'state':False}),
            base_objects.Leg({'id':-2, 'state':False}),
            base_objects.Leg({'id':2, 'state':True}),
            base_objects.Leg({'id':-2, 'state':True}),
            base_objects.Leg({'id':21, 'state':True}),
            base_objects.Leg({'id':21, 'state':True})])

        myproc = base_objects.Process({'legs':myleglist,
                                       'orders':{'QED': 0},
                                       'model':self.mymodel})

        self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

        matrix_element = color_ordered_helas_objects.COHelasMatrixElement(\
            self.myamplitude, gen_color=3, optimization=1)

        process_exporter = color_ordered_export_v4.ProcessExporterFortranCOSA()

        goal_flow_f = \
    ["""      COMPLEX*16 FUNCTION FLOW1(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 1
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=11)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=16)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_3(W(1,4),W(1,3),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_3(W(1,1),W(1,2),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,6),W(1,8),G,ZERO, ZERO, W(1,9))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,10))
      CALL L1_2(W(1,4),W(1,8),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_1(W(1,3),W(1,5),GQQ,ZERO, ZERO, W(1,12))
      CALL L1_2(W(1,1),W(1,6),GQQ,ZERO, ZERO, W(1,13))
      CALL L1_3(W(1,13),W(1,2),GQQ,ZERO, ZERO, W(1,14))
      CALL L1_2(W(1,13),W(1,5),GQQ,ZERO, ZERO, W(1,15))
      CALL L1_2(W(1,1),W(1,10),GQQ,ZERO, ZERO, W(1,16))
C     Amplitude(s) for diagram number 1
      CALL VVVV4_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(1))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(2))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,5),W(1,7),W(1,9),G,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,10),W(1,7),W(1,8),G,AMP(4))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,11),W(1,12),W(1,6),GQQ,AMP(5))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,4),W(1,12),W(1,9),GQQ,AMP(6))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,11),W(1,3),W(1,10),GQQ,AMP(7))
C     Amplitude(s) for diagram number 7
      CALL L1_0(W(1,5),W(1,7),W(1,14),G,AMP(8))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,15),W(1,2),W(1,7),GQQ,AMP(9))
C     Amplitude(s) for diagram number 9
      CALL L1_0(W(1,4),W(1,12),W(1,14),GQQ,AMP(10))
C     Amplitude(s) for diagram number 10
      CALL L1_0(W(1,16),W(1,2),W(1,7),GQQ,AMP(11))
      JAMP(1)=+1./2.*(-AMP(1)+AMP(2)-AMP(3)+AMP(4)-AMP(5)+IMAG1*AMP(6)
     $ -IMAG1*AMP(7)-IMAG1*AMP(8)-AMP(9)-AMP(10)-IMAG1*AMP(11))

      FLOW1=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info(),
     """      COMPLEX*16 FUNCTION FLOW2(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 2
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=12)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=16)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_3(W(1,4),W(1,3),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_3(W(1,1),W(1,2),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,5),W(1,8),G,ZERO, ZERO, W(1,9))
      CALL L1_1(W(1,6),W(1,8),G,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,3),W(1,6),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_2(W(1,4),W(1,5),GQQ,ZERO, ZERO, W(1,12))
      CALL L1_1(W(1,2),W(1,5),GQQ,ZERO, ZERO, W(1,13))
      CALL L1_2(W(1,1),W(1,6),GQQ,ZERO, ZERO, W(1,14))
      CALL L1_3(W(1,14),W(1,2),GQQ,ZERO, ZERO, W(1,15))
      CALL L1_3(W(1,1),W(1,13),GQQ,ZERO, ZERO, W(1,16))
C     Amplitude(s) for diagram number 1
      CALL VVVV4_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(1))
      CALL VVVV3_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(2))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,6),W(1,7),W(1,9),G,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,5),W(1,7),W(1,10),G,AMP(4))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,4),W(1,11),W(1,9),GQQ,AMP(5))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,12),W(1,11),W(1,8),GQQ,AMP(6))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,12),W(1,3),W(1,10),GQQ,AMP(7))
C     Amplitude(s) for diagram number 7
      CALL L1_0(W(1,14),W(1,13),W(1,7),GQQ,AMP(8))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,5),W(1,7),W(1,15),G,AMP(9))
C     Amplitude(s) for diagram number 9
      CALL L1_0(W(1,12),W(1,3),W(1,15),GQQ,AMP(10))
C     Amplitude(s) for diagram number 10
      CALL L1_0(W(1,6),W(1,7),W(1,16),G,AMP(11))
C     Amplitude(s) for diagram number 11
      CALL L1_0(W(1,4),W(1,11),W(1,16),GQQ,AMP(12))
      JAMP(1)=+1./2.*(+AMP(1)+AMP(2)+AMP(3)+AMP(4)-IMAG1*AMP(5)-AMP(6)
     $ +IMAG1*AMP(7)-AMP(8)+IMAG1*AMP(9)-AMP(10)-IMAG1*AMP(11)-AMP(12))

      FLOW2=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info(),
     """      COMPLEX*16 FUNCTION FLOW3(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 3
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=11)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=15)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_3(W(1,4),W(1,3),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_3(W(1,1),W(1,2),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,5),W(1,8),G,ZERO, ZERO, W(1,9))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,3),W(1,8),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_2(W(1,4),W(1,6),GQQ,ZERO, ZERO, W(1,12))
      CALL L1_1(W(1,2),W(1,5),GQQ,ZERO, ZERO, W(1,13))
      CALL L1_3(W(1,1),W(1,13),GQQ,ZERO, ZERO, W(1,14))
      CALL L1_2(W(1,1),W(1,7),GQQ,ZERO, ZERO, W(1,15))
C     Amplitude(s) for diagram number 1
      CALL VVVV3_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(1))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(2))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,6),W(1,7),W(1,9),G,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,10),W(1,7),W(1,8),G,AMP(4))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,12),W(1,11),W(1,5),GQQ,AMP(5))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,12),W(1,3),W(1,9),GQQ,AMP(6))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,4),W(1,11),W(1,10),GQQ,AMP(7))
C     Amplitude(s) for diagram number 7
      CALL L1_0(W(1,6),W(1,7),W(1,14),G,AMP(8))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,15),W(1,13),W(1,6),GQQ,AMP(9))
C     Amplitude(s) for diagram number 9
      CALL L1_0(W(1,12),W(1,3),W(1,14),GQQ,AMP(10))
C     Amplitude(s) for diagram number 10
      CALL L1_0(W(1,15),W(1,2),W(1,10),GQQ,AMP(11))
      JAMP(1)=+1./2.*(-AMP(1)-AMP(2)-AMP(3)-AMP(4)-AMP(5)-IMAG1*AMP(6)
     $ -IMAG1*AMP(7)+IMAG1*AMP(8)-AMP(9)-AMP(10)-IMAG1*AMP(11))

      FLOW3=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info(),
     """      COMPLEX*16 FUNCTION FLOW4(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 4
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=2 QED=0 WEIGHTED=4 singlet_QCD=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=5)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=13)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_2(W(1,1),W(1,6),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_2(W(1,7),W(1,5),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_3(W(1,4),W(1,2),GQQ,ZERO, ZERO, W(1,9))
      CALL L1_1(W(1,3),W(1,5),GQQ,ZERO, ZERO, W(1,10))
      CALL L1_2(W(1,1),W(1,9),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,12))
      CALL L1_2(W(1,1),W(1,12),GQQ,ZERO, ZERO, W(1,13))
C     Amplitude(s) for diagram number 1
      CALL L1_0(W(1,8),W(1,3),W(1,9),GQQ,AMP(1))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,7),W(1,10),W(1,9),GQQ,AMP(2))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,11),W(1,10),W(1,6),GQQ,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,11),W(1,3),W(1,12),GQQ,AMP(4))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,13),W(1,3),W(1,9),GQQ,AMP(5))
      JAMP(1)=+1./2.*(-AMP(1)-AMP(2)-AMP(3)-IMAG1*AMP(4)-IMAG1*AMP(5))

      FLOW4=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info(),
     """      COMPLEX*16 FUNCTION FLOW5(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 5
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=2 QED=0 WEIGHTED=4 singlet_QCD=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=4)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=12)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_2(W(1,1),W(1,6),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_3(W(1,7),W(1,3),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,2),W(1,5),GQQ,ZERO, ZERO, W(1,9))
      CALL L1_2(W(1,4),W(1,5),GQQ,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,3),W(1,6),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_3(W(1,1),W(1,11),GQQ,ZERO, ZERO, W(1,12))
C     Amplitude(s) for diagram number 1
      CALL L1_0(W(1,4),W(1,9),W(1,8),GQQ,AMP(1))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,10),W(1,2),W(1,8),GQQ,AMP(2))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,4),W(1,9),W(1,12),GQQ,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,10),W(1,2),W(1,12),GQQ,AMP(4))
      JAMP(1)=+1./2.*(-AMP(1)-AMP(2)-AMP(3)-AMP(4))

      FLOW5=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info(),
     """      COMPLEX*16 FUNCTION FLOW6(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 6
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=2 QED=0 WEIGHTED=4 singlet_QCD=2
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=5)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=12)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/1,-1,1,-1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL OXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL OXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL IXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_3(W(1,1),W(1,3),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_2(W(1,4),W(1,7),GQQ,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,2),W(1,5),GQQ,ZERO, ZERO, W(1,9))
      CALL L1_2(W(1,4),W(1,6),GQQ,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,2),W(1,7),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,12))
C     Amplitude(s) for diagram number 1
      CALL L1_0(W(1,8),W(1,9),W(1,6),GQQ,AMP(1))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,10),W(1,9),W(1,7),GQQ,AMP(2))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,10),W(1,11),W(1,5),GQQ,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,4),W(1,11),W(1,12),GQQ,AMP(4))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,8),W(1,2),W(1,12),GQQ,AMP(5))
      JAMP(1)=+1./2.*(-AMP(1)-AMP(2)-AMP(3)-IMAG1*AMP(4)-IMAG1*AMP(5))

      FLOW6=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info()]
     

        for iflow, color_flow in enumerate(matrix_element.get('color_flows')):

            process_exporter.write_co_flow_v4(\
                writers.FortranWriter(self.give_pos('test')),
                color_flow,
                self.myfortranmodel)

            #print open(self.give_pos('test')).read()
            self.assertFileContains('test', goal_flow_f[iflow])

        goal_matrix_f = \
"""C     -------------------------
      SUBROUTINE SMATRIX(P,ANS)
C     -------------------------
C     
C     Generated by MadGraph 5 v. 1.4.5, 2012-04-20
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX
      INTEGER IHEL,IDEN, I
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/72/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
          T=MATRIX(P ,NHEL(1,IHEL))
          ANS=ANS+T
          IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
            GOODHEL(IHEL)=.TRUE.
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END

C     ------------------------------
      REAL*8 FUNCTION MATRIX(P,NHEL)
C     ------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external momenta P(0:3,NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NPERMS
      PARAMETER (NPERMS=4)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,PERM(NEXTERNAL)
      COMPLEX*16 ZTEMP
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 ONEPERM
      EXTERNAL ONEPERM

      ZTEMP = (0.D0,0.D0)
      DO I=1,NPERMS
        CALL GETPERM(I,PERM)
        ZTEMP=ZTEMP+ONEPERM(P,NHEL,PERM)
      ENDDO
      MATRIX=REAL(ZTEMP)

      RETURN
      END

C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM(P,NHEL,PM)
C     --------------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NFLOWS
      PARAMETER (NFLOWS=19)
      INTEGER    NPERMS
      PARAMETER (NPERMS=4)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NFLOWS)
      INTEGER PERMS(NEXTERNAL,NPERMS),IFERM(NPERMS),PERM(NEXTERNAL)
      DATA (PERMS(I,   1),I=1,6) / 1, 2, 3, 4, 5, 6/
      DATA (PERMS(I,   2),I=1,6) / 1, 2, 3, 4, 6, 5/
      DATA (PERMS(I,   3),I=1,6) / 1, 3, 2, 4, 5, 6/
      DATA (PERMS(I,   4),I=1,6) / 1, 3, 2, 4, 6, 5/
      DATA IFERM/ 1, 1,-1,-1/
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 FLOW1,FLOW2,FLOW3,FLOW4,FLOW5,FLOW6
      EXTERNAL FLOW1,FLOW2,FLOW3,FLOW4,FLOW5,FLOW6
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,1))
      ENDDO
      JAMP(1)=IFERM(1)*FLOW1(P,NHEL,PERM)
      JAMP(2)=IFERM(1)*FLOW3(P,NHEL,PERM)
      JAMP(3)=IFERM(1)*FLOW4(P,NHEL,PERM)
      JAMP(4)=IFERM(1)*FLOW2(P,NHEL,PERM)
      JAMP(5)=IFERM(1)*FLOW5(P,NHEL,PERM)
      JAMP(6)=IFERM(1)*FLOW6(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,2))
      ENDDO
      JAMP(7)=IFERM(2)*FLOW1(P,NHEL,PERM)
      JAMP(8)=IFERM(2)*FLOW2(P,NHEL,PERM)
      JAMP(9)=IFERM(2)*FLOW3(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,3))
      ENDDO
      JAMP(10)=IFERM(3)*FLOW1(P,NHEL,PERM)
      JAMP(11)=IFERM(3)*FLOW2(P,NHEL,PERM)
      JAMP(12)=IFERM(3)*FLOW3(P,NHEL,PERM)
      JAMP(13)=IFERM(3)*FLOW4(P,NHEL,PERM)
      JAMP(14)=IFERM(3)*FLOW5(P,NHEL,PERM)
      JAMP(15)=IFERM(3)*FLOW6(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,4))
      ENDDO
      JAMP(16)=IFERM(4)*FLOW2(P,NHEL,PERM)
      JAMP(17)=IFERM(4)*FLOW3(P,NHEL,PERM)
      JAMP(18)=IFERM(4)*FLOW5(P,NHEL,PERM)
      JAMP(19)=IFERM(4)*FLOW6(P,NHEL,PERM)
      ZTEMP = (0.D0,0.D0)
      ZTEMP = ZTEMP+1/9D0*JAMP(1)*DCONJG(144D0*(JAMP(1))+48D0*(JAMP(3)
     $ +JAMP(10)+JAMP(11)+JAMP(12))+18D0*(JAMP(2)-JAMP(7)+JAMP(9))
     $ +16D0*(JAMP(13)+JAMP(14)+JAMP(15)))
      ZTEMP = ZTEMP+1/9D0*JAMP(4)*DCONJG(144D0*(JAMP(4))+48D0*(JAMP(5)
     $ +JAMP(10)+JAMP(11)+JAMP(16)+JAMP(17))+18D0*(JAMP(8))+16D0
     $ *(JAMP(13)+JAMP(14)+JAMP(18)+JAMP(19)))
      ZTEMP = ZTEMP+1/9D0*JAMP(2)*DCONJG(144D0*(JAMP(2))+48D0*(JAMP(6)
     $ +JAMP(10)+JAMP(12)+JAMP(16))+18D0*(JAMP(1)+JAMP(7)-JAMP(9))
     $ +16D0*(JAMP(13)+JAMP(15)+JAMP(18)))
      ZTEMP = ZTEMP+1/9D0*JAMP(3)*DCONJG(48D0*(JAMP(1))+16D0*(JAMP(3)
     $ +JAMP(10)+JAMP(11)+JAMP(12)))
      ZTEMP = ZTEMP+1/9D0*JAMP(5)*DCONJG(48D0*(JAMP(4))+16D0*(JAMP(5)
     $ +JAMP(10)+JAMP(11)+JAMP(16)+JAMP(17)))
      ZTEMP = ZTEMP+1/9D0*JAMP(6)*DCONJG(48D0*(JAMP(2))+16D0*(JAMP(6)
     $ +JAMP(10)+JAMP(12)+JAMP(16)))
      ONEPERM=ZTEMP

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM(IPERM,PERM)
C     ------------------------------------
C     
C     Gives permutation number IPERM. 
C     Return value is the fermion factor due to PERM
C     
C     Process: u u~ > u u~ g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
C     
C     ARGUMENTS 
C     
      INTEGER IPERM,PERM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IFLAG
      LOGICAL OK
      INTEGER COMP(NEXTERNAL)
      DATA COMP/1,2,2,3,4,4/
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=I
      ENDDO
      I=1
      DO WHILE(I.LT.IPERM)
        CALL IPNEXT(PERM,NEXTERNAL,IFLAG)
        OK=.TRUE.
        DO J=1,NEXTERNAL
          IF(COMP(PERM(J)).NE.COMP(J))THEN
            OK=.FALSE.
            EXIT
          ENDIF
        ENDDO
        IF(OK) I=I+1
      ENDDO
      END

""" % misc.get_pkg_info()

        process_exporter.write_matrix_element_v4(\
            writers.FortranWriter(self.give_pos('test')),
            matrix_element)

        #print open(self.give_pos('test')).read()
        self.assertFileContains('test', goal_matrix_f)

    def test_export_co_matrix_element_v4_standalone_gg_uxugg(self):
        """Test the result of exporting g g > u u~ g g to standalone"""

        # Create the amplitude
        myleglist = base_objects.LegList([\
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':-2, 'state':True}),
            base_objects.Leg({'id':2, 'state':True}),
            base_objects.Leg({'id':21, 'state':True}),
            base_objects.Leg({'id':21, 'state':True})])

        myproc = base_objects.Process({'legs':myleglist,
                                       'orders':{'QED': 0},
                                       'model':self.mymodel})

        self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

        matrix_element = color_ordered_helas_objects.COHelasMatrixElement(\
            self.myamplitude, gen_color=3, optimization=3)

        process_exporter = color_ordered_export_v4.ProcessExporterFortranCOSA()

        goal_flow_f = \
    ["""      COMPLEX*16 FUNCTION FLOW1(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 1
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: g g > u~ u g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=16)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=23)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/-1,-1,-1,1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL VXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL VXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL IXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL OXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_3(W(1,3),W(1,4),GQQ,ZERO, ZERO, W(1,7))
      CALL L1_1(W(1,2),W(1,1),G,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,9))
      CALL L1_2(W(1,3),W(1,6),GQQ,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,4),W(1,1),GQQ,ZERO, ZERO, W(1,11))
      CALL L1_1(W(1,5),W(1,2),G,ZERO, ZERO, W(1,12))
      CALL L1_1(W(1,4),W(1,8),GQQ,ZERO, ZERO, W(1,13))
      CALL L1_1(W(1,11),W(1,2),GQQ,ZERO, ZERO, W(1,14))
      CALL SUMF2(1.*IMAG1,W(1,13),1.*ONE,W(1,14),W(1,15))
      CALL L1_3(W(1,3),W(1,11),GQQ,ZERO, ZERO, W(1,16))
      CALL L1_1(W(1,7),W(1,1),G,ZERO, ZERO, W(1,17))
      CALL SUMV2(1.*ONE,W(1,16),-1.*IMAG1,W(1,17),W(1,18))
      CALL L1_1(W(1,5),W(1,8),G,ZERO, ZERO, W(1,19))
      CALL L1_1(W(1,12),W(1,1),G,ZERO, ZERO, W(1,20))
      CALL VVVV4_1(W(1,5),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,21))
      CALL VVVV1_1(W(1,5),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,22))
      CALL SUMV4(-2.*ONE,W(1,19),-2.*ONE,W(1,20),2.*ONE,W(1,21),
     $ -2.*ONE,W(1,22),W(1,23))
C     Amplitude(s) for diagram number 1
      CALL VVVV3_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(1))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(2))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,6),W(1,7),W(1,23),G,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,9),W(1,7),W(1,8),G,AMP(4))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,10),W(1,15),W(1,5),GQQ,AMP(5))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,10),W(1,4),W(1,23),GQQ,AMP(6))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,3),W(1,15),W(1,9),GQQ,AMP(7))
C     Amplitude(s) for diagram number 7
      CALL L1_0(W(1,6),W(1,12),W(1,18),G,AMP(8))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,10),W(1,11),W(1,12),GQQ,AMP(9))
C     Amplitude(s) for diagram number 9
      CALL L1_0(W(1,9),W(1,2),W(1,18),G,AMP(10))
C     Amplitude(s) for diagram number 10
      CALL VVVV4_0(W(1,6),W(1,5),W(1,2),W(1,18),GG,AMP(11))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,2),W(1,18),GG,AMP(12))
C     Amplitude(s) for diagram number 11
      CALL VVVV3_0(W(1,6),W(1,7),W(1,12),W(1,1),GG,AMP(13))
      CALL VVVV1_0(W(1,6),W(1,7),W(1,12),W(1,1),GG,AMP(14))
C     Amplitude(s) for diagram number 12
      CALL VVVV3_0(W(1,9),W(1,7),W(1,2),W(1,1),GG,AMP(15))
      CALL VVVV1_0(W(1,9),W(1,7),W(1,2),W(1,1),GG,AMP(16))
      JAMP(1)=+IMAG1*AMP(1)+IMAG1*AMP(2)-1./2.*IMAG1*AMP(3)+IMAG1
     $ *AMP(4)+AMP(5)+1./2.*AMP(6)+IMAG1*AMP(7)-AMP(8)+IMAG1*AMP(9)
     $ -AMP(10)+AMP(11)-AMP(12)+IMAG1*AMP(13)+IMAG1*AMP(14)+IMAG1
     $ *AMP(15)+IMAG1*AMP(16)

      FLOW1=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info()]
     

        for iflow, color_flow in enumerate(matrix_element.get('color_flows')):

            process_exporter.write_co_flow_v4(\
                writers.FortranWriter(self.give_pos('test')),
                color_flow,
                self.myfortranmodel)

            #print open(self.give_pos('test')).read()
            self.assertFileContains('test', goal_flow_f[iflow])

        goal_matrix_f = \
"""C     -------------------------
      SUBROUTINE SMATRIX(P,ANS)
C     -------------------------
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: g g > u~ u g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX
      INTEGER IHEL,IDEN, I
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/512/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
          T=MATRIX(P ,NHEL(1,IHEL))
          ANS=ANS+T
          IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
            GOODHEL(IHEL)=.TRUE.
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END

C     ------------------------------
      REAL*8 FUNCTION MATRIX(P,NHEL)
C     ------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external momenta P(0:3,NEXTERNAL)
C     
C     Process: g g > u~ u g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NPERMS
      PARAMETER (NPERMS=24)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,PERM(NEXTERNAL)
      COMPLEX*16 ZTEMP
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 ONEPERM
      EXTERNAL ONEPERM

      ZTEMP = (0.D0,0.D0)
      DO I=1,NPERMS
        CALL GETPERM(I,PERM)
        ZTEMP=ZTEMP+ONEPERM(P,NHEL,PERM)
      ENDDO
      MATRIX=REAL(ZTEMP)

      RETURN
      END

C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM(P,NHEL,PM)
C     --------------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g g > u~ u g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NFLOWS
      PARAMETER (NFLOWS=10)
      INTEGER    NPERMS
      PARAMETER (NPERMS=10)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NFLOWS)
      INTEGER PERMS(NEXTERNAL,NPERMS),IFERM(NPERMS),PERM(NEXTERNAL)
      DATA (PERMS(I,   1),I=1,6) / 1, 2, 3, 4, 5, 6/
      DATA (PERMS(I,   2),I=1,6) / 1, 2, 3, 4, 6, 5/
      DATA (PERMS(I,   3),I=1,6) / 1, 5, 3, 4, 2, 6/
      DATA (PERMS(I,   4),I=1,6) / 1, 6, 3, 4, 5, 2/
      DATA (PERMS(I,   5),I=1,6) / 2, 1, 3, 4, 5, 6/
      DATA (PERMS(I,   6),I=1,6) / 5, 2, 3, 4, 1, 6/
      DATA (PERMS(I,   7),I=1,6) / 5, 6, 3, 4, 1, 2/
      DATA (PERMS(I,   8),I=1,6) / 5, 6, 3, 4, 2, 1/
      DATA (PERMS(I,   9),I=1,6) / 6, 2, 3, 4, 5, 1/
      DATA (PERMS(I,  10),I=1,6) / 6, 5, 3, 4, 1, 2/
      DATA IFERM/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 FLOW1
      EXTERNAL FLOW1
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,1))
      ENDDO
      JAMP(1)=IFERM(1)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,2))
      ENDDO
      JAMP(2)=IFERM(2)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,3))
      ENDDO
      JAMP(3)=IFERM(3)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,4))
      ENDDO
      JAMP(4)=IFERM(4)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,5))
      ENDDO
      JAMP(5)=IFERM(5)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,6))
      ENDDO
      JAMP(6)=IFERM(6)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,7))
      ENDDO
      JAMP(7)=IFERM(7)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,8))
      ENDDO
      JAMP(8)=IFERM(8)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,9))
      ENDDO
      JAMP(9)=IFERM(9)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,10))
      ENDDO
      JAMP(10)=IFERM(10)*FLOW1(P,NHEL,PERM)
      ZTEMP = (0.D0,0.D0)
      ZTEMP = ZTEMP+1/54D0*JAMP(1)*DCONJG(512D0*(JAMP(1))+80D0
     $ *(JAMP(4)+JAMP(6))+71D0*(JAMP(7))+64D0*(-JAMP(2)-JAMP(3)
     $ -JAMP(5))+62D0*(JAMP(8)+JAMP(9)+JAMP(10)))
      ONEPERM=ZTEMP

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM(IPERM,PERM)
C     ------------------------------------
C     
C     Gives permutation number IPERM. 
C     Return value is the fermion factor due to PERM
C     
C     Process: g g > u~ u g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
C     
C     ARGUMENTS 
C     
      INTEGER IPERM,PERM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IFLAG
      LOGICAL OK
      INTEGER COMP(NEXTERNAL)
      DATA COMP/1,1,2,3,1,1/
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=I
      ENDDO
      I=1
      DO WHILE(I.LT.IPERM)
        CALL IPNEXT(PERM,NEXTERNAL,IFLAG)
        OK=.TRUE.
        DO J=1,NEXTERNAL
          IF(COMP(PERM(J)).NE.COMP(J))THEN
            OK=.FALSE.
            EXIT
          ENDIF
        ENDDO
        IF(OK) I=I+1
      ENDDO
      END

""" % misc.get_pkg_info()

        process_exporter.write_matrix_element_v4(\
            writers.FortranWriter(self.give_pos('test')),
            matrix_element)

        #print open(self.give_pos('test')).read()
        self.assertFileContains('test', goal_matrix_f)

    def test_export_co_matrix_element_v4_madevent_gg_4g(self):
        """Test the result of exporting g g > g g g g to MadEvent"""

        # Create the amplitude
        myleglist = base_objects.LegList([\
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':21, 'state':True}),
            base_objects.Leg({'id':21, 'state':True}),
            base_objects.Leg({'id':21, 'state':True}),
            base_objects.Leg({'id':21, 'state':True})])

        myproc = base_objects.Process({'legs':myleglist,
                                       'orders':{'QED': 0},
                                       'model':self.mymodel})

        self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

        matrix_element = color_ordered_helas_objects.COHelasMatrixElement(\
            self.myamplitude, gen_color=3, optimization=3,
            gen_periferal_diagrams = True)

        process_exporter = color_ordered_export_v4.ProcessExporterFortranCOME()

        goal_flow_f = \
    ["""      COMPLEX*16 FUNCTION FLOW1(P,NHEL,IP)
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns color ordered amplitude for color flow 1
C     for the point with external momenta P(0:3, NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=32)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=27)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1, ONE
      PARAMETER (IMAG1=(0D0,1D0),ONE=(1D0,0D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),IP(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL)
      COMPLEX*16 AMP(NGRAPHS), JAMP(1)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      DATA IC/-1,-1,1,1,1,1/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      CALL VXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL VXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL VXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL VXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_1(W(1,4),W(1,3),G,ZERO, ZERO, W(1,7))
      CALL L1_1(W(1,2),W(1,1),G,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,9))
      CALL L1_1(W(1,5),W(1,4),G,ZERO, ZERO, W(1,10))
      CALL L1_1(W(1,3),W(1,2),G,ZERO, ZERO, W(1,11))
      CALL L1_1(W(1,6),W(1,1),G,ZERO, ZERO, W(1,12))
      CALL L1_1(W(1,6),W(1,8),G,ZERO, ZERO, W(1,13))
      CALL L1_1(W(1,2),W(1,12),G,ZERO, ZERO, W(1,14))
      CALL VVVV3_1(W(1,6),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,15))
      CALL VVVV1_1(W(1,6),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,16))
      CALL SUMV4(2.*ONE,W(1,13),2.*ONE,W(1,14),2.*ONE,W(1,15),2.
     $ *ONE,W(1,16),W(1,17))
      CALL L1_1(W(1,3),W(1,8),G,ZERO, ZERO, W(1,18))
      CALL L1_1(W(1,11),W(1,1),G,ZERO, ZERO, W(1,19))
      CALL VVVV4_1(W(1,3),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,20))
      CALL VVVV1_1(W(1,3),W(1,2),W(1,1),GG,ZERO, ZERO, W(1,21))
      CALL SUMV4(-2.*ONE,W(1,18),-2.*ONE,W(1,19),2.*ONE,W(1,20),
     $ -2.*ONE,W(1,21),W(1,22))
      CALL L1_1(W(1,5),W(1,12),G,ZERO, ZERO, W(1,23))
      CALL L1_1(W(1,9),W(1,1),G,ZERO, ZERO, W(1,24))
      CALL VVVV4_1(W(1,6),W(1,5),W(1,1),GG,ZERO, ZERO, W(1,25))
      CALL VVVV3_1(W(1,6),W(1,5),W(1,1),GG,ZERO, ZERO, W(1,26))
      CALL SUMV4(-2.*ONE,W(1,23),2.*ONE,W(1,24),-2.*ONE,W(1,25),
     $ -2.*ONE,W(1,26),W(1,27))
C     Amplitude(s) for diagram number 1
      CALL VVVV4_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(1))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,7),W(1,8),GG,AMP(2))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,5),W(1,7),W(1,17),G,AMP(3))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,9),W(1,7),W(1,8),G,AMP(4))
C     Amplitude(s) for diagram number 4
      CALL VVVV4_0(W(1,6),W(1,10),W(1,3),W(1,8),GG,AMP(5))
      CALL VVVV1_0(W(1,6),W(1,10),W(1,3),W(1,8),GG,AMP(6))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,6),W(1,10),W(1,22),G,AMP(7))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,10),W(1,3),W(1,17),G,AMP(8))
C     Amplitude(s) for diagram number 7
      CALL VVVV4_0(W(1,9),W(1,4),W(1,3),W(1,8),GG,AMP(9))
      CALL VVVV1_0(W(1,9),W(1,4),W(1,3),W(1,8),GG,AMP(10))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,9),W(1,4),W(1,22),G,AMP(11))
C     Amplitude(s) for diagram number 9
      CALL VVVV4_0(W(1,5),W(1,4),W(1,3),W(1,17),GG,AMP(12))
      CALL VVVV1_0(W(1,5),W(1,4),W(1,3),W(1,17),GG,AMP(13))
C     Amplitude(s) for diagram number 10
      CALL VVVV4_0(W(1,6),W(1,5),W(1,4),W(1,22),GG,AMP(14))
      CALL VVVV1_0(W(1,6),W(1,5),W(1,4),W(1,22),GG,AMP(15))
C     Amplitude(s) for diagram number 11
      CALL VVVV4_0(W(1,5),W(1,4),W(1,11),W(1,12),GG,AMP(16))
      CALL VVVV1_0(W(1,5),W(1,4),W(1,11),W(1,12),GG,AMP(17))
C     Amplitude(s) for diagram number 12
      CALL L1_0(W(1,4),W(1,11),W(1,27),G,AMP(18))
C     Amplitude(s) for diagram number 13
      CALL L1_0(W(1,10),W(1,11),W(1,12),G,AMP(19))
C     Amplitude(s) for diagram number 14
      CALL VVVV4_0(W(1,5),W(1,7),W(1,2),W(1,12),GG,AMP(20))
      CALL VVVV1_0(W(1,5),W(1,7),W(1,2),W(1,12),GG,AMP(21))
C     Amplitude(s) for diagram number 15
      CALL L1_0(W(1,7),W(1,2),W(1,27),G,AMP(22))
C     Amplitude(s) for diagram number 16
      CALL VVVV4_0(W(1,10),W(1,3),W(1,2),W(1,12),GG,AMP(23))
      CALL VVVV1_0(W(1,10),W(1,3),W(1,2),W(1,12),GG,AMP(24))
C     Amplitude(s) for diagram number 17
      CALL VVVV4_0(W(1,4),W(1,3),W(1,2),W(1,27),GG,AMP(25))
      CALL VVVV1_0(W(1,4),W(1,3),W(1,2),W(1,27),GG,AMP(26))
C     Amplitude(s) for diagram number 18
      CALL VVVV4_0(W(1,6),W(1,10),W(1,11),W(1,1),GG,AMP(27))
      CALL VVVV1_0(W(1,6),W(1,10),W(1,11),W(1,1),GG,AMP(28))
C     Amplitude(s) for diagram number 19
      CALL VVVV4_0(W(1,9),W(1,4),W(1,11),W(1,1),GG,AMP(29))
      CALL VVVV1_0(W(1,9),W(1,4),W(1,11),W(1,1),GG,AMP(30))
C     Amplitude(s) for diagram number 20
      CALL VVVV4_0(W(1,9),W(1,7),W(1,2),W(1,1),GG,AMP(31))
      CALL VVVV1_0(W(1,9),W(1,7),W(1,2),W(1,1),GG,AMP(32))
      JAMP(1)=-2*AMP(1)+2*AMP(2)-AMP(3)+2*AMP(4)-2*AMP(5)+2*AMP(6)
     $ -AMP(7)-AMP(8)-2*AMP(9)+2*AMP(10)-AMP(11)+AMP(12)-AMP(13)
     $ +AMP(14)-AMP(15)+2*AMP(16)-2*AMP(17)-AMP(18)-2*AMP(19)
     $ +2*AMP(20)-2*AMP(21)-AMP(22)+2*AMP(23)-2*AMP(24)+AMP(25)
     $ -AMP(26)-2*AMP(27)+2*AMP(28)-2*AMP(29)+2*AMP(30)-2*AMP(31)
     $ +2*AMP(32)

      FLOW1=JAMP(1)
      RETURN
      END

""" % misc.get_pkg_info()]
     

        for iflow, color_flow in enumerate(matrix_element.get('color_flows')):

            process_exporter.write_co_flow_v4(\
                writers.FortranWriter(self.give_pos('test')),
                color_flow,
                self.myfortranmodel)

            #print open(self.give_pos('test')).read()
            self.assertFileContains('test', goal_flow_f[iflow])

        goal_matrix_f = \
"""C     -------------------------
      SUBROUTINE SMATRIX(P,ANS)
C     -------------------------
C     
C     Generated by MadGraph 5 v. %(version)s, %(date)s
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Color-ordered MadGraph for Madevent Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'genps.inc'
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=51)
      INTEGER    NDIAGS
      PARAMETER (NDIAGS=51)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB)
      INTEGER    NPERMS
      PARAMETER (NPERMS=120)
      INTEGER    NFLOWS
      PARAMETER (NFLOWS=1)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, MATRIX
      REAL*8 R,SUMHEL,TS(NCOMB)
      INTEGER I,IDEN
      INTEGER IPROC,II
      LOGICAL GOODHEL(NCOMB)
      REAL*8 HWGT,XTOT,XTRY,XREJ,XR,YFRAC(0:NCOMB),FRAC
      INTEGER IDUM,IPERM,NGOOD,IGOOD(NCOMB),JHEL,J,JJ
      REAL*8 TOTFACT,FACTNUL
      INTEGER JHELBK,NUPPER
      LOGICAL FIRSTFLOW
      REAL*8   RVEC
      SAVE FRAC,FACTNUL
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS),JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,         JAMP2

      CHARACTER*101        HEL_BUFF
      COMMON/TO_HELICITY/  HEL_BUFF

      REAL*8 POL(2)
      COMMON/TO_POLARIZATION/ POL

      INTEGER          ISUM_HEL
      LOGICAL                    MULTI_CHANNEL
      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL

      INTEGER PERM(NEXTERNAL)
      COMMON/TO_COPERM/PERM

      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
      DATA NTRY,IDUM /0,-1/
      DATA XTRY, XREJ, NGOOD /0,0,0/
      SAVE YFRAC, IGOOD, JHEL
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/6144/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1

      IF (MULTI_CHANNEL) THEN
        DO I=1,NDIAGS
          AMP2(I)=0D0
        ENDDO
        JAMP2(0)=NFLOWS
        DO I=1,INT(JAMP2(0))
          JAMP2(I)=0D0
        ENDDO
      ENDIF
      ANS = 0D0
      WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
      DO I=1,NCOMB
        TS(I)=0D0
      ENDDO
      IF(NTRY.EQ.1)THEN
C       Calculate summed permutation weight
        TOTFACT=NPERMS
        FRAC=MIN(NUMPERMS,NPERMS)*1D0/TOTFACT
        FACTNUL=1D0
C       Correction factor for events where no permutation is selected
        DO IPERM=1,NPERMS
          FACTNUL=FACTNUL*(1-FRAC)
        ENDDO
        FACTNUL=1-FACTNUL
      ENDIF
C     Decide between helicity all-sum case (ISUM_HEL=0) and partial
C      sum case
      JHELBK=JHEL
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LE. MAXTRIES) THEN
        NUPPER=NCOMB
        HWGT=1.
      ELSE
        NUPPER=ISUM_HEL
        HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
      ENDIF
 10   FIRSTFLOW=.TRUE.
C     Start loop over permutations
      DO IPERM=1,NPERMS
        JHEL=JHELBK
C       Decide whether to use this permutation
        CALL RANMAR(RVEC)
        IF(RVEC.GT.FRAC) CYCLE
C       Get permutation
        CALL GETPERM(IPERM,PERM)
C       Start loop over helicities
        DO J=1,NUPPER
          IF (ISUM_HEL .EQ. 0 .OR. NTRY .LE. MAXTRIES) THEN
            I=J
            IF(.NOT.GOODHEL(I) .AND. NTRY .GT. MAXTRIES) CYCLE
          ELSE
            JHEL=JHEL+1
            IF (JHEL .GT. NGOOD) JHEL=1
            I = IGOOD(JHEL)
          ENDIF
          T=MATRIX(P,NHEL(1,I),PERM,FIRSTFLOW)
          DO JJ=1,NINCOMING
            IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ)))
     $       ) THEN
              T=T*ABS(POL(JJ))
            ELSE IF(POL(JJ).NE.1D0)THEN
              T=T*(2D0-ABS(POL(JJ)))
            ENDIF
          ENDDO
          ANS=ANS+T*HWGT/FRAC*FACTNUL
          TS(I)=T*HWGT/FRAC*FACTNUL
        ENDDO
        FIRSTFLOW=.FALSE.
C       print *,'ANS: ',ANS
        IF(NTRY.LE.MAXTRIES)THEN
          DO I=1,NCOMB
            IF (.NOT.GOODHEL(I) .AND. (TS(I).GT.ANS*LIMHEL/NCOMB)) THEN
              GOODHEL(I)=.TRUE.
              NGOOD = NGOOD +1
              IGOOD(NGOOD) = I
              PRINT *,'Added good helicity ',I,TS(I)/ANS,' in event '
     $         ,NTRY(IMIRROR)
            ENDIF
          ENDDO
        ENDIF
        IF(NTRY.EQ.MAXTRIES)THEN
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
        ENDIF
        IF (ISUM_HEL .EQ. 1) THEN
          WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
        ENDIF
      ENDDO
C     If no permutation was chosen, repeat
      IF(FIRSTFLOW) GOTO 10
      IF (ISUM_HEL .NE. 1) THEN
        CALL RANMAR(RVEC)
        R=RVEC*ANS
        SUMHEL=0D0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+TS(I)
          IF(R.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
            GOTO 20
          ENDIF
        ENDDO
 20     CONTINUE
      ENDIF
      IF (MULTI_CHANNEL) THEN
        XTOT=0D0
        DO I=1,NDIAGS
          XTOT=XTOT+AMP2(I)
        ENDDO
C       print *,'XTOT,AMP2: ',XTOT,AMP2(SUBDIAG(1))
        IF (XTOT.NE.0D0) THEN
          ANS=ANS*AMP2(MAPCONFIG(ICONFIG))/XTOT
        ELSE
          ANS=0D0
        ENDIF
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END

C     --------------------------------------------------------
      REAL*8 FUNCTION MATRIX(P,NHEL,PERM,FIRSTFLOW)
C     --------------------------------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external momenta P(0:3,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=51)
      INCLUDE 'genps.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INCLUDE 'coupl.inc'
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=38)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL),PERM(NEXTERNAL)
      LOGICAL FIRSTFLOW
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IC(NEXTERNAL),IP(NEXTERNAL)
      COMPLEX*16 ZTEMP
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 W(18,NWAVEFUNCS)
      DATA IC/-1,-1,1,1,1,1/
      DATA IP/1,2,3,4,5,6/
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 ONEPERM
      EXTERNAL ONEPERM

C     Skip calculating amp2s for PS if secondary flows
      IF(.NOT.FIRSTFLOW) GOTO 10

      CALL VXXXXX(P(0,IP(1)),ZERO,NHEL(IP(1)),1*IC(IP(1)),W(1,1))
      CALL VXXXXX(P(0,IP(2)),ZERO,NHEL(IP(2)),1*IC(IP(2)),W(1,2))
      CALL VXXXXX(P(0,IP(3)),ZERO,NHEL(IP(3)),1*IC(IP(3)),W(1,3))
      CALL VXXXXX(P(0,IP(4)),ZERO,NHEL(IP(4)),1*IC(IP(4)),W(1,4))
      CALL VXXXXX(P(0,IP(5)),ZERO,NHEL(IP(5)),1*IC(IP(5)),W(1,5))
      CALL VXXXXX(P(0,IP(6)),ZERO,NHEL(IP(6)),1*IC(IP(6)),W(1,6))
      CALL L1_1(W(1,6),W(1,1),G,ZERO, ZERO, W(1,7))
      CALL L1_1(W(1,5),W(1,7),G,ZERO, ZERO, W(1,8))
      CALL L1_1(W(1,3),W(1,2),G,ZERO, ZERO, W(1,9))
C     Amplitude(s) for diagram number 1
      CALL L1_0(W(1,4),W(1,9),W(1,8),G,AMP(1))
      CALL L1_1(W(1,5),W(1,4),G,ZERO, ZERO, W(1,10))
C     Amplitude(s) for diagram number 2
      CALL L1_0(W(1,10),W(1,9),W(1,7),G,AMP(2))
      CALL L1_1(W(1,4),W(1,3),G,ZERO, ZERO, W(1,11))
C     Amplitude(s) for diagram number 3
      CALL L1_0(W(1,11),W(1,2),W(1,8),G,AMP(3))
      CALL L1_1(W(1,6),W(1,5),G,ZERO, ZERO, W(1,12))
      CALL L1_1(W(1,12),W(1,1),G,ZERO, ZERO, W(1,13))
C     Amplitude(s) for diagram number 4
      CALL L1_0(W(1,4),W(1,9),W(1,13),G,AMP(4))
C     Amplitude(s) for diagram number 5
      CALL L1_0(W(1,11),W(1,2),W(1,13),G,AMP(5))
      CALL L1_1(W(1,5),W(1,1),G,ZERO, ZERO, W(1,14))
      CALL L1_1(W(1,6),W(1,14),G,ZERO, ZERO, W(1,15))
C     Amplitude(s) for diagram number 6
      CALL L1_0(W(1,4),W(1,9),W(1,15),G,AMP(6))
      CALL L1_1(W(1,6),W(1,4),G,ZERO, ZERO, W(1,16))
C     Amplitude(s) for diagram number 7
      CALL L1_0(W(1,16),W(1,9),W(1,14),G,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL L1_0(W(1,11),W(1,2),W(1,15),G,AMP(8))
      CALL L1_1(W(1,4),W(1,7),G,ZERO, ZERO, W(1,17))
      CALL L1_1(W(1,5),W(1,3),G,ZERO, ZERO, W(1,18))
C     Amplitude(s) for diagram number 9
      CALL L1_0(W(1,18),W(1,2),W(1,17),G,AMP(9))
      CALL L1_1(W(1,16),W(1,1),G,ZERO, ZERO, W(1,19))
C     Amplitude(s) for diagram number 10
      CALL L1_0(W(1,5),W(1,9),W(1,19),G,AMP(10))
C     Amplitude(s) for diagram number 11
      CALL L1_0(W(1,18),W(1,2),W(1,19),G,AMP(11))
      CALL L1_1(W(1,4),W(1,1),G,ZERO, ZERO, W(1,20))
      CALL L1_1(W(1,6),W(1,20),G,ZERO, ZERO, W(1,21))
C     Amplitude(s) for diagram number 12
      CALL L1_0(W(1,5),W(1,9),W(1,21),G,AMP(12))
C     Amplitude(s) for diagram number 13
      CALL L1_0(W(1,12),W(1,9),W(1,20),G,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL L1_0(W(1,18),W(1,2),W(1,21),G,AMP(14))
      CALL L1_1(W(1,4),W(1,14),G,ZERO, ZERO, W(1,22))
      CALL L1_1(W(1,6),W(1,3),G,ZERO, ZERO, W(1,23))
C     Amplitude(s) for diagram number 15
      CALL L1_0(W(1,23),W(1,2),W(1,22),G,AMP(15))
      CALL L1_1(W(1,10),W(1,1),G,ZERO, ZERO, W(1,24))
C     Amplitude(s) for diagram number 16
      CALL L1_0(W(1,6),W(1,9),W(1,24),G,AMP(16))
C     Amplitude(s) for diagram number 17
      CALL L1_0(W(1,23),W(1,2),W(1,24),G,AMP(17))
      CALL L1_1(W(1,5),W(1,20),G,ZERO, ZERO, W(1,25))
C     Amplitude(s) for diagram number 18
      CALL L1_0(W(1,23),W(1,2),W(1,25),G,AMP(18))
      CALL L1_1(W(1,4),W(1,2),G,ZERO, ZERO, W(1,26))
C     Amplitude(s) for diagram number 19
      CALL L1_0(W(1,3),W(1,26),W(1,8),G,AMP(19))
C     Amplitude(s) for diagram number 20
      CALL L1_0(W(1,18),W(1,26),W(1,7),G,AMP(20))
C     Amplitude(s) for diagram number 21
      CALL L1_0(W(1,3),W(1,26),W(1,13),G,AMP(21))
C     Amplitude(s) for diagram number 22
      CALL L1_0(W(1,3),W(1,26),W(1,15),G,AMP(22))
C     Amplitude(s) for diagram number 23
      CALL L1_0(W(1,23),W(1,26),W(1,14),G,AMP(23))
      CALL L1_1(W(1,3),W(1,7),G,ZERO, ZERO, W(1,27))
C     Amplitude(s) for diagram number 24
      CALL L1_0(W(1,10),W(1,2),W(1,27),G,AMP(24))
      CALL L1_1(W(1,23),W(1,1),G,ZERO, ZERO, W(1,28))
C     Amplitude(s) for diagram number 25
      CALL L1_0(W(1,5),W(1,26),W(1,28),G,AMP(25))
      CALL L1_1(W(1,3),W(1,1),G,ZERO, ZERO, W(1,29))
      CALL L1_1(W(1,6),W(1,29),G,ZERO, ZERO, W(1,30))
C     Amplitude(s) for diagram number 26
      CALL L1_0(W(1,5),W(1,26),W(1,30),G,AMP(26))
C     Amplitude(s) for diagram number 27
      CALL L1_0(W(1,12),W(1,26),W(1,29),G,AMP(27))
C     Amplitude(s) for diagram number 28
      CALL L1_0(W(1,10),W(1,2),W(1,30),G,AMP(28))
      CALL L1_1(W(1,3),W(1,14),G,ZERO, ZERO, W(1,31))
C     Amplitude(s) for diagram number 29
      CALL L1_0(W(1,16),W(1,2),W(1,31),G,AMP(29))
      CALL L1_1(W(1,18),W(1,1),G,ZERO, ZERO, W(1,32))
C     Amplitude(s) for diagram number 30
      CALL L1_0(W(1,6),W(1,26),W(1,32),G,AMP(30))
      CALL L1_1(W(1,5),W(1,29),G,ZERO, ZERO, W(1,33))
C     Amplitude(s) for diagram number 31
      CALL L1_0(W(1,16),W(1,2),W(1,33),G,AMP(31))
      CALL L1_1(W(1,5),W(1,2),G,ZERO, ZERO, W(1,34))
C     Amplitude(s) for diagram number 32
      CALL L1_0(W(1,3),W(1,34),W(1,17),G,AMP(32))
C     Amplitude(s) for diagram number 33
      CALL L1_0(W(1,11),W(1,34),W(1,7),G,AMP(33))
C     Amplitude(s) for diagram number 34
      CALL L1_0(W(1,3),W(1,34),W(1,19),G,AMP(34))
C     Amplitude(s) for diagram number 35
      CALL L1_0(W(1,3),W(1,34),W(1,21),G,AMP(35))
C     Amplitude(s) for diagram number 36
      CALL L1_0(W(1,23),W(1,34),W(1,20),G,AMP(36))
C     Amplitude(s) for diagram number 37
      CALL L1_0(W(1,4),W(1,34),W(1,28),G,AMP(37))
C     Amplitude(s) for diagram number 38
      CALL L1_0(W(1,4),W(1,34),W(1,30),G,AMP(38))
C     Amplitude(s) for diagram number 39
      CALL L1_0(W(1,16),W(1,34),W(1,29),G,AMP(39))
      CALL L1_1(W(1,3),W(1,20),G,ZERO, ZERO, W(1,35))
C     Amplitude(s) for diagram number 40
      CALL L1_0(W(1,12),W(1,2),W(1,35),G,AMP(40))
      CALL L1_1(W(1,11),W(1,1),G,ZERO, ZERO, W(1,36))
C     Amplitude(s) for diagram number 41
      CALL L1_0(W(1,6),W(1,34),W(1,36),G,AMP(41))
      CALL L1_1(W(1,4),W(1,29),G,ZERO, ZERO, W(1,37))
C     Amplitude(s) for diagram number 42
      CALL L1_0(W(1,12),W(1,2),W(1,37),G,AMP(42))
      CALL L1_1(W(1,6),W(1,2),G,ZERO, ZERO, W(1,38))
C     Amplitude(s) for diagram number 43
      CALL L1_0(W(1,3),W(1,38),W(1,22),G,AMP(43))
C     Amplitude(s) for diagram number 44
      CALL L1_0(W(1,11),W(1,38),W(1,14),G,AMP(44))
C     Amplitude(s) for diagram number 45
      CALL L1_0(W(1,3),W(1,38),W(1,24),G,AMP(45))
C     Amplitude(s) for diagram number 46
      CALL L1_0(W(1,3),W(1,38),W(1,25),G,AMP(46))
C     Amplitude(s) for diagram number 47
      CALL L1_0(W(1,18),W(1,38),W(1,20),G,AMP(47))
C     Amplitude(s) for diagram number 48
      CALL L1_0(W(1,4),W(1,38),W(1,32),G,AMP(48))
C     Amplitude(s) for diagram number 49
      CALL L1_0(W(1,4),W(1,38),W(1,33),G,AMP(49))
C     Amplitude(s) for diagram number 50
      CALL L1_0(W(1,10),W(1,38),W(1,29),G,AMP(50))
C     Amplitude(s) for diagram number 51
      CALL L1_0(W(1,5),W(1,38),W(1,36),G,AMP(51))

      DO I=1,NGRAPHS
        AMP2(I)=AMP2(I)+AMP(I)*DCONJG(AMP(I))
      ENDDO

 10   ZTEMP = (0.D0,0.D0)
      ZTEMP=ZTEMP+ONEPERM(P,NHEL,PERM)
      MATRIX=REAL(ZTEMP)

      RETURN
      END

C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM(P,NHEL,PM)
C     --------------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INTEGER    NJAMPS
      PARAMETER (NJAMPS=24)
      INTEGER    NPERMS
      PARAMETER (NPERMS=24)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NJAMPS)
      INTEGER PERMS(NEXTERNAL,NPERMS),IFERM(NPERMS),PERM(NEXTERNAL)
      DATA (PERMS(I,   1),I=1,6) / 1, 2, 3, 4, 5, 6/
      DATA (PERMS(I,   2),I=1,6) / 1, 2, 3, 5, 4, 6/
      DATA (PERMS(I,   3),I=1,6) / 1, 2, 4, 3, 5, 6/
      DATA (PERMS(I,   4),I=1,6) / 1, 2, 5, 4, 3, 6/
      DATA (PERMS(I,   5),I=1,6) / 1, 3, 2, 4, 5, 6/
      DATA (PERMS(I,   6),I=1,6) / 1, 4, 3, 2, 5, 6/
      DATA (PERMS(I,   7),I=1,6) / 1, 4, 5, 2, 3, 6/
      DATA (PERMS(I,   8),I=1,6) / 1, 4, 5, 3, 2, 6/
      DATA (PERMS(I,   9),I=1,6) / 1, 5, 3, 4, 2, 6/
      DATA (PERMS(I,  10),I=1,6) / 1, 5, 4, 2, 3, 6/
      DATA (PERMS(I,  11),I=1,6) / 2, 1, 3, 4, 5, 6/
      DATA (PERMS(I,  12),I=1,6) / 2, 3, 4, 5, 1, 6/
      DATA (PERMS(I,  13),I=1,6) / 3, 2, 1, 4, 5, 6/
      DATA (PERMS(I,  14),I=1,6) / 3, 4, 1, 2, 5, 6/
      DATA (PERMS(I,  15),I=1,6) / 3, 4, 2, 1, 5, 6/
      DATA (PERMS(I,  16),I=1,6) / 3, 4, 5, 2, 1, 6/
      DATA (PERMS(I,  17),I=1,6) / 4, 2, 3, 1, 5, 6/
      DATA (PERMS(I,  18),I=1,6) / 4, 3, 1, 2, 5, 6/
      DATA (PERMS(I,  19),I=1,6) / 4, 5, 2, 3, 1, 6/
      DATA (PERMS(I,  20),I=1,6) / 4, 5, 3, 1, 2, 6/
      DATA (PERMS(I,  21),I=1,6) / 5, 1, 2, 3, 4, 6/
      DATA (PERMS(I,  22),I=1,6) / 5, 2, 3, 4, 1, 6/
      DATA (PERMS(I,  23),I=1,6) / 5, 3, 4, 1, 2, 6/
      DATA (PERMS(I,  24),I=1,6) / 5, 4, 1, 2, 3, 6/
      DATA IFERM/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
     $ , 1, 1, 1, 1, 1, 1/
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 FLOW1
      EXTERNAL FLOW1
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,1))
      ENDDO
      JAMP(1)=IFERM(1)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,2))
      ENDDO
      JAMP(2)=IFERM(2)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,3))
      ENDDO
      JAMP(3)=IFERM(3)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,4))
      ENDDO
      JAMP(4)=IFERM(4)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,5))
      ENDDO
      JAMP(5)=IFERM(5)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,6))
      ENDDO
      JAMP(6)=IFERM(6)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,7))
      ENDDO
      JAMP(7)=IFERM(7)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,8))
      ENDDO
      JAMP(8)=IFERM(8)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,9))
      ENDDO
      JAMP(9)=IFERM(9)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,10))
      ENDDO
      JAMP(10)=IFERM(10)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,11))
      ENDDO
      JAMP(11)=IFERM(11)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,12))
      ENDDO
      JAMP(12)=IFERM(12)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,13))
      ENDDO
      JAMP(13)=IFERM(13)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,14))
      ENDDO
      JAMP(14)=IFERM(14)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,15))
      ENDDO
      JAMP(15)=IFERM(15)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,16))
      ENDDO
      JAMP(16)=IFERM(16)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,17))
      ENDDO
      JAMP(17)=IFERM(17)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,18))
      ENDDO
      JAMP(18)=IFERM(18)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,19))
      ENDDO
      JAMP(19)=IFERM(19)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,20))
      ENDDO
      JAMP(20)=IFERM(20)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,21))
      ENDDO
      JAMP(21)=IFERM(21)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,22))
      ENDDO
      JAMP(22)=IFERM(22)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,23))
      ENDDO
      JAMP(23)=IFERM(23)*FLOW1(P,NHEL,PERM)
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,24))
      ENDDO
      JAMP(24)=IFERM(24)*FLOW1(P,NHEL,PERM)
      ZTEMP = (0.D0,0.D0)
      ZTEMP = ZTEMP+1/648D0*JAMP(1)*DCONJG(3641D0*(JAMP(1))+572D0
     $ *(JAMP(4)+JAMP(6)+JAMP(13)+JAMP(16)+JAMP(22)+JAMP(24))
     $ +554D0*(JAMP(7)+JAMP(14))+454D0*(-JAMP(2)-JAMP(3)-JAMP(5)
     $ -JAMP(11)-JAMP(12)-JAMP(21))+428D0*(JAMP(8)+JAMP(9)+JAMP(10)
     $ +JAMP(15)+JAMP(17)+JAMP(18)+JAMP(19)+JAMP(20)+JAMP(23)))
      ONEPERM=ZTEMP

      JAMP2(1)=JAMP2(1)+JAMP(1)*DCONJG(JAMP(1))

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM(IPERM,PERM)
C     ------------------------------------
C     
C     Gives permutation number IPERM. 
C     Return value is the fermion factor due to PERM
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
C     
C     ARGUMENTS 
C     
      INTEGER IPERM,PERM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IFLAG
      LOGICAL OK
      INTEGER COMP(NEXTERNAL)
      DATA COMP/1,2,2,2,2,2/
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=I
      ENDDO
      I=1
      DO WHILE(I.LT.IPERM)
        CALL IPNEXT(PERM,NEXTERNAL,IFLAG)
        OK=.TRUE.
        DO J=1,NEXTERNAL
          IF(COMP(PERM(J)).NE.COMP(J))THEN
            OK=.FALSE.
            EXIT
          ENDIF
        ENDDO
        IF(OK) I=I+1
      ENDDO
      END


""" % misc.get_pkg_info()

        process_exporter.write_matrix_element_v4(\
            writers.FortranWriter(self.give_pos('test')),
            matrix_element,
            self.myfortranmodel)

        #print open(self.give_pos('test')).read()
        self.assertFileContains('test', goal_matrix_f)

