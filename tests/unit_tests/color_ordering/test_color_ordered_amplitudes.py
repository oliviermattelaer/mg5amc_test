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

"""Unit test library for the various base objects of the core library"""

import copy
from fractions import Fraction
import itertools
import logging
import math
import os
import unittest

import madgraph.core.base_objects as base_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.color_algebra as color
import madgraph.color_ordering.color_ordered_amplitudes as \
       color_ordered_amplitudes
import madgraph.color_ordering.color_ordered_export_v4 as \
       color_ordered_export_v4
import madgraph.iolibs.drawing_eps as draw
import madgraph.iolibs.export_v4 as export_v4

#===============================================================================
# ColorOrderedAmplitudeTest
#===============================================================================
class ColorOrderedAmplitudeTest(unittest.TestCase):
    """Test class for all functions related to the diagram generation"""

    def setUp(self):

        self.mypartlist = base_objects.ParticleList()
        self.myinterlist = base_objects.InteractionList()
        self.mymodel = base_objects.Model()
        self.myprocess = base_objects.Process()

        self.ref_dict_to0 = {}
        self.ref_dict_to1 = {}

        self.mycolorflow = color_ordered_amplitudes.ColorOrderedFlow()
        self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude()

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
                      'color': [color.ColorString([color.f(0, 1, -1),
                                                   color.f(2, 3, -1)]),
                                color.ColorString([color.f(2, 0, -1),
                                                   color.f(1, 3, -1)]),
                                color.ColorString([color.f(1, 2, -1),
                                                   color.f(0, 3, -1)])],
                      'lorentz':['gggg1', 'gggg2', 'gggg3'],
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
                      'color': [],
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
                      'color': [],
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

        self.myinterlist.append(base_objects.Interaction({
                      'id': 8,
                      'particles': base_objects.ParticleList(\
                                            [s,
                                             antis, \
                                             gamma]),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQED'},
                      'orders':{'QED':1}}))

        # Coupling of e to gamma

        self.myinterlist.append(base_objects.Interaction({
                      'id': 9,
                      'particles': base_objects.ParticleList(\
                                            [e, \
                                             antie, \
                                             gamma]),
                      'color': [],
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
                      'color': [],
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

        self.ref_dict_to0 = self.myinterlist.generate_ref_dict()[0]
        self.ref_dict_to1 = self.myinterlist.generate_ref_dict()[1]

    def test_color_ordered_model(self):
        """Test the ColorOrderedModel based on this model"""
        co_model = color_ordered_amplitudes.ColorOrderedModel(self.mymodel)
        singlet = base_objects.Particle({'name': 'g',
                                         'antiname': 'g',
                                         'spin': 3,
                                         'color': 1,
                                         'charge': 0.00,
                                         'mass': 'zero',
                                         'width': 'zero',
                                         'pdg_code': 4,
                                         'texname': 'g',
                                         'antitexname': 'g',
                                         'line': 'curly',
                                         'propagating': True,
                                         'is_part': True,
                                         'self_antipart': True})
        self.assertTrue(singlet in co_model.get('particles'))

        color_string = color.ColorString([color.T(0,1)])
        color_string.Nc_power = -1

        singlet_interaction1 = base_objects.Interaction(\
            {'id': 12,
             'particles': base_objects.ParticleList([co_model.get_particle(2),
                                                     co_model.get_particle(-2),
                                                     singlet]),
             'color': [color_string],
             'lorentz': ['L1'],
             'couplings': {(0, 0): 'SFACT*GQQ'},
             'orders': {'singlet_QCD': 1}})
        self.assertTrue(singlet_interaction1 in co_model.get('interactions'))

        singlet_interaction2 = base_objects.Interaction(\
            {'id': 13,
             'particles': base_objects.ParticleList([co_model.get_particle(1),
                                                     co_model.get_particle(-1),
                                                     singlet]),
             'color': [color_string],
             'lorentz': ['L1'],
             'couplings': {(0, 0): 'SFACT*GQQ'},
             'orders': {'singlet_QCD': 1}})
        self.assertTrue(singlet_interaction2 in co_model.get('interactions'))
        singlet_interaction3 = base_objects.Interaction(\
            {'id': 14,
             'particles': base_objects.ParticleList([co_model.get_particle(3),
                                                     co_model.get_particle(-3),
                                                     singlet]),
             'color': [color_string],
             'lorentz': ['L1'],
             'couplings': {(0, 0): 'SFACT*GQQ'},
             'orders': {'singlet_QCD': 1}})
        self.assertTrue(singlet_interaction3 in co_model.get('interactions'))

    def test_color_ordered_gluons(self):
        """Test the number of color ordered diagrams gg>ng with n up to 4"""

        goal_ndiags = [3, 10, 38,  154, 654,  2871,  12925]
        goal_nperms = [6, 24, 120, 720, 5040, 40320, 362880]

        # Time for 7 gluons: 46 s

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range (2, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            # Test the process after setup
            mycoleglist = base_objects.LegList([\
                color_ordered_amplitudes.ColorOrderedLeg(leg) for leg in myleglist])

            for i, leg in enumerate(mycoleglist):
                leg.set('color_ordering', {0: (i+1, i+1)})

            mycoproc = base_objects.Process({'legs':mycoleglist,
                                           'orders':{'QCD':ngluon},
                                           'model':self.mymodel})

            mycolorflow = self.myamplitude.get('color_flows')[0]

            self.assertEqual(mycolorflow.get('process'),
                             mycoproc)

            # Call generate_diagram and output number of diagrams
            ndiags = len(mycolorflow.get('diagrams'))

            #print "Number of diagrams for %d gluons: %d, cmp %d" % (ngluon,
            #                                                        ndiags,
            #                                                        goal_ndiags[ngluon-2])
            #print self.myamplitude.get('diagrams').nice_string()

            self.assertEqual(len(mycolorflow.get('permutations')),
                             goal_nperms[ngluon-2])

            self.assertEqual(ndiags, goal_ndiags[ngluon - 2])

    def test_color_ordered_uux_nglue(self):
        """Test the number of color flows and diagrams generated for uu~>gg
        """

        goal_nflows = [1, 1, 1, 1]
        goal_ndiagrams = [2, 6, 21, 81]
        goal_nperms = [2, 6, 24, 120]

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range (2, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon, 'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluon-2])
            mycolorflow = self.myamplitude.get('color_flows')[0]
            self.assertEqual(len(mycolorflow.\
                                 get('diagrams')), goal_ndiagrams[ngluon-2])
            self.assertEqual(len(mycolorflow.get('permutations')),
                             goal_nperms[ngluon-2])

    def test_color_ordered_uux_uuxng(self):
        """Test the number of color flows and diagrams for uu~>uu~+ng
        """
        goal_ndiags = [[1, 1], [3, 2, 3, 2], [10, 5, 11, 4, 10, 5],[37, 16, 41, 10, 41, 10, 37, 16]]
        goal_nflows = [2, 4, 6, 8]
        goal_nperms = [[2] * 2, [2] * 4, [4] * 6, [12] * 8]

        for ngluons in range(0, 3):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+2, 'QED': 0}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            perms = self.myamplitude.get('color_flows')[0].get('permutations')
            perms_0 = [base_objects.reorder_permutation(perm, perms[0]) for \
                       perm in perms]
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                self.assertEqual(len(flow.get('diagrams')),
                             goal_ndiags[ngluons][iflow])
                self.assertEqual(len(flow.get('permutations')),
                                 goal_nperms[ngluons][iflow])
                # Check if all perms are equal
                perms = flow.get('permutations')
                perms = [base_objects.reorder_permutation(perm, flow.get('permutations')[0]) for \
                       perm in flow.get('permutations')]
                self.assertEqual(perms, perms_0)
            
    def test_color_ordered_uux_ddxng(self):
        """Test the number of color flows and diagrams for uu~>dd~+ng
        """
        goal_ndiags = [[1, 1], [3, 3, 2, 2], [10, 11, 10, 5, 4, 5],
                       [37, 41, 41, 37, 16, 10, 10, 16]]
        goal_nflows = [2, 4, 6, 8]
        goal_nperms = [[1]*2, [1]*4, [2]*6, [6]*8]

        for ngluons in range(0, 3):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':1,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-1,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+2, 'QED': 0}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()
            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                  l in c.get('process').get('legs')]
            #    print c.nice_string()

            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                self.assertEqual(len(flow.get('diagrams')),
                                 goal_ndiags[ngluons][iflow])
                self.assertEqual(len(flow.get('permutations')),
                                 goal_nperms[ngluons][iflow])

    def test_color_ordered_uux_ddxssxng(self):
        """Test the number of color flows and diagrams for uu~>dd~ssx+ng
        """
        goal_ndiags = [[4, 4, 4, 6, 4, 4],
                       [16, 16, 14, 8, 16, 14, 14, 14, 14, 14, 8, 14, 16, 16, 8, 14, 16, 14],
                       [63, 69, 63, 50, 28, 20, 69, 69, 56, 28, 63, 50, 50, 56, 50, 38, 32, 38, 28, 28, 32, 32, 20, 38, 63, 69, 63, 20, 28, 50, 69, 69, 28, 56, 63, 50]]
        
        goal_nflows = [6, 18, 36]
        #goal_nflows = []
        #goal_ndiags = []
        for ngluons in range(0, 3):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':1,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-1,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':3,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-3,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+4, 'QED': 0}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()
            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('id'),
            #                                   l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            #goal_nflows.append(len(self.myamplitude.get('color_flows')))
            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            #diags = []
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                #plot = draw.MultiEpsDiagramDrawer(flow.get('diagrams'),
                #                                  "test.eps",
                #                                  model=self.mymodel)
                #plot.draw()
                #diags.append(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')))
                self.assertEqual(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')),
                            goal_ndiags[ngluons][iflow])
            #goal_ndiags.append(diags)
            #print goal_nflows
            #print "goal_ndiags = ",goal_ndiags

    def test_color_ordered_uux_uuxddxng(self):
        """Test number of color flows and diagrams for uu~>uu~dd~+ng
        """
        goal_ndiags =  [[4, 4, 4, 6, 4, 4],
                        [16, 14, 16, 14, 14, 14, 8, 14, 16, 8, 14, 14, 16, 8, 16, 14, 16, 14],
                        [63, 50, 69, 56, 63, 50, 50, 38, 28, 32, 20, 38, 69, 28, 69, 28, 56, 32, 28, 32, 63, 20, 50, 38, 63, 20, 69, 28, 63, 50, 69, 28, 69, 56, 63, 50]]
        
        goal_nflows = [6, 18, 36]
        #goal_nflows = []
        #goal_ndiags = []
        goal_nperms = [[2]*6,[2]*18,[4]*36]
        
        for ngluons in range(0, 1):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':1,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-1,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+4, 'QED': 0}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()
            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            #goal_nflows.append(len(self.myamplitude.get('color_flows')))
            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            #diags = []
            #perms = []
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                #diags.append(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')))
                self.assertEqual(len(flow.get('diagrams')),
                                 goal_ndiags[ngluons][iflow])
                #perms.append(len(flow.get('permutations')))
                self.assertEqual(len(flow.get('permutations')),
                                 goal_nperms[ngluons][iflow])
            #goal_ndiags.append(diags)
            #goal_nperms.append(perms)
            #print goal_nflows
            #print "goal_ndiags = ",goal_ndiags
            #print "goal_nperms = ",goal_nperms

    def test_color_ordered_ddx_uuxuuxng(self):
        """Test number of color flows and diagrams for dd~>uu~uu~+ng
        """
        goal_ndiags =  [[4, 4, 4, 4, 4, 6],
                        [16, 14, 16, 8, 16, 14, 16, 14, 16, 14, 8, 14, 14, 14, 16, 8, 14, 14],
                        [63, 50, 69, 28, 63, 20, 69, 56, 69, 28, 63, 50, 63, 50, 69, 56, 63, 50, 20, 38, 28, 32, 50, 38, 69, 28, 69, 28, 28, 32, 56, 32, 63, 20, 50, 38]]
        
        goal_nflows = [6, 18, 36]
        goal_nperms = [[2]*6,[2]*18,[4]*36]        
        #goal_nflows = []
        #goal_ndiags = []

        for ngluons in range(0, 1):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':1,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-1,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+4, 'QED': 0}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()
            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            #goal_nflows.append(len(self.myamplitude.get('color_flows')))
            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            #diags = []
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                #diags.append(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')))
                self.assertEqual(len(flow.get('diagrams')),
                             goal_ndiags[ngluons][iflow])
                self.assertEqual(len(flow.get('permutations')),
                                 goal_nperms[ngluons][iflow])
            #goal_ndiags.append(diags)
            #print goal_nflows
            #print "goal_ndiags = ",goal_ndiags

    def test_color_ordered_uux_epem_nglue(self):
        """Test the number of color flows and diagrams generated for uu~>e+e-ng
        """

        goal_nflows = [1, 1, 1, 1, 1]
        goal_ndiagrams = [1, 2, 5, 16, 58]
        goal_perms = [[[1, 2, 3, 4]],
                      [[1, 2, 3, 4, 5]],
                      [[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 6, 5]],
                      [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 7, 6], [1, 2, 3, 4, 6, 5, 7], [1, 2, 3, 4, 6, 7, 5], [1, 2, 3, 4, 7, 5, 6], [1, 2, 3, 4, 7, 6, 5]]]
        # Test 0 to 4 gluons in the final state
        for ngluon in range (0, 4):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':-11}),
                base_objects.Leg({'id':11})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon, 'QED': 2},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluon])
            self.assertEqual(len(self.myamplitude.get('color_flows')[0].\
                                 get('diagrams')), goal_ndiagrams[ngluon])
            self.assertEqual(self.myamplitude.get('color_flows')[0].get('permutations'),
                             goal_perms[ngluon])
                             
            
    def test_color_ordered_uux_uuxepemng(self):
        """Test color flows and diagrams for uu~>uu~e+e-+ng with n up to 2
        """
        goal_ndiags = [[4, 4], [14, 10, 14, 10], [50, 28, 56, 24, 50, 28]]
        goal_nflows = [2, 4, 6]

        for ngluons in range(0, 2):

            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':11,
                                             'state':True}))
            myleglist.append(base_objects.Leg({'id':-11,
                                             'state':True}))
            myleglist.extend([base_objects.Leg({'id':21,
                                                 'state':True})] * ngluons)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel,
                                           'orders':{'QCD':ngluons+2, 'QED': 2}})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()
            #for c in self.myamplitude.get('color_flows'):
            #    print "color flow process: ",[(l.get('number'), l.get('color_ordering')) for \
            #                                   l in c.get('process').get('legs')]
            #    print c.nice_string()

            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluons])
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                self.assertEqual(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')),
                             goal_ndiags[ngluons][iflow])
                
    def test_color_ordered_gg_h_nglue(self):
        """Test the number of color ordered diagrams gg>h+ng with n up to 3"""

        mypartlist = base_objects.ParticleList()
        myinterlist = base_objects.InteractionList()
        mymodel = base_objects.Model()
        myprocess = base_objects.Process()
        myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude()

        # A gluon
        mypartlist.append(base_objects.Particle({'name':'g',
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

        # A Higgs
        mypartlist.append(base_objects.Particle({'name':'h',
                      'antiname':'h',
                      'spin':1,
                      'color':1,
                      'mass':'zero',
                      'width':'zero',
                      'texname':'H',
                      'antitexname':'H',
                      'line':'wavy',
                      'charge':0.,
                      'pdg_code':25,
                      'propagating':True,
                      'is_part':True,
                      'self_antipart':True}))

        # 3 gluon vertex
        myinterlist.append(base_objects.Interaction({
                      'id': 1,
                      'particles': base_objects.ParticleList(\
                                            [mypartlist[0]] * 3),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'G'},
                      'orders':{'QCD':1}}))

        # 4 gluon vertex
        myinterlist.append(base_objects.Interaction({
                      'id': 2,
                      'particles': base_objects.ParticleList(\
                                            [mypartlist[0]] * 4),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'G^2'},
                      'orders':{'QCD':2}}))

        # Gluon couplings to Higgs
        myinterlist.append(base_objects.Interaction({
                      'id': 3,
                      'particles': base_objects.ParticleList(\
                                            [mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[1]]),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'HIG':1}}))

        myinterlist.append(base_objects.Interaction({
                      'id': 4,
                      'particles': base_objects.ParticleList(\
                                            [mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[1]]),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'HIG':1, 'QCD':1}}))

        myinterlist.append(base_objects.Interaction({
                      'id': 5,
                      'particles': base_objects.ParticleList(\
                                            [mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[0], \
                                             mypartlist[1]]),
                      'color': [],
                      'lorentz':['L1'],
                      'couplings':{(0, 0):'GQQ'},
                      'orders':{'HIG':1, 'QCD':2}}))

        mymodel.set('particles', mypartlist)
        mymodel.set('interactions', myinterlist)

        goal_ndiags = [1, 4, 19, 90]

        # Test 0-4 gluons in the final state
        for ngluon in range (0, 4):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.append(base_objects.Leg({'id':25,
                                                'state':True}))

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon, 'HIG':1},
                                           'model':mymodel})

            myamplitude.set('process', myproc)

            myamplitude.generate_diagrams()

            # Test the process after setup
            mycoleglist = base_objects.LegList([\
                color_ordered_amplitudes.ColorOrderedLeg(leg) for leg in myleglist])

            for i, leg in enumerate(mycoleglist):
                leg.set('color_ordering', {0: (i+1, i+1)})

            mycolorflow = myamplitude.get('color_flows')[0]

            # Call generate_diagram and output number of diagrams
            ndiags = len(mycolorflow.get('diagrams'))

            #print mycolorflow.get('process').nice_string()
            #print "Number of diagrams for %d gluons: %d, cmp %d" % (ngluon,
            #                                                        ndiags,
            #                                                        goal_ndiags[ngluon])
            #print myamplitude.get('diagrams').nice_string()

            self.assertEqual(ndiags, goal_ndiags[ngluon])

    def test_color_ordered_multi_process(self):
        """Test the number of color ordered diagrams gg>ng with n up to 4"""

        goal_ndiags = [3, 10]

        ngluon = 2

        # Create the amplitude
        myleglist = base_objects.MultiLegList([base_objects.MultiLeg({'ids':[21],
                                          'state':False})] * 2)

        myleglist.extend([base_objects.MultiLeg({'ids':[21],
                                          'state':True})] * ngluon)

        myproc = base_objects.ProcessDefinition({'legs':myleglist,
                                       'orders':{'QCD':ngluon},
                                       'model':self.mymodel})

        mymultiproc = color_ordered_amplitudes.ColorOrderedMultiProcess(
            myproc)

        self.assertEqual(len(mymultiproc.get('amplitudes')), 1)

        myamplitude = mymultiproc.get('amplitudes')[0]

        self.assertEqual(len(myamplitude.get('color_flows')), 1)
        
        mycolorflow = myamplitude.get('color_flows')[0]

        # Call generate_diagram and output number of diagrams
        ndiags = len(mycolorflow.get('diagrams'))

        self.assertEqual(ndiags, goal_ndiags[ngluon - 2])

#===============================================================================
# COHelasMatrixElementTest
#===============================================================================
class COHelasMatrixElementTest(unittest.TestCase):
    """Test class for functions related to the B-G matrix elements"""

    mypartlist = base_objects.ParticleList()
    myinterlist = base_objects.InteractionList()
    mymodel = base_objects.Model()
    myprocess = base_objects.Process()

    ref_dict_to0 = {}
    ref_dict_to1 = {}

    mycolorflow = color_ordered_amplitudes.ColorOrderedFlow()
    myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude()

    def setUp(self):

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

        self.ref_dict_to0 = self.myinterlist.generate_ref_dict()[0]
        self.ref_dict_to1 = self.myinterlist.generate_ref_dict()[1]

    def test_matrix_element_gluons(self):
        """Test the matrix element for all-gluon amplitudes"""

        # Time for gg>7g: 1980 s
        goal_wavefunctions = [6, 10, 27, 49, 96, 153]
        goal_amplitudes = [4, 15, 32, 70, 120, 210]
        goal_matrix =  [{(0, 1): (Fraction(-1, 3), 0),
                         (0, 0): (Fraction(19, 6), 0),
                         (0, 5): (Fraction(2, 3), 0),
                         (0, 4): (Fraction(-1, 3), 0),
                         (0, 3): (Fraction(-1, 3), 0),
                         (0, 2): (Fraction(-1, 3), 0)},
                        {(0, 20): (0, 0), (0, 14): (Fraction(17, 27), 0),
                         (0, 7): (0, 0), (0, 16): (0, 0), (0, 10): (0, 0),
                         (0, 3): (0, 0), (0, 21): (Fraction(17, 27), 0),
                         (0, 15): (0, 0), (0, 17): (Fraction(17, 27), 0),
                         (0, 11): (0, 0), (0, 4): (0, 0),
                         (0, 22): (Fraction(17, 27), 0),
                         (0, 0): (Fraction(455, 108), 0),
                         (0, 18): (Fraction(-29, 54), 0),
                         (0, 12): (0, 0), (0, 5): (Fraction(17, 27), 0),
                         (0, 8): (0, 0), (0, 23): (0, 0),
                         (0, 1): (Fraction(-29, 54), 0),
                         (0, 19): (0, 0), (0, 13): (0, 0),
                         (0, 6): (Fraction(-29, 54), 0),
                         (0, 9): (Fraction(-29, 54), 0),
                         (0, 2): (Fraction(-29, 54), 0)},
                        {(0, 86): (0, 0), (0, 20): (0, 0), (0, 89): (0, 0),
                         (0, 7): (0, 0), (0, 76): (0, 0), (0, 10): (0, 0),
                         (0, 61): (0, 0), (0, 98): (0, 0), (0, 32): (0, 0),
                         (0, 83): (0, 0), (0, 17): (Fraction(107, 162), 0),
                         (0, 70): (0, 0), (0, 4): (0, 0), (0, 73): (0, 0),
                         (0, 55): (0, 0), (0, 58): (0, 0), (0, 111): (0, 0),
                         (0, 45): (0, 0), (0, 80): (Fraction(107, 162), 0),
                         (0, 30): (0, 0), (0, 67): (0, 0),
                         (0, 1): (Fraction(-227, 324), 0),
                         (0, 118): (0, 0), (0, 52): (0, 0), (0, 39): (0, 0),
                         (0, 108): (0, 0), (0, 42): (0, 0),
                         (0, 93): (Fraction(107, 162), 0),
                         (0, 27): (0, 0), (0, 64): (0, 0),
                         (0, 14): (Fraction(143, 162), 0),
                         (0, 115): (0, 0), (0, 49): (0, 0), (0, 102): (0, 0),
                         (0, 36): (0, 0), (0, 105): (Fraction(143, 162), 0),
                         (0, 87): (0, 0), (0, 21): (Fraction(107, 162), 0),
                         (0, 90): (0, 0), (0, 24): (Fraction(-227, 324), 0),
                         (0, 77): (0, 0), (0, 11): (0, 0),
                         (0, 112): (Fraction(107, 162), 0),
                         (0, 62): (Fraction(107, 162), 0),
                         (0, 99): (0, 0), (0, 33): (Fraction(-227, 324), 0),
                         (0, 84): (Fraction(107, 162), 0),
                         (0, 18): (0, 0), (0, 71): (0, 0),
                         (0, 5): (Fraction(143, 162), 0),
                         (0, 74): (0, 0), (0, 8): (0, 0), (0, 59): (0, 0),
                         (0, 96): (Fraction(-227, 324), 0),
                         (0, 46): (0, 0), (0, 81): (0, 0), (0, 31): (0, 0),
                         (0, 68): (0, 0), (0, 2): (Fraction(-227, 324), 0),
                         (0, 119): (0, 0), (0, 53): (0, 0), (0, 56): (0, 0),
                         (0, 109): (0, 0), (0, 43): (0, 0),
                         (0, 94): (Fraction(107, 162), 0),
                         (0, 28): (0, 0), (0, 65): (Fraction(143, 162), 0),
                         (0, 15): (0, 0), (0, 116): (0, 0), (0, 50): (0, 0),
                         (0, 103): (0, 0), (0, 37): (0, 0), (0, 106): (0, 0),
                         (0, 40): (0, 0), (0, 22): (Fraction(107, 162), 0),
                         (0, 91): (0, 0), (0, 25): (0, 0), (0, 78): (0, 0),
                         (0, 12): (0, 0), (0, 113): (0, 0), (0, 63): (0, 0),
                         (0, 100): (0, 0), (0, 34): (0, 0), (0, 85): (0, 0),
                         (0, 19): (0, 0), (0, 88): (0, 0),
                         (0, 6): (Fraction(-227, 324), 0),
                         (0, 75): (0, 0), (0, 9): (0, 0),
                         (0, 60): (Fraction(277, 324), 0),
                         (0, 97): (0, 0), (0, 47): (0, 0), (0, 82): (0, 0),
                         (0, 16): (Fraction(277, 324), 0),
                         (0, 69): (0, 0), (0, 3): (0, 0), (0, 72): (0, 0),
                         (0, 54): (Fraction(143, 162), 0),
                         (0, 57): (0, 0), (0, 110): (0, 0), (0, 44): (0, 0),
                         (0, 95): (0, 0), (0, 29): (0, 0), (0, 66): (0, 0),
                         (0, 0): (Fraction(3641, 648), 0),
                         (0, 117): (0, 0), (0, 51): (0, 0), (0, 38): (0, 0),
                         (0, 107): (0, 0), (0, 41): (0, 0), (0, 23): (0, 0),
                         (0, 92): (0, 0), (0, 26): (0, 0), (0, 79): (0, 0),
                         (0, 13): (0, 0), (0, 114): (Fraction(143, 162), 0),
                         (0, 48): (0, 0), (0, 101): (0, 0), (0, 35): (0, 0),
                         (0, 104): (0, 0)}]

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(2,5):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon},
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            #print "Generated diagrams for ", myproc.nice_string()

            matrix_element = color_ordered_amplitudes.COHelasMatrixElement(\
                self.myamplitude, gen_color=3, optimization=3)

            #print "Generated matrix element"

            helas_flow = matrix_element.get('color_flows')[0]
            #print "permutations: ",helas_flow.get('permutations')
            if ngluon < 5:
                self.assertEqual(matrix_element.get('color_matrix').\
                                 col_matrix_fixed_Nc, goal_matrix[ngluon-2])

            # Check that we get back the correct number of diagrams
            # when we do get_base_amplitude. This is a very powerful
            # check that the BG currents are correct.
            diagrams = self.myamplitude.get('color_flows')[0].get('diagrams')
            base_amplitude = helas_flow.get_base_amplitude().get('diagrams')
            self.assertEqual(len(diagrams),len(base_amplitude))

            diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                            for d in diagrams]
            base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                         for d in base_amplitude]
            base_tag_copy = copy.copy(base_tags)
            for tag in diagram_tags:
                self.assertTrue(tag in base_tags)
                base_tags.remove(tag)
            self.assertEqual(base_tags, [])

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(helas_flow))
            #print "For ",ngluon," FS gluons, there are ",\
            #      len(helas_flow.get_all_amplitudes()),' amplitudes and ',\
            #      len(helas_flow.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            # Test that wavefunctions that have been combined are
            # not in any amplitudes
            comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                            for w in helas_flow.get_all_wavefunctions()\
                            if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
            amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                           for amp in helas_flow.get_all_amplitudes()],\
                          []))
            self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

            # Test that all wavefunctions are used
            wf_numbers = [w.get('number') for w in \
                          helas_flow.get_all_wavefunctions()]
            mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                           for a in helas_flow.get_all_amplitudes() + \
                                     helas_flow.get_all_wavefunctions()],\
                          []))
            self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

            #print helas_flow.get_base_amplitude().nice_string()

            #print "wavefunctions: ", len(helas_flow.get_all_wavefunctions())
            self.assertEqual(len(helas_flow.get_all_wavefunctions()),
                             goal_wavefunctions[ngluon-2])
            #print "amplitudes: ", len(helas_flow.get_all_amplitudes())
            self.assertEqual(len(helas_flow.get_all_amplitudes()),
                             goal_amplitudes[ngluon-2])
            
            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(helas_flow))

    def test_matrix_element_gluons_non_optimized(self):
        """Test the matrix element for all-gluon (non-optimized)"""

        # Time for gg>7g: 314 s

        goal_wavefunctions = [6, 10, 24, 42, 108, 189]
        goal_amplitudes = [4, 15, 68, 322, 1608, 8283]

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(2,5):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon},
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            matrix_element = color_ordered_amplitudes.COHelasMatrixElement(\
                self.myamplitude, gen_color=False, optimization=1)

            mycolorflow = matrix_element.get('color_flows')[0]

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(mycolorflow))
            #print "For ",ngluon," FS gluons, there are ",\
            #      len(mycolorflow.get_all_amplitudes()),' amplitudes and ',\
            #      len(mycolorflow.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            self.assertEqual(len(mycolorflow.get_all_wavefunctions()),
                             goal_wavefunctions[ngluon-2])
            self.assertEqual(len(mycolorflow.get_all_amplitudes()),
                             goal_amplitudes[ngluon-2])
            
            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(mycolorflow))

            #print mycolorflow.get_base_amplitude().nice_string()

    def test_matrix_element_photons(self):
        """Test the matrix element for uu~>na"""

        goal_amplitudes = [2, 6, 12, 30, 60, 140]
        goal_wavefunctions = [6, 11, 32, 77, 190, 429]
        #goal_amplitudes = []
        #goal_wavefunctions = []

        # Test 2, 3, 4 and 5 photons in the final state
        for nphoton in range(2, 5):

            # Create the amplitude
            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))

            myleglist.extend([base_objects.Leg({'id':22,
                                              'state':True})] * nphoton)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':0},
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            mycolorflow = self.myamplitude.get('color_flows')[0]

            matrix_element = color_ordered_amplitudes.COHelasFlow(\
                mycolorflow, gen_color=False, optimization = 3)

            # Check that we get back the correct number of diagrams
            # when we do get_base_amplitude. This is a very powerful
            # check that the BG currents are correct.
            diagrams = mycolorflow.get('diagrams')
            base_amplitude = matrix_element.get_base_amplitude().get('diagrams')
            self.assertEqual(len(diagrams),len(base_amplitude))

            diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                            for d in diagrams]
            base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                         for d in base_amplitude]
            base_tag_copy = copy.copy(base_tags)
            for tag in diagram_tags:
                self.assertTrue(tag in base_tags)
                base_tags.remove(tag)
            self.assertEqual(base_tags, [])

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(matrix_element))
            #print "For ",nphoton," FS photons, there are ",\
            #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
            #      len(matrix_element.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            # Test that wavefunctions that have been combined are
            # not in any amplitudes
            comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                            for w in matrix_element.get_all_wavefunctions()\
                            if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
            amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                           for amp in matrix_element.get_all_amplitudes()],\
                          []))
            self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

            # Test that all wavefunctions are used
            wf_numbers = [w.get('number') for w in \
                          matrix_element.get_all_wavefunctions()]
            mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                           for a in matrix_element.get_all_amplitudes() + \
                                     matrix_element.get_all_wavefunctions()],\
                          []))
            self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

            self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                             goal_wavefunctions[nphoton-2])
            self.assertEqual(len(matrix_element.get_all_amplitudes()),
                             goal_amplitudes[nphoton-2])

            #goal_wavefunctions.append(len(matrix_element.get_all_wavefunctions()))
            #goal_amplitudes.append(len(matrix_element.get_all_amplitudes()))
            #print "goal_wavefunctions = ",goal_wavefunctions
            #print "goal_amplitudes = ",goal_amplitudes

            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

    def test_matrix_element_photons_non_optimized(self):
        """Test the matrix element for uu~>na (non-optimized)"""

        goal_amplitudes = [2, 6, 24, 120, 720, 5040]
        goal_wavefunctions = [6, 11, 26, 57, 200, 527]
        #goal_amplitudes = []
        #goal_wavefunctions = []

        # Test 2, 3, 4 and 5 photons in the final state
        for nphoton in range(2, 5):

            # Create the amplitude
            myleglist = base_objects.LegList()

            myleglist.append(base_objects.Leg({'id':2,
                                             'state':False}))
            myleglist.append(base_objects.Leg({'id':-2,
                                             'state':False}))

            myleglist.extend([base_objects.Leg({'id':22,
                                              'state':True})] * nphoton)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            mycolorflow = self.myamplitude.get('color_flows')[0]

            matrix_element = color_ordered_amplitudes.COHelasFlow(\
                mycolorflow, gen_color=False, optimization = 1)

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(matrix_element))
            #print "For ",nphoton," FS photons, there are ",\
            #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
            #      len(matrix_element.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                             goal_wavefunctions[nphoton-2])
            self.assertEqual(len(matrix_element.get_all_amplitudes()),
                             goal_amplitudes[nphoton-2])

            #goal_wavefunctions.append(len(matrix_element.get_all_wavefunctions()))
            #goal_amplitudes.append(len(matrix_element.get_all_amplitudes()))
            #print "goal_wavefunctions = ",goal_wavefunctions
            #print "goal_amplitudes = ",goal_amplitudes

            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))


    def test_matrix_element_uux_nglue(self):
        """Test color ordered matrix elements for uu~>ng
        """

        goal_wavefunctions =  [6, 10, 21, 41, 70, 121]
        goal_amplitudes =  [2, 7, 18, 38, 76, 130]
        #goal_amplitudes = []
        #goal_wavefunctions = []

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(2, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            matrix_element = color_ordered_amplitudes.COHelasMatrixElement(\
                self.myamplitude, gen_color=False, optimization=3)

            mycolorflow = matrix_element.get('color_flows')[0]

            # Check that we get back the correct number of diagrams
            # when we do get_base_amplitude. This is a very powerful
            # check that the BG currents are correct.
            diagrams = self.myamplitude.get('color_flows')[0].get('diagrams')
            base_amplitude = mycolorflow.get_base_amplitude().get('diagrams')
            self.assertEqual(len(diagrams),len(base_amplitude))

            diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                            for d in diagrams]
            base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                         for d in base_amplitude]
            base_tag_copy = copy.copy(base_tags)
            for tag in diagram_tags:
                self.assertTrue(tag in base_tags)
                base_tags.remove(tag)
            self.assertEqual(base_tags, [])

            for cf in matrix_element.get('color_flows')[1:]:
                self.assertEqual(cf.get('permutations'),
                                 mycolorflow.get('permutations'))

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(matrix_element))

            for d in mycolorflow.get('diagrams'):
                self.assertTrue(len(d.get('amplitudes')) > 0)

            # Test that wavefunctions that have been combined are
            # not in any amplitudes
            comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                            for w in mycolorflow.get_all_wavefunctions()\
                            if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
            amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                           for amp in mycolorflow.get_all_amplitudes()],\
                          []))
            self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

            # Test that all wavefunctions are used
            wf_numbers = [w.get('number') for w in \
                          mycolorflow.get_all_wavefunctions()]
            mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                           for a in mycolorflow.get_all_amplitudes() + \
                                     mycolorflow.get_all_wavefunctions()],\
                          []))
            self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

            #print "For ",ngluon," FS gluons, there are ",\
            #      len(mycolorflow.get_all_amplitudes()),' amplitudes and ',\
            #      len(mycolorflow.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            self.assertEqual(len(mycolorflow.get_all_amplitudes()),
                             goal_amplitudes[ngluon-2])
            self.assertEqual(len(mycolorflow.get_all_wavefunctions()),
                             goal_wavefunctions[ngluon-2])

            #goal_wavefunctions.append(len(mycolorflow.get_all_wavefunctions()))
            #goal_amplitudes.append(len(mycolorflow.get_all_amplitudes()))
            #print "goal_wavefunctions = ",goal_wavefunctions
            #print "goal_amplitudes = ",goal_amplitudes

            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(mycolorflow))

    def test_matrix_element_uux_nglue_non_optimized(self):
        """Test color ordered matrix elements for uu~>ng (non-optimized)
        """

        goal_amplitudes = [2, 7, 28, 126]
        goal_wavefunctions = [6, 10, 18, 34]

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(2, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            matrix_element = color_ordered_amplitudes.COHelasMatrixElement(\
                self.myamplitude, gen_color=False, optimization=1)

            mycolorflow = matrix_element.get('color_flows')[0]

            for cf in matrix_element.get('color_flows')[1:]:
                self.assertEqual(cf.get('permutations'),
                                 mycolorflow.get('permutations'))

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(matrix_element))

            for d in mycolorflow.get('diagrams'):
                self.assertTrue(len(d.get('amplitudes')) > 0)

            # Test that wavefunctions that have been combined are
            # not in any amplitudes
            comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                            for w in mycolorflow.get_all_wavefunctions()\
                            if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
            amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                           for amp in mycolorflow.get_all_amplitudes()],\
                          []))
            self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

            # Test that all wavefunctions are used
            wf_numbers = [w.get('number') for w in \
                          mycolorflow.get_all_wavefunctions()]
            mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                           for a in mycolorflow.get_all_amplitudes() + \
                                     mycolorflow.get_all_wavefunctions()],\
                          []))
            self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

            #print "For ",ngluon," FS gluons, there are ",\
            #      len(mycolorflow.get_all_amplitudes()),' amplitudes and ',\
            #      len(mycolorflow.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            self.assertEqual(len(mycolorflow.get_all_amplitudes()),
                             goal_amplitudes[ngluon-2])
            self.assertEqual(len(mycolorflow.get_all_wavefunctions()),
                             goal_wavefunctions[ngluon-2])

            # Test JAMP (color amplitude) output
            #print '\n'.join(export_v4.get_JAMP_lines(mycolorflow))

    def test_matrix_element_gu_gunglue(self):
        """Test color ordered matrix elements for uu~>ng
        """

        goal_amplitudes = [2, 7, 16, 38]
        goal_wavefunctions = [6, 10, 23, 41]

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(0, 4):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':21, 'state':False}),
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':21, 'state':True}),
                base_objects.Leg({'id':2, 'state':True})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            matrix_element = color_ordered_amplitudes.COHelasMatrixElement(\
                self.myamplitude, gen_color=False, optimization=3)

            mycolorflow = matrix_element.get('color_flows')[0]

            # Check that we get back the correct number of diagrams
            # when we do get_base_amplitude. This is a very powerful
            # check that the BG currents are correct.
            diagrams = self.myamplitude.get('color_flows')[0].get('diagrams')
            base_amplitude = mycolorflow.get_base_amplitude().get('diagrams')
            self.assertEqual(len(diagrams),len(base_amplitude))

            diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                            for d in diagrams]
            base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                         for d in base_amplitude]
            base_tag_copy = copy.copy(base_tags)
            for tag in diagram_tags:
                self.assertTrue(tag in base_tags)
                base_tags.remove(tag)
            self.assertEqual(base_tags, [])

            for cf in matrix_element.get('color_flows')[1:]:
                self.assertEqual(cf.get('permutations'),
                                 mycolorflow.get('permutations'))

            #print "\n".join(\
            #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
            #    get_matrix_element_calls(matrix_element))

            for d in mycolorflow.get('diagrams'):
                self.assertTrue(len(d.get('amplitudes')) > 0)

            # Test that wavefunctions that have been combined are
            # not in any amplitudes
            comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                            for w in mycolorflow.get_all_wavefunctions()\
                            if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
            amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                           for amp in mycolorflow.get_all_amplitudes()],\
                          []))
            self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

            # Test that all wavefunctions are used
            wf_numbers = [w.get('number') for w in \
                          mycolorflow.get_all_wavefunctions()]
            mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                           for a in mycolorflow.get_all_amplitudes() + \
                                     mycolorflow.get_all_wavefunctions()],\
                          []))
            self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

            #print "For ",ngluon," FS gluons, there are ",\
            #      len(mycolorflow.get_all_amplitudes()),' amplitudes and ',\
            #      len(mycolorflow.get_all_wavefunctions()),\
            #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
            #      ' diagrams'

            self.assertEqual(len(mycolorflow.get_all_amplitudes()),
                             goal_amplitudes[ngluon])
            self.assertEqual(len(mycolorflow.get_all_wavefunctions()),
                             goal_wavefunctions[ngluon])

            # Test JAMP (color amplitude) output
            #exporter = color_ordered_export_v4.ProcessExporterFortranCO()
            #print '\n'.join(exporter.get_JAMP_lines(mycolorflow))

    def test_matrix_element_uux_ddxng(self):
        """Test color flow matrix elements for uu~>dd~ng
        """

        goal_amplitudes = [[1, 1],[3,3,2,2],[8,8,9,4,2,5],
                           [20, 21, 21, 20, 8, 3, 3, 8]]
        goal_wavefunctions = [[5, 5],[9,9,8,8],[18,18,16,14,13,12],
                              [32, 32, 32, 32, 25, 21, 21, 25]]
        goal_nflows = [2, 4, 6, 8]

        #goal_wavefunctions = []
        #goal_amplitudes = []
        #goal_nflows = []

        # Test 0-3 gluons in the final state
        for ngluon in range(0, 3):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            #diags=[]
            #amps=[]
            #wfs=[]

            for iflow, mycolorflow in \
                enumerate(self.myamplitude.get('color_flows')):

                matrix_element = color_ordered_amplitudes.COHelasFlow(\
                    mycolorflow, gen_color=False, optimization = 3)

                #print "\n".join(\
                #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
                #    get_matrix_element_calls(matrix_element))
                # Test JAMP (color amplitude) output
                #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

                for d in matrix_element.get('diagrams'):
                    self.assertTrue(len(d.get('amplitudes')) > 0)

                # Check that we get back the correct number of diagrams
                # when we do get_base_amplitude. This is a very powerful
                # check that the BG currents are correct.
                diagrams = mycolorflow.get('diagrams')
                base_amplitude = matrix_element.get_base_amplitude().get('diagrams')
                diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                                for d in diagrams]
                base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                             for d in base_amplitude]
                base_tag_copy = copy.copy(base_tags)
                
                #print 'diagrams: ',diagrams.nice_string()
                #print 'matrix_element: ',base_amplitude.nice_string()

                for tag in diagram_tags:
                    self.assertTrue(tag in base_tags)
                    base_tags.remove(tag)
                    
                self.assertEqual(base_tags, [])

                # Test that wavefunctions that have been combined are
                # not in any amplitudes
                comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                                for w in matrix_element.get_all_wavefunctions()\
                                if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
                amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                               for amp in matrix_element.get_all_amplitudes()],\
                              []))
                self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

                # Test that all wavefunctions are used
                wf_numbers = [w.get('number') for w in \
                              matrix_element.get_all_wavefunctions()]
                mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                               for a in matrix_element.get_all_amplitudes() + \
                                         matrix_element.get_all_wavefunctions()],\
                              []))
                self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

                #print "For ",ngluon," FS gluons, ",iflow+1,"th flow, there are ",\
                #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
                #      len(matrix_element.get_all_wavefunctions()),\
                #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
                #      ' diagrams'

                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))

                self.assertEqual(len(matrix_element.get_all_amplitudes()),
                                 goal_amplitudes[ngluon][iflow])
                self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                                 goal_wavefunctions[ngluon][iflow])
                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))

            #goal_wavefunctions.append(wfs)
            #goal_amplitudes.append(amps)

            #print goal_nflows
            #print goal_wavefunctions
            #print goal_amplitudes

    def test_matrix_element_uux_ddxng_non_optimized(self):
        """Test color flow matrix element for uu~>dd~ng witout optimization
        """

        goal_amplitudes = [[1, 1], [3, 3, 2, 2], [11, 12, 11, 5, 4, 5], [44, 50, 50, 44, 18, 10, 10, 18]]
        goal_wavefunctions = [[5, 5], [9, 9, 8, 8], [16, 16, 15, 13, 12, 12], [33, 33, 34, 30, 26, 19, 18, 24]]
        goal_ndiags = [[1, 1], [3, 3, 2, 2], [10, 11, 10, 5, 4, 5], [35, 41, 41, 35, 16, 10, 10, 16]]
        goal_nflows = [2, 4, 6, 8]

        #goal_wavefunctions = []
        #goal_amplitudes = []
        #goal_ndiags = []
        #goal_nflows = []
        # Test 0, 1, 2 and 3 gluons in the final state
        for ngluon in range(0, 3):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            #goal_nflows.append(len(self.myamplitude.get('color_flows')))
            #diags=[]
            #amps=[]
            #wfs=[]

            for iflow, mycolorflow in \
                enumerate(self.myamplitude.get('color_flows')):

                matrix_element = color_ordered_amplitudes.COHelasFlow(\
                    mycolorflow, gen_color=False, optimization = 1)

                # Check that we get back the correct number of diagrams
                # when we do get_base_amplitude. This is a very powerful
                # check that the BG currents are correct.
                self.assertEqual(len(matrix_element.get_base_amplitude().get('diagrams')),
                                 len(mycolorflow.get('diagrams')))

                #print "\n".join(\
                #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
                #    get_matrix_element_calls(matrix_element))

                for d in matrix_element.get('diagrams'):
                    self.assertTrue(len(d.get('amplitudes')) > 0)

                #print "For ",ngluon," FS gluons, there are ",\
                #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
                #      len(matrix_element.get_all_wavefunctions()),\
                #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
                #      ' diagrams'

                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))
                self.assertEqual(len(matrix_element.get_all_amplitudes()),
                                 goal_amplitudes[ngluon][iflow])
                self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                                 goal_wavefunctions[ngluon][iflow])

                # Test JAMP (color amplitude) output
                #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

            #goal_wavefunctions.append(wfs)
            #goal_amplitudes.append(amps)
            #goal_ndiags.append(diags)

            #print goal_nflows
            #print goal_wavefunctions
            #print goal_amplitudes
            #print goal_ndiags


    def test_matrix_element_uux_ddxssxng(self):
        """Test color flow matrix element for uu~>dd~ss~ng
        """

        goal_wavefunctions =  [[12, 13, 13, 15, 12, 13],
                               [22, 22, 24, 18, 22, 24, 24, 24, 22, 22, 18, 22, 22, 22, 18, 24, 22, 24],
                               [43, 42, 41, 47, 36, 28, 44, 42, 52, 36, 43, 47, 47, 52, 47, 37, 36, 38, 36, 36, 36, 39, 28, 38, 43, 44, 43, 27, 33, 44, 42, 42, 33, 49, 41, 44]]
        goal_amplitudes =  [[4, 4, 4, 6, 4, 4],
                            [11, 11, 10, 4, 11, 10, 10, 10, 10, 10, 4, 10, 11, 11, 4, 10, 11, 10],
                            [22, 24, 24, 19, 7, 3, 22, 24, 19, 7, 22, 19, 19, 19, 19, 14, 13, 13, 7, 7, 13, 10, 3, 13, 22, 22, 22, 4, 10, 22, 24, 24, 10, 22, 24, 22]]
        goal_ndiags =  [[4, 4, 4, 6, 4, 4],
                        [16, 16, 14, 8, 16, 14, 14, 14, 14, 14, 8, 14, 16, 16, 8, 14, 16, 14],
                        [63, 69, 63, 50, 28, 20, 69, 69, 56, 28, 63, 50, 50, 56, 50, 38, 32, 38, 28, 28, 32, 32, 20, 38, 63, 69, 63, 20, 28, 50, 69, 69, 28, 56, 63, 50]]

        goal_nflows = [6, 18, 36]

        #goal_wavefunctions = []
        #goal_amplitudes = []
        #goal_ndiags = []
        
        # Test 0-2 gluons in the final state
        for ngluon in range(0, 1):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1}),
                base_objects.Leg({'id':3}),
                base_objects.Leg({'id':-3})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            #goal_nflows.append(len(self.myamplitude.get('color_flows')))
            #diags=[]
            #amps=[]
            #wfs=[]
            for iflow, mycolorflow in \
                enumerate(self.myamplitude.get('color_flows')):
                #print mycolorflow.nice_string()
                
                matrix_element = color_ordered_amplitudes.COHelasFlow(\
                    mycolorflow, gen_color=False, optimization = 3)

                #print "\n".join(\
                #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
                #    get_matrix_element_calls(matrix_element))
                # Test JAMP (color amplitude) output
                #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

                # Check that we get back the correct number of diagrams
                # when we do get_base_amplitude. This is a very powerful
                # check that the BG currents are correct.
                diagrams = mycolorflow.get('diagrams')
                base_amplitude = matrix_element.get_base_amplitude().get('diagrams')
                diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                                for d in diagrams]
                base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                             for d in base_amplitude]
                base_tag_copy = copy.copy(base_tags)
                
                for tag in diagram_tags:
                    self.assertTrue(tag in base_tags)
                    base_tags.remove(tag)
                    
                self.assertEqual(base_tags, [])

                # Test that wavefunctions that have been combined are
                # not in any amplitudes
                comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                                for w in matrix_element.get_all_wavefunctions()\
                                if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
                amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                               for amp in matrix_element.get_all_amplitudes()],\
                              []))
                self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

                # Test that all wavefunctions are used
                wf_numbers = [w.get('number') for w in \
                              matrix_element.get_all_wavefunctions()]
                mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                               for a in matrix_element.get_all_amplitudes() + \
                                         matrix_element.get_all_wavefunctions()],\
                              []))
                self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

                for d in matrix_element.get('diagrams'):
                    self.assertTrue(len(d.get('amplitudes')) > 0)

                #print "For ",ngluon," FS gluons, ",iflow+1,"th flow, there are ",\
                #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
                #      len(matrix_element.get_all_wavefunctions()),\
                #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
                #      ' diagrams'

                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))

                self.assertEqual(len(matrix_element.get_all_amplitudes()),
                                 goal_amplitudes[ngluon][iflow])
                self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                                 goal_wavefunctions[ngluon][iflow])

            #goal_wavefunctions.append(wfs)
            #goal_amplitudes.append(amps)
            #goal_ndiags.append(diags)

            #print 'goal_wavefunctions = ',goal_wavefunctions
            #print 'goal_amplitudes = ',goal_amplitudes
            #print 'goal_ndiags = ',goal_ndiags
                              
    def test_matrix_element_uux_uuxuuxng(self):
        """Test color flow matrix element for uu~>uu~uu~ng
        """
        goal_wavefunctions =  [[18, 27, 15],
                               [35, 48, 22, 35, 48, 22, 35, 48, 22],
                               [69, 87, 37, 74, 97, 36, 67, 85, 38, 74, 97, 36, 72, 97, 39, 67, 85, 38]]
        goal_amplitudes =  [[8, 12, 6],
                            [22, 24, 10, 22, 24, 10, 22, 24, 10],
                            [44, 42, 14, 46, 36, 13, 46, 44, 13, 46, 36, 13, 48, 36, 10, 46, 44, 13]]
        goal_ndiags =  [[8, 12, 6],
                        [32, 36, 14, 32, 36, 14, 32, 36, 14],
                        [126, 120, 38, 138, 112, 32, 126, 120, 38, 138, 112, 32, 138, 112, 32, 126, 120, 38]]
        goal_nflows =  [3, 9, 18]

        #goal_wavefunctions = []
        #goal_amplitudes = []
        #goal_ndiags = []
        #goal_nflows = []
        
        # Test 0-2 gluons in the final state
        for ngluon in range(0, 1):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':2}),
                base_objects.Leg({'id':-2}),
                base_objects.Leg({'id':2}),
                base_objects.Leg({'id':-2})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluon])
            #goal_nflows.append(len(self.myamplitude.get('color_flows')))

            #diags=[]
            #amps=[]
            #wfs=[]
            for iflow, mycolorflow in \
                enumerate(self.myamplitude.get('color_flows')):
                #print mycolorflow.nice_string()
                
                matrix_element = color_ordered_amplitudes.COHelasFlow(\
                    mycolorflow, gen_color=False, optimization = 3)

                #print "\n".join(\
                #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
                #    get_matrix_element_calls(matrix_element))
                # Test JAMP (color amplitude) output
                #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

                # Check that we get back the correct number of diagrams
                # when we do get_base_amplitude. This is a very powerful
                # check that the BG currents are correct.
                diagrams = mycolorflow.get('diagrams')
                base_amplitude = matrix_element.get_base_amplitude().get('diagrams')
                diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                                for d in diagrams]
                base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                             for d in base_amplitude]
                base_tag_copy = copy.copy(base_tags)
                
                for tag in diagram_tags:
                    self.assertTrue(tag in base_tags)
                    base_tags.remove(tag)
                    
                self.assertEqual(base_tags, [])

                # Test that wavefunctions that have been combined are
                # not in any amplitudes
                comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                                for w in matrix_element.get_all_wavefunctions()\
                                if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
                amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                               for amp in matrix_element.get_all_amplitudes()],\
                              []))
                self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

                # Test that all wavefunctions are used
                wf_numbers = [w.get('number') for w in \
                              matrix_element.get_all_wavefunctions()]
                mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                               for a in matrix_element.get_all_amplitudes() + \
                                         matrix_element.get_all_wavefunctions()],\
                              []))
                self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

                for d in matrix_element.get('diagrams'):
                    self.assertTrue(len(d.get('amplitudes')) > 0)

                #print "For ",ngluon," FS gluons, ",iflow+1,"th flow, there are ",\
                #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
                #      len(matrix_element.get_all_wavefunctions()),\
                #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
                #      ' diagrams'

                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))

                self.assertEqual(len(matrix_element.get_all_amplitudes()),
                                 goal_amplitudes[ngluon][iflow])
                self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                                 goal_wavefunctions[ngluon][iflow])

            #goal_wavefunctions.append(wfs)
            #goal_amplitudes.append(amps)
            #goal_ndiags.append(diags)

            #print 'goal_wavefunctions = ',goal_wavefunctions
            #print 'goal_amplitudes = ',goal_amplitudes
            #print 'goal_ndiags = ',goal_ndiags
            #print 'goal_nflows = ',goal_nflows
                              
    def test_matrix_element_uux_uuxddxng(self):
        """Test color flow matrix element for uu~>uu~dd~ng
        """

        goal_wavefunctions =  [[12, 13, 13, 15, 12, 13],
                               [22, 24, 22, 24, 24, 22, 18, 22, 22, 18, 24, 22, 22, 18, 22, 24, 22, 24],
                               [43, 47, 42, 52, 41, 47, 47, 37, 36, 36, 28, 38, 44, 36, 42, 36, 52, 36, 36, 39, 43, 28, 47, 38, 43, 27, 44, 33, 43, 44, 42, 33, 42, 49, 41, 44]]
        goal_amplitudes =  [[4, 4, 4, 6, 4, 4],
                            [11, 10, 11, 10, 10, 10, 4, 10, 11, 4, 10, 10, 11, 4, 11, 10, 11, 10],
                            [22, 19, 24, 19, 24, 19, 19, 14, 7, 13, 3, 13, 22, 7, 24, 7, 19, 13, 7, 10, 22, 3, 19, 13, 22, 4, 22, 10, 22, 22, 24, 10, 24, 22, 24, 22]]
        goal_ndiags =  [[4, 4, 4, 6, 4, 4],
                        [16, 14, 16, 14, 14, 14, 8, 14, 16, 8, 14, 14, 16, 8, 16, 14, 16, 14],
                        [63, 50, 69, 56, 63, 50, 50, 38, 28, 32, 20, 38, 69, 28, 69, 28, 56, 32, 28, 32, 63, 20, 50, 38, 63, 20, 69, 28, 63, 50, 69, 28, 69, 56, 63, 50]]
        goal_nflows =  [6, 18, 36]
        #goal_wavefunctions = []
        #goal_amplitudes = []
        #goal_ndiags = []
        #goal_nflows = []
        
        # Test 0-2 gluons in the final state
        for ngluon in range(0, 1):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':2}),
                base_objects.Leg({'id':-2}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            self.myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            self.assertEqual(len(self.myamplitude.get('color_flows')),
                             goal_nflows[ngluon])
            #goal_nflows.append(len(self.myamplitude.get('color_flows')))

            #diags=[]
            #amps=[]
            #wfs=[]
            for iflow, mycolorflow in \
                enumerate(self.myamplitude.get('color_flows')):
                #print mycolorflow.nice_string()
                
                matrix_element = color_ordered_amplitudes.COHelasFlow(\
                    mycolorflow, gen_color=False, optimization = 3)

                #print "\n".join(\
                #    color_ordered_export_v4.COFortranUFOHelasCallWriter(self.mymodel).\
                #    get_matrix_element_calls(matrix_element))
                # Test JAMP (color amplitude) output
                #print '\n'.join(export_v4.get_JAMP_lines(matrix_element))

                # Check that we get back the correct number of diagrams
                # when we do get_base_amplitude. This is a very powerful
                # check that the BG currents are correct.
                diagrams = mycolorflow.get('diagrams')
                base_amplitude = matrix_element.get_base_amplitude().get('diagrams')
                diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                                for d in diagrams]
                base_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                             for d in base_amplitude]
                base_tag_copy = copy.copy(base_tags)
                
                for tag in diagram_tags:
                    self.assertTrue(tag in base_tags)
                    base_tags.remove(tag)
                    
                self.assertEqual(base_tags, [])

                # Test that wavefunctions that have been combined are
                # not in any amplitudes
                comb_wfs = set(sum([[m.get('number') for m in w.get('mothers')] \
                                for w in matrix_element.get_all_wavefunctions()\
                                if isinstance(w, color_ordered_amplitudes.BGHelasCurrent)], []))
                amp_wfs = set(sum([[m.get('number') for m in amp.get('mothers')] \
                               for amp in matrix_element.get_all_amplitudes()],\
                              []))
                self.assertFalse(any([num in amp_wfs for num in comb_wfs]))

                # Test that all wavefunctions are used
                wf_numbers = [w.get('number') for w in \
                              matrix_element.get_all_wavefunctions()]
                mother_numbers = set(sum([[m.get('number') for m in a.get('mothers')] \
                               for a in matrix_element.get_all_amplitudes() + \
                                         matrix_element.get_all_wavefunctions()],\
                              []))
                self.assertTrue(all([num in mother_numbers for num in wf_numbers]))

                for d in matrix_element.get('diagrams'):
                    self.assertTrue(len(d.get('amplitudes')) > 0)

                #print "For ",ngluon," FS gluons, ",iflow+1,"th flow, there are ",\
                #      len(matrix_element.get_all_amplitudes()),' amplitudes and ',\
                #      len(matrix_element.get_all_wavefunctions()),\
                #      ' wavefunctions for ', len(mycolorflow.get('diagrams')),\
                #      ' diagrams'

                #diags.append(len(mycolorflow.get('diagrams')))
                #wfs.append(len(matrix_element.get_all_wavefunctions()))
                #amps.append(len(matrix_element.get_all_amplitudes()))

                self.assertEqual(len(matrix_element.get_all_amplitudes()),
                                 goal_amplitudes[ngluon][iflow])
                self.assertEqual(len(matrix_element.get_all_wavefunctions()),
                                 goal_wavefunctions[ngluon][iflow])

            #goal_wavefunctions.append(wfs)
            #goal_amplitudes.append(amps)
            #goal_ndiags.append(diags)

            #print 'goal_wavefunctions = ',goal_wavefunctions
            #print 'goal_amplitudes = ',goal_amplitudes
            #print 'goal_ndiags = ',goal_ndiags
            #print 'goal_nflows = ',goal_nflows
                              
    def test_color_matrix_uux_ddxng(self):
        """Test color flow matrix elements for uu~>dd~ng
        """

        # Test 0-3 gluons in the final state
        for ngluon in range(0, 3):

            # Create the non-optimized amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED': 0},
                                           'model':self.mymodel})

            myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)

            matrix_element1 = color_ordered_amplitudes.COHelasMatrixElement(\
                myamplitude, gen_color=3, optimization=1)

            matrix_element3 = color_ordered_amplitudes.COHelasMatrixElement(\
                myamplitude, gen_color=3, optimization=3)

            for i,(cf1,cf3) in enumerate(zip(matrix_element1.get('color_flows'),
                                             matrix_element3.get('color_flows'))):
                self.assertEqual(cf1.get('color_string'),
                                 cf3.get('color_string'))
            self.assertEqual(matrix_element1.get('color_basis'),
                             matrix_element3.get('color_basis'))


#===============================================================================
# TestDiagramTag
#===============================================================================
class TestDiagramTag(unittest.TestCase):
    """Test class for the DiagramTag class"""


    def setUp(self):

        self.mypartlist = base_objects.ParticleList()
        self.myinterlist = base_objects.InteractionList()
        self.mymodel = base_objects.Model()

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
    
    def test_order_diagram_tag_uux_nglue(self):
        """Test the OrderDiagramTag for u u~ > n g
        """

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range(2, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon, 'QED': 0},
                                           'model':self.mymodel})

            myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            diagrams = myamplitude.get('color_flows')[0].get('diagrams')
            diagram_tags = [color_ordered_amplitudes.OrderDiagramTag(d) \
                            for d in diagrams]

            #print myamplitude.get('process').nice_string()
            
            for i,(d,dtag) in enumerate(zip(diagrams, diagram_tags)):
                #print '%3r: ' % (i+1),d.nice_string()
                #print 'new: ',dtag.diagram_from_tag(self.mymodel).nice_string()
                # Check that the resulting diagram is recreated in the same way
                # from the diagram tag (by checking the diagram tag)
                self.assertEqual(dtag,
                                 color_ordered_amplitudes.OrderDiagramTag(\
                                     dtag.diagram_from_tag(self.mymodel)))
