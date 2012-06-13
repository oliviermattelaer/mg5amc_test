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
import madgraph.iolibs.save_load_object as save_load_object

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
        self.mymodel.set('order_hierarchy', {'QCD': 1, 'QED': 2})

        self.ref_dict_to0 = self.myinterlist.generate_ref_dict()[0]
        self.ref_dict_to1 = self.myinterlist.generate_ref_dict()[1]

    def test_color_ordered_model(self):
        """Test the ColorOrderedModel based on this model"""
        co_model = color_ordered_amplitudes.ColorOrderedModel(self.mymodel)
        self.assertEqual(co_model.get('order_hierarchy'),
                         {'QCD': 1, 'singlet_QCD': 1, 'QED': 2})
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
             'couplings': {(0, 0): 'GQQ'},
             'orders': {'singlet_QCD': 1}})

        self.assertTrue(singlet_interaction1 in co_model.get('interactions'))

        singlet_interaction2 = base_objects.Interaction(\
            {'id': 13,
             'particles': base_objects.ParticleList([co_model.get_particle(1),
                                                     co_model.get_particle(-1),
                                                     singlet]),
             'color': [color_string],
             'lorentz': ['L1'],
             'couplings': {(0, 0): 'GQQ'},
             'orders': {'singlet_QCD': 1}})
        self.assertTrue(singlet_interaction2 in co_model.get('interactions'))
        singlet_interaction3 = base_objects.Interaction(\
            {'id': 14,
             'particles': base_objects.ParticleList([co_model.get_particle(3),
                                                     co_model.get_particle(-3),
                                                     singlet]),
             'color': [color_string],
             'lorentz': ['L1'],
             'couplings': {(0, 0): 'GQQ'},
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
        goal_ndiags = [[1, 1], [3, 3, 2, 2], [10, 11, 10, 5, 4, 5], [37, 41, 41, 37, 16, 10, 10, 16]]
        goal_nflows = [2, 4, 6, 8]
        goal_nperms = [[2] * 2, [2] * 4, [4] * 6, [12] * 8]
        
        #goal_ndiags = []

        for ngluons in range(0, 4):

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
            #diags = []
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                #diags.append(len(flow.get('diagrams')))
                self.assertEqual(len(flow.get('diagrams')),
                             goal_ndiags[ngluons][iflow])
                self.assertEqual(len(flow.get('permutations')),
                                 goal_nperms[ngluons][iflow])
                # Check if all perms are equal
                perms = flow.get('permutations')
                perms = [base_objects.reorder_permutation(perm, flow.get('permutations')[0]) for \
                       perm in flow.get('permutations')]
                self.assertEqual(perms, perms_0)
            #goal_ndiags.append(diags)
        #print 'goal_ndiags = ',goal_ndiags

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
        goal_ndiags = [[4, 4, 4, 4, 4, 6], 
                       [16, 16, 16, 16, 16, 16, 14, 8, 14, 14, 14, 8, 8, 14, 14, 14, 14, 14], 
                       [63, 69, 63, 69, 69, 63, 63, 69, 63, 69, 69, 63, 50, 28, 20, 56, 28, 50, 50, 56, 50, 28, 28, 20, 20, 28, 50, 28, 56, 50, 38, 32, 38, 32, 32, 38]]
        
        goal_nflows = [6, 18, 36]
        #goal_nflows = []
        #goal_ndiags = []
        for ngluons in range(0, 2):

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
        goal_ndiags =  [[4, 4, 4, 4, 4, 6], 
                        [16, 16, 16, 16, 16, 16, 14, 14, 14, 8, 8, 14, 8, 14, 14, 14, 14, 14], 
                        [63, 69, 63, 69, 69, 63, 63, 69, 63, 69, 69, 63, 50, 56, 50, 50, 28, 20, 28, 28, 56, 28, 20, 50, 20, 28, 50, 28, 56, 50, 38, 32, 38, 32, 32, 38]]
        
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
        goal_ndiags = [[4, 4], [14, 14, 10, 10], [50, 56, 50, 28, 24, 28]]
        goal_nflows = [2, 4, 6]

        #goal_ndiags = []

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
            #diags = []
            for iflow, flow in enumerate(self.myamplitude.get('color_flows')):
                #diags.append(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')))
                self.assertEqual(len(self.myamplitude.get('color_flows')[iflow].get('diagrams')),
                             goal_ndiags[ngluons][iflow])
                
            #goal_ndiags.append(diags)
        #print "goal_ndiags = ",goal_ndiags

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

    def test_periferal_diagrams_gluons(self):
        """Test periferal diagrams for gg>ng"""

        goal_ndiags = [12, 66, 180, 990]

        # Time for 6 gluons: 5 min

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range (3, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(\
                                 include_all_t = False, tch_depth = 1,
                                 identify_depth = 10)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "allperiferal%i.eps" % ngluon,
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 3])
        # print goal_ndiags

    def test_periferal_diagrams_gluons_tch_10_id_1(self):
        """Test periferal diagrams for gg>ng"""

        goal_ndiags = [6, 12, 20, 30]

        # Time for 6 gluons: 5 min

        # Test 3, 4 and 5 gluons in the final state
        for ngluon in range (3, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(\
                            include_all_t = True, tch_depth = 10,
                            identify_depth = 1)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "allperiferal101%i.eps" % ngluon, 
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 3])
        # print goal_ndiags

    def test_periferal_diagrams_gluons_tch_2_id_1(self):
        """Test periferal diagrams for gg>ng"""

        goal_ndiags = [6, 24, 80, 300]

        # Time for 6 gluons: 5 min

        # Test 3, 4 and 5 gluons in the final state
        for ngluon in range (3, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(\
                            include_all_t = True, tch_depth = 2,
                            identify_depth = 1)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "allperiferal21%i.eps" % ngluon, 
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 3])
        # print goal_ndiags

    def test_periferal_diagrams_gluons_tch_2_id_2(self):
        """Test periferal diagrams for gg>ng"""

        goal_ndiags = [6, 36, 240, 1530]

        # Time for 6 gluons: 5 min

        # Test 3, 4 and 5 gluons in the final state
        for ngluon in range (3, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(\
                            include_all_t = True, tch_depth = 2,
                            identify_depth = 2)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "allperiferal22%i.eps" % ngluon, 
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 3])
        # print goal_ndiags

    def test_periferal_diagrams_gluons_tch_3_id_2(self):
        """Test periferal diagrams for gg>ng"""

        goal_ndiags = [6, 24, 120, 720]

        # Time for 6 gluons: 5 min

        # Test 3, 4 and 5 gluons in the final state
        for ngluon in range (3, 6):

            # Create the amplitude
            myleglist = base_objects.LegList([base_objects.Leg({'id':21,
                                              'state':False})] * 2)

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(\
                            include_all_t = True, tch_depth = 3,
                            identify_depth = 2)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "allperiferal22%i.eps" % ngluon, 
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 3])
        # print goal_ndiags

    def test_periferal_diagrams_uux_ddxng(self):
        """Test periferal diagrams for uu~>dd~+ng"""

        goal_ndiags = []

        # Time for 6 gluons: 5 min

        # Test 2, 3, 4 and 5 gluons in the final state
        for ngluon in range (1, 4):

            # Create the amplitude
            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':2, 'state':False}),
                base_objects.Leg({'id':-2, 'state':False}),
                base_objects.Leg({'id':1}),
                base_objects.Leg({'id':-1})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * ngluon)

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QED':0},
                                           'model':self.mymodel})

            self.myamplitude.set('process', myproc)

            self.myamplitude.generate_diagrams()

            diagrams, tch_depth = \
                      self.myamplitude.get_periferal_diagrams_from_flows(identify_depth=1)

            print diagrams.nice_string()
            plot = draw.MultiEpsDiagramDrawer(diagrams,
                                              "uuxddxperiferal%i.eps" % ngluon,
                                              model=self.mymodel)
            plot.draw()
            goal_ndiags.append(len(diagrams))
            self.assertEqual(len(diagrams), goal_ndiags[ngluon - 1])
        print goal_ndiags

    def test_periferal_diagrams_gu_wpgd(self):
        """Test periferal diagrams for uu~>dd~+ng"""

        # Time for 6 gluons: 5 min

        # Create the amplitude
        myleglist = base_objects.LegList([\
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':2, 'state':False}),
            base_objects.Leg({'id':24}),
            base_objects.Leg({'id':21}),
            base_objects.Leg({'id':1})])

        myproc = base_objects.Process({'legs':myleglist,
                                       'orders':{'WEIGHTED':4},
                                       'model':self.mymodel})

        self.myamplitude.set('process', myproc)

        self.myamplitude.generate_diagrams()

        for i,cf in enumerate(self.myamplitude.get('color_flows')):
            plot = draw.MultiEpsDiagramDrawer(cf.get('diagrams'),
                                              "test1_%d.eps" % (i+1),
                                              model=self.mymodel)
            plot.draw()

        diagrams, tch_depth = \
                  self.myamplitude.get_periferal_diagrams_from_flows(True,
                                                                     10,
                                                                     1)

        print diagrams.nice_string()
        self.assertEqual(len(diagrams), 3)

    def test_periferal_diagrams_gdx_wpgux(self):
        """Test periferal diagrams for uu~>dd~+ng"""

        # Time for 6 gluons: 5 min

        # Create the amplitude
        myleglist = base_objects.LegList([\
            base_objects.Leg({'id':21, 'state':False}),
            base_objects.Leg({'id':-1, 'state':False}),
            base_objects.Leg({'id':24}),
            base_objects.Leg({'id':21}),
            base_objects.Leg({'id':-2})])

        myproc = base_objects.Process({'legs':myleglist,
                                       'orders':{'WEIGHTED':4},
                                       'model':self.mymodel})

        self.myamplitude.set('process', myproc)

        self.myamplitude.generate_diagrams()

        for i,cf in enumerate(self.myamplitude.get('color_flows')):
            plot = draw.MultiEpsDiagramDrawer(cf.get('diagrams'),
                                              "test2_%i.eps" % (i+1),
                                              model=self.mymodel)
            plot.draw()

        diagrams, tch_depth = \
                  self.myamplitude.get_periferal_diagrams_from_flows(True,
                                                                     10,
                                                                     1)

        print diagrams.nice_string()
        self.assertEqual(len(diagrams), 3)

#===============================================================================
# TestOrderDiagramTag
#===============================================================================
class TestOrderDiagramTag(unittest.TestCase):
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

#===============================================================================
# TestPeriferalDiagramTag
#===============================================================================
class TestPeriferalDiagramTag(unittest.TestCase):
    """Test class for the PeriferalDiagramTag class"""


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
    
    def test_periferal_diagram_tag_gg_nglue(self):
        """Test the PeriferalDiagramTag for g g > n g
        """

        goal_periferals = [[1, 2, 4, 5, 7],
                           [3, 13],
                           [5, 6, 8, 9, 19, 20, 43, 44, 46, 47, 57, 58, 81, 82]]
        maxgluons = 4

        # Test 3, 4 and 5 gluons in the final state
        for ngluon in range(0, maxgluons-2):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':21, 'state':False}),
                base_objects.Leg({'id':21, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * (ngluon + 3))

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon+3, 'QED': 0},
                                           'model':self.mymodel})

            myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            diagrams = myamplitude.get('color_flows')[0].get('diagrams')
            diagram_tags = [color_ordered_amplitudes.PeriferalDiagramTag(d) \
                            for d in diagrams]
            # plot = draw.MultiEpsDiagramDrawer(diagrams,
            #                                   "alldiags%i.eps" % ngluon,
            #                                   model=self.mymodel)
            # plot.draw()

            periferals = []
            # periferal_diagrams = base_objects.DiagramList()
            for i, (d, dtag) in enumerate(zip(diagrams, diagram_tags)):
                # print d.nice_string()
                # print i+1,dtag.check_periferal_diagram(self.mymodel)
                check_res = dtag.check_periferal_diagram(self.mymodel)
                if check_res:
                    periferals.append(i+1)
            #         periferal_diagrams.append(d)
            # plot = draw.MultiEpsDiagramDrawer(periferal_diagrams,
            #                                   "periferals%i.eps" % ngluon,
            #                                   model=self.mymodel)
            # plot.draw()
            
            # goal_periferals.append(periferals)
            self.assertEqual(periferals, goal_periferals[ngluon])

        # print 'goal_periferals = ',goal_periferals

    def test_periferal_diagram_tag_gg_nglue_order2(self):
        """Test the PeriferalDiagramTag for g g > n g with order=2
        """

        goal_periferals = [[1, 2, 4, 5, 7],
                           [2, 3, 5, 6, 8, 12, 13, 15, 16, 18, 22, 24, 25, 28],
                           [5, 6, 8, 9, 19, 20, 43, 44, 46, 47, 57, 58, 81, 82],
                           [9, 14, 15, 19, 20, 22, 30, 55, 57, 60, 64, 65, 70,
                            98, 163, 168, 169, 173, 174, 176, 184, 209, 211,
                            214, 218, 219, 224, 252, 315, 321, 322, 325, 349,
                            352, 353, 420]]

        goal_pass = [[4, 5, 7],
                     [12, 13, 16, 25, 28],
                     [44, 47, 58, 82],
                     [169, 173, 211, 352, 420]]

        maxgluons = 4

        # Test 3,4 and 5 gluons in the final state
        for ngluon in range(0, maxgluons-2):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':21, 'state':False}),
                base_objects.Leg({'id':21, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * (ngluon + 3))

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon+3, 'QED': 0},
                                           'model':self.mymodel})

            myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            diagrams = myamplitude.get('color_flows')[0].get('diagrams')
            diagram_tags = [color_ordered_amplitudes.PeriferalDiagramTag(d) \
                            for d in diagrams]
            # plot = draw.MultiEpsDiagramDrawer(diagrams,
            #                                   "alldiags%i.eps" % ngluon,
            #                                   model=self.mymodel)
            # plot.draw()

            periferals = []
            pass_restrictions = []
            periferal_diagrams = base_objects.DiagramList()
            for i, (d, dtag) in enumerate(zip(diagrams, diagram_tags)):
                check_res = dtag.check_periferal_diagram(self.mymodel,
                                                            order = 2)
                # print d.nice_string()
                # print i+1,check_res
                if not check_res: continue
                periferals.append(i+1)
                if dtag.pass_restrictions(self.mymodel):
                    pass_restrictions.append(i+1)
            #         periferal_diagrams.append(d)
            # plot = draw.MultiEpsDiagramDrawer(periferal_diagrams,
            #                                   "pass%i.eps" % ngluon,
            #                                   model=self.mymodel)
            # plot.draw()
            
            # goal_periferals.append(periferals)
            # goal_pass.append(pass_restrictions)
            self.assertEqual(periferals, goal_periferals[ngluon])

        # print 'goal_periferals = ',goal_periferals
        # print 'goal_pass = ',goal_pass

    def test_periferal_diagram_tag_gg_nglue_order2_include_all_t(self):
        """Test the PeriferalDiagramTag for g g > n g with order=2 including t
        """

        goal_periferals = []

        goal_pass = []

        maxgluons = 6

        # Test 3,4 and 5 gluons in the final state
        for ngluon in range(0, maxgluons-2):

            # Create the amplitude
            myleglist = base_objects.LegList([\
                base_objects.Leg({'id':21, 'state':False}),
                base_objects.Leg({'id':21, 'state':False})])

            myleglist.extend([base_objects.Leg({'id':21,
                                              'state':True})] * (ngluon + 3))

            myproc = base_objects.Process({'legs':myleglist,
                                           'orders':{'QCD':ngluon+3, 'QED': 0},
                                           'model':self.mymodel})

            myamplitude = color_ordered_amplitudes.ColorOrderedAmplitude(myproc)
            diagrams = myamplitude.get('color_flows')[0].get('diagrams')
            diagram_tags = [color_ordered_amplitudes.PeriferalDiagramTag(d) \
                            for d in diagrams]
            # plot = draw.MultiEpsDiagramDrawer(diagrams,
            #                                   "alldiags%i.eps" % ngluon,
            #                                   model=self.mymodel)
            # plot.draw()

            periferals = []
            pass_restrictions = []
            periferal_diagrams = base_objects.DiagramList()
            pass_diagrams = base_objects.DiagramList()
            for i, (d, dtag) in enumerate(zip(diagrams, diagram_tags)):
                check_res = dtag.check_periferal_diagram(self.mymodel,
                                                         order = 2,
                                                         include_all_t = True)
                # print d.nice_string()
                # print i+1,check_res
                if not check_res: continue
                periferals.append(i+1)
                periferal_diagrams.append(d)
                if dtag.pass_restrictions(self.mymodel):
                    pass_restrictions.append(i+1)
                    pass_diagrams.append(d)
            plot = draw.MultiEpsDiagramDrawer(periferal_diagrams,
                                              "periferal%i.eps" % ngluon,
                                              model=self.mymodel)
            plot.draw()
            plot = draw.MultiEpsDiagramDrawer(pass_diagrams,
                                              "pass%i.eps" % ngluon,
                                              model=self.mymodel)
            plot.draw()
            
            goal_periferals.append(periferals)
            goal_pass.append(pass_restrictions)
            self.assertEqual(periferals, goal_periferals[ngluon])

        print 'goal_periferals = ',goal_periferals
        print 'goal_pass = ',goal_pass

