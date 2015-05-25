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
 
"""Classes for generation of color-ordered diagrams.

ColorOrderedAmplitude keeps track of all color flows, and
ColorOrderedFlow performs the diagram generation. ColorOrderedLeg has
extra needed flags. See documentation for the functions
get_combined_legs and get_combined_vertices for the algorithm used to
generate only color-ordered diagrams.

ColorOrderedModel adds color singlets coupling only to triplets for
each color octet in the model.
"""
import array
import copy
import fractions
import itertools
import logging

import models.import_ufo as import_ufo
import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.color_amp as color_amp
import madgraph.core.diagram_generation as diagram_generation
from madgraph import MadGraph5Error
import madgraph.various.misc as misc

logger = logging.getLogger('madgraph.color_ordered_amplitudes')

#===============================================================================
# DiagramTag class to order diagram vertices in most optimal way
#===============================================================================

class OrderDiagramTagChainLink(diagram_generation.DiagramTagChainLink):
    """Chain link class for OrderDiagramTag with different ordering"""

    def get_flat_links(self):
        """Return a sorted, flattened list of links"""

        if self.end_link:
            return [self.links[0]]
        return sorted(sum([l.get_flat_links() for l in self.links],[]))

    def __lt__(self, other):
        """Compare self with other in the order:
        1. depth 2. measure of flat links
        3. len(links) 4. vertex id"""

        if self == other:
            return False

        if self.depth != other.depth:
            return self.depth < other.depth

        return self.get_flat_links() < other.get_flat_links()

        if len(self.links) != len(other.links):
            return len(self.links) < len(other.links)

        if self.vertex_id != other.vertex_id:
            return self.vertex_id < other.vertex_id

class OrderDiagramTag(diagram_generation.DiagramTag):
    """DiagramTag daughter class to order the vertices in all diagrams
    to ensure that all diagrams are correctly combined into BG
    currents for color ordered amplitudes. For regular amplitudes,
    this gives the most optimal ordering for speed of matrix element
    calculation."""

    link_class = OrderDiagramTagChainLink

    @staticmethod
    def link_from_leg(leg, model):
        """Returns the info needed to recreate a leg."""

        # Include also onshell, since this specifies forbidden s-channel
        return [((leg.get('number'), leg.get('id'), leg.get('state'),
                  leg.get('onshell')),
                  leg.get('number'))]
        
    @staticmethod
    def leg_from_link(link):
        """Return a leg from a link"""

        if link.end_link:
            # This is an external leg, info in links
            return base_objects.Leg({'number':link.links[0][0][0],
                                     'id':link.links[0][0][1],
                                     'state':link.links[0][0][2],
                                     'onshell':link.links[0][0][3]})

        # This shouldn't happen
        assert False

#===============================================================================
# DiagramTag class to select and manipulate periferal diagrams for
# color ordered phase space integration
#===============================================================================

class PeriferalDiagramTagChainLink(OrderDiagramTagChainLink):
    """Chain link class for PeriferalDiagramTag with special functions"""

    def check_periferal_legs(self, model, order, nallowed, 
                             include_all_t = False, amp_link = False):
        """Check if this is a periferal subdiagram (to the order given).
        amp_link is True for a link corresponding to an amplitude,
        False for a link corresponding to a wavefunction.

        Periferal means that at most one (massless) external leg
        attaches to a vertex with depth > order. Decay-type vertices
        (where a non-zero mass parameter appears only once) are not
        included in this condition. Returns the number of external
        legs connecting to a forbidden vertex plus a sorted list of
        tuples of the external numbers of the vertices,
        e.g. [(1,3),(2,4),5] for a t-channel 5-particle diagram, or an
        empty list if the diagram is not periferal."""

        # Algorithm: Go through the links and collect external
        # particles. If more than 1 external particle couples to a
        # vertex with depth > order, return []

        # If this is an external link, return leg number (as an array)
        if self.end_link:
            return 0, [self.links[0][1]]

        # Remove any 4-particle vertices
        if len(self.links) > (3 if amp_link else 2):
            return 0, []

        # If this vertex has depth 1, return a tuple of sorted external legs
        if self.depth == 1:
            return 0, [tuple(sorted([l.links[0][1] for l in self.links]))]

        # Otherwise, go through mother links
        leglist = []
        nexternal = 0
        for link in self.links:
            n, legs = link.check_periferal_legs(model, order, nallowed,
                                                include_all_t)
            # If any link returns an empty list, we already have a
            # non-periferal subdiagram
            if not legs:
                return nexternal, []
            nexternal += n
            if nexternal > nallowed and not include_all_t:
                return nexternal, []
            leglist.extend(legs)

        # If the depth of this link is > order, check if we have
        # external legs

        depth = self.depth
        if amp_link:
            # For the last link, we need to use the depth of the
            # shortest non-external leg
            depth = min([l.depth for l in self.links if not l.end_link]) + 1

        decay_vertex = False

        if depth > order:
            # A decay vertex has one non-zero mass parameter coming
            # only once
            masses = [p.get('mass') for p in \
                      model.get_interaction(self.vertex_id[0]).get('particles') \
                      if p.get('mass').lower() != 'zero']
            for mass in set(masses):
                if len([m for m in masses if m == mass]) == 1:
                    decay_vertex = True
                    break
            # Now check if there is any external leg attached
            if not decay_vertex:
                n = len([l for l in self.links if l.end_link])
                nexternal += n
                # Fail if external legs and depth is > order + 1, or
                # if nexternal > nallowed
                if nexternal > nallowed or (n > 0 and depth > order + 1):
                    if not include_all_t or n != 1:
                        return nexternal, []

        # For amp vertex or deep vertices that are not decay vertices,
        # return the sorted leglist
        if amp_link or depth > order and not decay_vertex:
            return nexternal, sorted(leglist)
        # For vertices with depth <= order and decay vertices, return a tuple
        return nexternal, [tuple(sorted(leglist))]

    def pass_restrictions(self, model, amp_link = False, tch_depth = 1):
        """Check if this subdiagram passes the following
        restrictions:

        - Maximum two external legs connecting to any s-channel
          propagator
        - No s-channel propagators connecting within tch_depth of
          initial state legs (depth 0 means that initial state legs
          can connect, 1 means s-channel can connect directly to initial
          state leg, 2 that initial-state leg can connect only to
          final-state external leg, etc)
        - tch_depth == 10 is short for only allowing t-channel diagrams
        """

        # If tch_depth == 0, we allow any diagrams
        if not tch_depth:
            return True
            
        # First check if any daughter links fail
        if any([not l.pass_restrictions(model, False, tch_depth) \
                for l in self.links if not l.end_link]):
            return False

        # Ensure model is accessible
        self.model = model
        
        # # A decay vertex has one non-zero mass parameter coming
        # # only once
        # decay_vertex = False
        # masses = [p.get('mass') for p in \
        #           model.get_interaction(self.vertex_id).get('particles') \
        #           if p.get('mass').lower() != 'zero']

        # for mass in set(masses):
        #     if len([m for m in masses if m == mass]) == 1:
        #         decay_vertex = True
        #         break

        # # Decay vertices can't fail
        # if decay_vertex:
        #     return True
        if self.depth == 1:
            # Only (1,2) is forbidden among 2-particle combinations
            return sorted(self.get_external_numbers()) != [1,2]

        # Note that self.depth > 1 from now on

        # We never allow long s-channel chains
        external_numbers = [sorted(l.get_external_numbers()) for l in \
                            self.links]
        sum_external = sum(external_numbers, [])

        # Fail if this is pure s-channel
        if 1 not in sum_external and 2 not in sum_external:
            return False

        # If this is amplitude, fail if any single link contains both
        # 1 and 2 (meaning that this is an s-channel vertex)
        if amp_link and any([ext[:2] == [1,2] for \
                             ext in external_numbers]):
            return False

        # Require an external parton if the second link is < tch_depth - 1
        depths = sorted([(l.depth,external_numbers[i]) for (i, l) \
                         in enumerate(self.links)])

        if depths[0][0] == 0 and depths[1][0] == 0:
            # This is final vertex in 2->2 process
            return [depths[0][1][0],depths[1][1][0]] != [1,2]

        # This checks tch_depth for s-channel mergings for wf-like vertices
        #if not amp_link and depths[1][0] < tch_depth and \
        #       (depths[0][0] == 0 and depths[0][1][0] <= 2 or \
        #        depths[0][0] > 0 and depths[1][0] < tch_depth-1):
        #    return False

        # This checks that final vertex is t-channel
        #if amp_link and depths[1][0] < tch_depth:
        #    if depths[0][0] == 0 and depths[0][1][0] <= 2 or \
        #       depths[0][0] > 0 and depths[1][0] < tch_depth-1:
        #        print 'fail 282'
        #        return True

        # Otherwise, return True
        return True

    def fill_comp_array(self, amp_link = False, identify_depth = 1):
        """Fills a comparison array which allows fast comparison
        between tags, containing [depth] (if not external) and
        [depth, external number] (if external)."""

        if self.end_link:
            return [[0] + self.get_external_numbers() + \
                        [self.links[0][0][1]]]

        res_array = []
        links = self.links
        left_links = []
        depth = self.depth
        if amp_link:
            # This is the amplitude (i.e., the original call)
            depths = sorted([(d,i) for (i,d) in \
                             enumerate([l.depth for l in self.links])])
            depth = depths[0][0]+depths[1][0]+1
            # links are the first two links (out of 3)
            links = [self.links[d[1]] for d in depths[:2]]
            # left_links is the last link (with largest depth)
            left_links.append(self.links[depths[2][1]])
        for link in links:
            # Recursive call
            res_array.extend(link.fill_comp_array(False, identify_depth))
        res_array.sort()
        if depth <= identify_depth:
            #print "depth: ",depth,"res_array: ",res_array
            res_array.insert(0, [depth])
            # Remove particle id, since not needed
            for i,entry in enumerate(res_array):
                if len(entry) == 3:
                    res_array[i].pop(2)
            res_array = [self.flatten(res_array)]
            #print "flattened res_array: ",res_array
        else:
            # Remove number and keep only id
            res_array.sort()
            res_array.insert(0, [depth])
            for i,entry in enumerate(res_array):
                if len(entry) == 3:
                    res_array[i].pop(1)
            res_array = [self.flatten(res_array)]
            #print "large depth: ",depth,"res_array: ",res_array
        for link in left_links:
            res_array.extend(link.fill_comp_array(False, identify_depth))
        #print "final res_array before sort: ",res_array
        res_array.sort()
        return res_array

    @staticmethod
    def flatten(array):
        """Flatten a list [N,[list],[list],...] to [N,list,list,...]"""
        depth = []
        if isinstance(array[0], int): depth = [array.pop(0)]
        return depth + sum(array,[])
    
class PeriferalDiagramTag(OrderDiagramTag):
    """DiagramTag daughter class to select only periferal diagrams, as
    needed for MadEvent phase space integration of color ordered
    amplitudes"""

    link_class = PeriferalDiagramTagChainLink

    def check_periferal_diagram(self, model, order = 1,
                                nallowed = -1,
                                include_all_t = False):
        """Check if this is a periferal diagram (to the order given).

        Periferal means that at most one (massless) external leg
        attaches to a vertex with depth > order. Vertices where all
        particles are different (decays) are not included in this
        condition. Returns the number of external legs attached to a
        vertex with depth > order, and a sorted list of tuples of the
        external numbers of the vertices, e.g. 1, [(1,3),(2,4),5] for a
        t-channel 5-particle diagram, or an empty list if the diagram
        is not periferal."""

        # Algorithm: Go through the links and collect external
        # particles. Allow 1 external particle coupling to a vertex if
        # order == 1, otherwise if any external couples to depth >
        # order, return []

        if nallowed < 0:
            nallowed = 1 if order == 1 else 0
                                             
        return self.tag.check_periferal_legs(model,
                                             order,
                                             nallowed,
                                             include_all_t,
                                             amp_link = True)[1]

    def pass_restrictions(self, model, tch_depth = 1):
        """Check if this periferal diagram passes the following
        restrictions:

        - No s-channel propagators with more than 2 external particles
        - If include_only_t, only pure t-channels allowed (no s-channels)
        - If not allow_group_12, no (1,2) combinations allowed
        - If not allow_12_sch, initial state legs must couple directly
          to external particles (not to s-channel propagators).
        """
        #return True
        return self.tag.pass_restrictions(model, True, tch_depth)

    def get_comp_array(self, identify_depth = 1):
        """Give a comparison array which allows fast comparison
        between tags"""

        comp_list = self.tag.fill_comp_array(amp_link = True,
                                             identify_depth = identify_depth)
        #print "comp_array: ",comp_list
        return array.array('i', sum(comp_list, []))

#===============================================================================
# ColorOrderedLeg
#===============================================================================
class ColorOrderedLeg(base_objects.Leg):
    """Leg with two additional flags: fermion line number, and color
    ordering number. These will disallow all non-correct orderings of
    colored fermions and gluons. So far, only color triplets and color
    octets are treated (i.e., no higher color representations or
    epsilon^{ijk})"""

    def default_setup(self):
        """Default values for all properties"""

        super(ColorOrderedLeg, self).default_setup()
        # Flag to keep track of color ordering
        # For colored particles: {Cycle: (upper, lower)}
        # For color singlets: {}
        self['color_ordering'] = {}
    
    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name == 'color_ordering':
            if not isinstance(value, dict):
                raise self.PhysicsObjectError, \
                        "%s is not a valid dict object" % str(value)
        else:
            super(ColorOrderedLeg, self).filter(name, value)

        return True

    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""

        return ['id', 'number', 'state', 'from_group',
                'color_ordering']

#===============================================================================
# ColorOrderedAmplitude
#===============================================================================
class ColorOrderedAmplitude(diagram_generation.Amplitude):
    """ColorOrderedAmplitude: Set of color flows. The setup_process
    subroutine generates ColorOrderedFlow objects corresponding to all
    unique color flows for the process.

    Note that if we have more than two color triplet lines in the
    process, then different color flows are generated corresponding to
    different numbers of singlet gluons.
    """

    def __init__(self, argument=None):
        """Allow color-ordered initialization with Process"""

        if isinstance(argument, base_objects.Process):
            super(ColorOrderedAmplitude, self).__init__()
            self.set('process', argument)
            self.generate_diagrams()
        elif argument != None:
            # call the mother routine
            super(ColorOrderedAmplitude, self).__init__(argument)
        else:
            # call the mother routine
            super(ColorOrderedAmplitude, self).__init__()

    def default_setup(self):
        """Add number of gluon flag to Amplitude members"""
        
        super(ColorOrderedAmplitude, self).default_setup()
        self['color_flows'] = ColorOrderedFlowList()

    def get(self, name):
        """Special get for diagrams"""

        if name == 'diagrams':
            return base_objects.DiagramList(sum([cf.get('diagrams') \
                                for cf in self.get('color_flows')],[]))
        return super(ColorOrderedAmplitude, self).get(name)

    def generate_diagrams(self):
        """Add singlets for each color octet in model. Setup
        ColorOrderedFlows corresponding to all color flows of this
        process. Each ordering of unique particles with color triplets
        pairwise combined as (3bar 3)...(3bar 3)... etc. corresponds
        to a unique color flow. If multiple triplet lines are present,
        one flow is generated for each number of singlet exchanges
        (specified by the coupling order singlet_QCD)"""

        process = self.get('process')

        # Check that this is a valid process
        if not process.check_valid_process():
            return base_objects.DiagramList()

        # Ensure that all legs are unique
        process.set('legs',
                    base_objects.LegList([copy.copy(l) for l in \
                                          process.get('legs')]))
        
        # Set leg numbers for the process
        for i, leg in enumerate(process.get('legs')):
            leg.set('number', i+1)

        legs = base_objects.LegList([copy.copy(l) for l in \
                                     process.get('legs')])

        if not isinstance(process.get('model'), ColorOrderedModel):
            process.set('model', ColorOrderedModel(process.get('model')))

        model = process.get('model')
        
        # Add color negative singlets to model corresponding to all
        # color octets, with only 3-3bar-8 interactions.
        # Color factor: -1/(2*Nc) * Id(1,2)

        # Reverse initial state legs to have all legs outgoing
        for i, leg in enumerate(legs):
            if not leg.get('state'):
                leg.set('id', model.get_particle(leg.get('id')).\
                        get_anti_pdg_code())

        # Create list of leg (number, id, color) to find unique color flows
        order_legs = [(l.get('number'), l.get('id'),
                       model.get_particle(l.get('id')).get_color()) \
                      for l in legs]

        # Identify particle types: Color singlets (can connect
        # anywhere), color triplets, anti-triplets and octets
        triplet_legs = [l for l in order_legs if l[2] == 3]
        anti_triplet_legs = [l for l in order_legs if l[2] == -3]
        octet_legs = [l for l in order_legs if abs(l[2]) == 8]
        singlet_legs = [l for l in order_legs if abs(l[2]) == 1]

        if len(triplet_legs) != len(anti_triplet_legs):
            # Need triplets in pairs for valid amplitude
            # (color triplet-epsilon^{ijk} is not allowed)
            raise MadGraph5Error, \
                  "# triplets != # anti-triplets in color ordered amplitude"
            

        if len(triplet_legs + anti_triplet_legs + octet_legs + singlet_legs) != \
               len(legs):
            raise MadGraph5Error, \
                  "Non-triplet/octet/singlet found in color ordered amplitude"
        
        # Determine all unique permutations of particles corresponding
        # to color flows. Valid color flows have color triplet pairs
        # organized as (3bar 3)...(3bar 3)...

        # Extract all unique permutations of legs
        leg_perms = []
        good_perm_ids = []
        same_flavor_perms = {}
        if anti_triplet_legs:
            last_leg = [anti_triplet_legs.pop(0)]
        else:
            last_leg = [octet_legs.pop(-1)]

        for perm in itertools.permutations(sorted(anti_triplet_legs + \
                                                  triplet_legs + \
                                                  octet_legs,
                                                  lambda l1, l2: l1[0] - l2[0])):
            # Remove permutations where the order of anti-triplets is changed
            anti_triplet_numbers = [p[0] for p in perm if p[2]==-3]
            if anti_triplet_numbers != sorted(anti_triplet_numbers):
                continue

            # Permutations with same flavor ordering as existing perm
            # are kept in same_flavor_perms
            perm_ids = array.array('i', [p[1] for p in perm])
            try:
                ind = good_perm_ids.index(perm_ids)
                same_flavor_perms[ind].append([p[0] for p in perm] + \
                                              [last_leg[0][0]] + \
                                              [l[0] for l in singlet_legs])
                continue
            except ValueError:
                pass
            
            # Remove permutations where a triplet and
            # anti-triplet are not next to each other and ordered as 
            # (3bar 3)...(3bar 3)...
            failed = False
            # trip is a list of [position, (number, id, color)] for
            # triplet and antitriplet legs
            trip = [(i,p) for (i,p) in enumerate(last_leg + list(perm)) \
                    if abs(p[2])==3]
            for i in range(0,len(trip),2):
                if trip[i][0] != trip[i+1][0]-1 or \
                       trip[i][1][2]-trip[i+1][1][2] != -6:
                    failed = True
                    break
            if failed:
                continue

            # Initiate list of permutations for this leglist
            same_flavor_perms[len(leg_perms)] = [[p[0] for p in perm] + \
                                                 [last_leg[0][0]] + \
                                                 [l[0] for l in singlet_legs]]
            # Add permutation ids for comparison above
            good_perm_ids.append(perm_ids)
            # Add permutation to accepted permutations
            leg_perms.append(perm)

        # Create color flow amplitudes corresponding to all resulting
        # permutations
        color_flows = ColorOrderedFlowList()
        colegs = base_objects.LegList([ColorOrderedLeg(l) for l in legs])

        # Set color flow flags (color_ordering) so that every
        # chain (3 8 8 .. 8 3bar) has a unique color ordering
        # group, and each entry in the chain has a unique color flow
        # flag {chain:(n,n)}

        ntriplets = len(triplet_legs)
        # Create the color ordered flows for all combinations
        # numbers of octet and singlet gluon
        for iflow in range(max(1, ntriplets)):
            for iperm, perm in enumerate(leg_perms):
                colegs = base_objects.LegList([ColorOrderedLeg(l) for l in legs])
                # Keep track of number of triplets
                ichain = 0
                ileg = 0
                # Set color ordering flags for all colored legs
                for perm_leg in list(perm) + last_leg:
                    leg = colegs[perm_leg[0]-1]
                    if perm_leg[2] == 3:
                        ichain += 1
                        ileg = 0
                    ileg += 1
                    leg.set('color_ordering', {ichain: (ileg, ileg)})

                # Restore initial state leg identities
                for leg in colegs:
                    if not leg.get('state'):
                        leg.set('id', model.get_particle(leg.get('id')).\
                                get_anti_pdg_code())

                # Create the color ordered process
                coprocess = copy.copy(process)
                coprocess.set('orders', copy.copy(process.get('orders')))
                coprocess.set('legs', colegs)
                if 'QCD' in process.get('orders'):
                    coprocess.get('orders')['QCD'] = \
                                        process.get('orders')['QCD'] - 2*iflow
                    if coprocess.get('orders')['QCD'] < 0:
                        continue
                coprocess.get('orders')['singlet_QCD'] = 2*iflow

                # Create and generate diagrams for the color flow
                flow = ColorOrderedFlow(coprocess)
                
                if flow.get('diagrams'):
                    # Set permutation information for this flow
                    # (relative to first permutation)
                    perms = same_flavor_perms[iperm]
                    flow.set('permutations',
                             [base_objects.reorder_permutation(perms[0], p) \
                              for p in perms])
                    # Add flow to list
                    color_flows.append(flow)
                    if not 'QCD' in process.get('orders'):
                        # Set coupling orders for the process
                        process.set('orders', coprocess.get('orders'))
                        process.get('orders')['QCD'] = \
                                          process.get('orders')['QCD'] + 2*iflow
                        process.get('orders')['singlet_QCD'] = 0
                    # Sort diagrams for this color flow using OrderDiagramTags
                    for idiag, diag in enumerate(flow.get('diagrams')):
                        flow.get('diagrams')[idiag] = \
                                      OrderDiagramTag(diag).\
                                                         diagram_from_tag(model)
            # If there are no diagrams for singlet_QCD=0, there are no
            # diagrams at all for this process. Cancel generation.
            if iflow == 0 and not color_flows:
                break

        self.set('color_flows', color_flows)
        return self.get('diagrams')

    def get_periferal_diagrams_from_flows(self, include_all_t = False,
                                          tch_depth = 1,
                                          identify_depth = 10):
        """Generate the periferal diagrams needed for efficient phase
        space integration from the diagrams of the color flows and
        their permutations, using the PeriferalDiagramTag.

        Use all depth-2 periferal diagrams from the basic flows, and
        then generate all diagrams from all permutations before
        applying the following rules:

        Rules:
        - External lines can couple only to vertices of max depth 2
        - No (1,2) (pure s-channel) combinations
        - No connection of external line to s-channel propagator

        Remember which flows contribute to each diagram, since this
        can be used to select flows for a given channel in the phase
        space integration."""

        process = self.get('color_flows')[0].get('process')
        nfinal = len([l for l in process.get('legs') if l.get('state')])
        model = process.get('model')
        all_gluons = not any([l.get('id') != 21 for l in process.get('legs')])

        basic_diagrams = []
        basic_tags = []
        all_diagrams = []
        all_tags = []
        failed_tags = []

        # Pick out the periferal diagrams from the basic flows
        idiag=0
        done = False
        while not done:
            for iflow, flow in enumerate(self.get('color_flows')):
                if flow.get('process').get('orders')['singlet_QCD'] > 0:
                    continue
                for diag in flow.get('diagrams'):
                    tag = PeriferalDiagramTag(diag)
                    if not tag.check_periferal_diagram(model, order = 2,
                                                       include_all_t = \
                                                           include_all_t):
                        # This diagram is not periferal
                        continue

                    idiag += 1
                    tag_array = tag.get_comp_array(identify_depth = \
                                                       identify_depth)
                    # If the diagram is not already represented, add it to
                    # basic_diagrams. If this is all-gluon amplitude, keep
                    # only t-channel diagrams.
                    if tag_array not in basic_tags and \
                       (not all_gluons or \
                        tag.pass_restrictions(model, tch_depth = tch_depth)):
                        # Append tag to basic_tag
                        basic_tags.append(tag_array)
                        # Append diagram to basic_diagrams
                        basic_diagrams.append(diag)
            
            # Check that there are some t-channel diagrams in the base flows
            if not any([PeriferalDiagramTag(diag).pass_restrictions(model, 
                                                    tch_depth = tch_depth) \
                            for diag in basic_diagrams]) and tch_depth > 0:
                tch_depth -= 1
            else:
                done = True

        assert(basic_diagrams)
        
        # Now go through all permutations to get the full set of diagrams
        first_perm = None
        iperm = 0
        comp_list = self.get_comp_list(process)
        for perm in itertools.permutations(range(1,len(process.get('legs'))+1)):
            if any([comp_list[perm[i]-1] != comp_list[i] for i in \
                    range(len(process.get('legs')))]): continue
            if not first_perm: first_perm = perm
            iperm = iperm + 1
            # Go through the basic diagrams and replace indices.
            perm_dict = dict(zip(first_perm, perm))
            for basic_diag in basic_diagrams:
                diag = basic_diag.renumber_legs(perm_dict,
                                              process.get('legs'))
                tag = PeriferalDiagramTag(diag)
                # Use get_comp_array for very fast comparison between tags
                tag_array = tag.get_comp_array(identify_depth = identify_depth)
                # Check if diagram already failed
                if tag_array in failed_tags:
                    continue

                # Check if the diagram is already represented,
                # otherwise append tag, diagram, and (flow, permutation)
                try:
                    index = all_tags.index(tag_array)
                except ValueError:
                    # Check if this diagram passes the rules to be used
                    # for phase space integration
                    if tag.pass_restrictions(model, tch_depth = tch_depth):
                        all_tags.append(tag_array)
                        all_diagrams.append(diag)
                    else:
                        failed_tags.append(tag_array)

        return base_objects.DiagramList(all_diagrams), tch_depth

    @staticmethod
    def get_comp_list(process):
        """Return a list of numbers that are identical for particles
        that are permutated"""

        model = process.get('model')
        part_states = [(model.get_particle(l.get('id')), l.get('state')) \
                       for l in process.get('legs')]
        pdg_colors = [(p[0].get_pdg_code(), p[0].get_color()) if p[1] else \
                     (p[0].get_anti_pdg_code(), p[0].get_anti_color()) \
                     for p in part_states]
        pdg_dict = {}
        comp_list = []
        comp_id = 0
        # For all-octet amplitudes, exclude the first particle
        octets_only = pdg_colors[0][1] == 8 and \
                      len(set([p[1] for p in pdg_colors])) == 1
        for ipart, (pdg, col) in enumerate(pdg_colors):
            if pdg in pdg_dict and col in [3, 8] and \
                   (not octets_only or ipart > 1):
                comp_list.append(pdg_dict[pdg])
            else:
                comp_id += 1
                comp_list.append(comp_id)
                pdg_dict[pdg] = comp_id

        return comp_list
        

#===============================================================================
# ColorOrderedMultiProcess
#===============================================================================
class ColorOrderedMultiProcess(diagram_generation.MultiProcess):
    """Version of MultiProcess which generates ColorOrderedAmplitudes"""

    amplitude_class = ColorOrderedAmplitude

    @staticmethod
    def cross_amplitude(amplitude, org_perm, new_perm):
        """Return the color ordered amplitude crossed with the
        permutation new_perm"""
        
        # Initiate new amplitude
        new_amp = copy.copy(amplitude)
        # Create dict from original leg numbers to new leg numbers
        perm_map = dict(zip(org_perm, new_perm))
        # New process with reordered legs
        process = amplitude.get('process').renumber_legs(perm_map)
        # Set process
        new_amp.set('process', process)
        # Cross all color flows
        color_flows = ColorOrderedFlowList()
        for flow in amplitude.get('color_flows'):
            color_flows.append(diagram_generation.MultiProcess.\
                               cross_amplitude(flow, org_perm, new_perm))
            color_flows[-1].reorder_permutations(perm_map)
        new_amp.set('color_flows', color_flows)
        # Make sure to reset mirror process
        new_amp.set('has_mirror_process', False)
        return new_amp

#===============================================================================
# ColorOrderedFlow
#===============================================================================
class ColorOrderedFlow(diagram_generation.Amplitude):
    """Initialize with a color ordered process (created by
    ColorOrderedAmplitude), then call generate_diagrams() to generate
    the set of color ordered diagrams for the color flow.
    """

    def __init__(self, argument=None):
        """Allow color-ordered initialization with Process"""

        if isinstance(argument, base_objects.Process):
            super(ColorOrderedFlow, self).__init__()
            self.set('process', argument)
            # Set max color order for all color ordering groups
            legs = argument.get('legs')
            groups = set(sum([l.get('color_ordering').keys() for l in legs],[]))
            self.set('max_color_orders', dict( \
                [(g, max([l.get('color_ordering')[g][0] for l \
                      in legs if g in l.get('color_ordering')]))\
                 for g in groups]))
            self.generate_diagrams()
        elif argument != None:
            # call the mother routine
            super(ColorOrderedFlow, self).__init__(argument)
        else:
            # call the mother routine
            super(ColorOrderedFlow, self).__init__()

    def default_setup(self):
        """Add number of color orderings and permutations to Amplitude members"""
        
        super(ColorOrderedFlow, self).default_setup()
        self['max_color_orders'] = {}
        self['permutations'] = []

    def get_combined_legs(self, legs, leg_vert_ids, number, state):
        """Determine if the combination of legs is valid with color
        ordering. Algorithm: A combination is valid if

        1) All particles with the same color ordering group have
        adjacent color flow numbers

        2) No color octet is dead-ending (if there are other color groups)

        3) Resulting legs have: a) no uncompleted groups (if color singlet),
           b) exactly one uncompleted group (if color triplet),
           c) no uncompleted group, and 1 or 2 groups (if color octet).
        """

        # Find all color ordering groups
        groups = set(sum([l.get('color_ordering').keys() for l in legs],[]))
        model = self.get('process').get('model')
        new_leg_colors = dict([(leg_id,
                                model.get_particle(leg_id).get('color')) \
                           for (leg_id, vert_id) in leg_vert_ids])
        old_leg_colors = dict([(leg_id,
                                model.get_particle(leg_id).get('color')) \
                               for leg_id in set([l.get('id') for l in legs])])
        color_orderings = dict([(leg_id, {}) for (leg_id, vert_id) in \
                                leg_vert_ids])

        for group in groups:
            # sort the legs with color ordering number
            color_legs = [(i, l.get('color_ordering')[group]) for \
                          (i, l) in enumerate(legs) \
                          if group in l.get('color_ordering')]

            # Check for unattached color octet legs
            if len(groups) >= 2 and len(color_legs) == 1 and \
               old_leg_colors[legs[color_legs[0][0]].get('id')] == 8 and \
                   len(legs[color_legs[0][0]].get('color_ordering')) == 1:
                return []            

            color_legs.sort(lambda l1, l2: l1[1][0] - l2[1][0])

            color_ordering = (color_legs[0][1][0], color_legs[-1][1][1])

            # Check that we don't try to combine legs with
            # non-adjacent color ordering (allowing to wrap around
            # cyclically)
            lastmax = color_legs[0][1][1]
            ngap = 0
            # First check if there is a gap between last and first
            if (color_legs[0][1][0] > 1 or \
                color_legs[-1][1][1] < self.get('max_color_orders')[group]) and \
                color_legs[-1][1][1] + 1 != color_legs[0][1][0]:
                ngap = 1
            # Check if there are gaps between other legs
            for leg in color_legs[1:]:
                if leg[1][0] != lastmax + 1:
                    ngap += 1
                    color_ordering = (leg[1][0], lastmax)
                if ngap == 2:
                    # For adjacent legs, only allow one gap
                    return []
                lastmax = leg[1][1]
            # Set color ordering for new legs
            for leg_id in color_orderings.keys():
                color_orderings[leg_id][group] = color_ordering
                
        # Check validity of resulting color orderings for the legs
        ileg = 0
        while ileg < len(leg_vert_ids):
            leg_id, vert_id = leg_vert_ids[ileg]
            if abs(new_leg_colors[leg_id]) == 1:
                # Color singlets need all groups to be completed
                if any([color_orderings[leg_id][group] != \
                        (1, self.get('max_color_orders')[group]) \
                        for group in groups]):
                    leg_vert_ids.remove((leg_id, vert_id))
                else:
                    ileg += 1
                color_orderings[leg_id] = {}
            elif abs(new_leg_colors[leg_id]) == 3:
                # Color triplets should have exactly one color ordering
                color_orderings[leg_id] = \
                    dict([(group, color_orderings[leg_id][group]) for \
                          group in groups if \
                          color_orderings[leg_id][group] != \
                           (1, self.get('max_color_orders')[group])])
                if len(color_orderings[leg_id].keys()) != 1:
                    leg_vert_ids.remove((leg_id, vert_id))                    
                else:
                    ileg += 1
            elif abs(new_leg_colors[leg_id]) == 8:
                # Color octets should have 1 or 2 valid orderings
                # valid_orderings is non-completed groups
                valid_orderings = [group for group in groups if \
                                   color_orderings[leg_id][group] != \
                                   (1, self.get('max_color_orders')[group])]
                if (len(valid_orderings) < 1 or \
                    len(valid_orderings) > 2):
                    leg_vert_ids.remove((leg_id, vert_id))
                else:
                    ileg += 1
                # Remove completed groups
                for group in groups:
                    if group not in valid_orderings:
                        del color_orderings[leg_id][group]
                
        # Return all legs that have valid color_orderings
        mylegs = [(ColorOrderedLeg({'id':leg_id,
                                   'number':number,
                                   'state':state,
                                   'from_group':True,
                                   'color_ordering': color_orderings[leg_id]}),
                   vert_id) for (leg_id, vert_id) in leg_vert_ids]

        return mylegs

    def get_combined_vertices(self, legs, vert_ids):
        """Check that all color-ordering groups are pairwise
        connected, i.e., we do not allow groups that are disjuct."""

        # Find all color ordering groups
        groups = set(sum([l.get('color_ordering').keys() for l in legs],[]))
        # Extract legs colored under each group
        group_legs = {}
        for group in groups:
            group_legs[group] = set([l.get('number') for l in legs \
                                     if group in l.get('color_ordering')])
        # Check that all groups are pair-wise connected by particles
        connected_groups = []
        if len(groups) > 1:
            for g1, g2 in itertools.combinations(groups, 2):
                if group_legs[g1].intersection(group_legs[g2]):
                    connected_groups.append((g1,g2))
            if len(groups) == 2 and len(connected_groups) != 1 or \
               len(groups) > 2 and len(connected_groups) != len(groups):
                return []

        return vert_ids

    def reorder_permutations(self, perm_map):
        """Reorder permutations according to perm_map."""

        # Replace numbers in the first permutation
        first_perm = [(v,i) for (i,v) in enumerate([perm_map[p] for p in \
                      self.get('permutations')[0]])]
        # Get the ordering to sort the first permutations
        first_perm_order = [i for (v, i) in sorted(first_perm)]
        new_perms = []
        for perm in self.get('permutations'):
            # First replace numbers
            replaced_perm = [perm_map[p] for p in perm]
            # Then reorder permutation according to first permutation order
            new_perms.append([replaced_perm[first_perm_order[i]] for i in \
                              range(len(perm))])
        # Replace old permutations with new ones
        self.set('permutations', new_perms)

    
#===============================================================================
# ColorOrderedFlowList
#===============================================================================
class ColorOrderedFlowList(diagram_generation.AmplitudeList):
    """List of ColorOrderedFlow objects
    """

    def is_valid_element(self, obj):
        """Test if object obj is a valid Amplitude for the list."""

        return isinstance(obj, ColorOrderedFlow)


#===============================================================================
# ColorOrderedModel
#===============================================================================
class ColorOrderedModel(import_ufo.RestrictModel):
    """When initiated with a Model, adds negative singlets for all octets."""
    
    # Customized constructor
    def __init__(self, *arguments):
        """Generate a ColorOrderedModel from a Model
        """

        super(ColorOrderedModel, self).__init__(*arguments)
        
        if len(arguments) == 1 and isinstance(arguments[0], base_objects.Model) \
               and not isinstance(arguments[0], ColorOrderedModel):
            self.add_color_singlets()

    def add_color_singlets(self):
        """Go through model and add color singlets for all color
        octets, with couplings only to triplet pairs. The color factor
        will be multiplied by -Nc/2 for every two singlet_QCD orders
        when an amplitude containing color singlets is created."""

        particles = copy.copy(self.get('particles'))
        interactions = copy.copy(self.get('interactions'))
        # Pick out all color octets
        octets = [p for p in particles if p.get('color') == 8]
        # Create corresponding color singlets with new pdg codes
        singlets = [copy.copy(o) for o in octets]
        pdgs = [p.get('pdg_code') for p in particles]
        free_pdgs = [i for i in \
                     range(1,len(self.get('particles'))+len(singlets)+1) if \
                     i not in pdgs]
        octet_singlet_dict = {}
        for i, singlet in enumerate(singlets):
            singlet.set('color', 1)
            singlet.set('pdg_code', free_pdgs[i])
            octet_singlet_dict[octets[i].get('pdg_code')] = \
                                                      singlet.get('pdg_code')
            if not singlet.get('self_antipart'):
                octet_singlet_dict[octets[i].get_anti_pdg_code()] = \
                                                 singlet.get_anti_pdg_code()

        particles.extend(singlets)
        self.set('particles', particles)

        # Now pick out all 8-3-3bar interactions
        octet_interactions = [inter for inter in interactions if \
                              sorted([p.get_color() for p in \
                                      inter.get('particles') if \
                                      p.get_color() != 1]) == [-3, 3, 8]]

        # Create corresponding singlet interactions
        singlet_interactions = [copy.copy(i) for i in octet_interactions]
        max_inter_id = max([i.get('id') for i in interactions])
        for inter in singlet_interactions:
            # Set unique interaction id
            max_inter_id += 1
            inter.set('id', max_inter_id)
            # Get particles and indices for triplet/antitriplet and octet
            parts = copy.copy(inter.get('particles'))
            colors = [p.get_color() for p in parts]
            trip_ind = colors.index(3)
            antitrip_ind = colors.index(-3)
            oct_ind = colors.index(8)
            # Replace octet with singlet in interaction
            parts[oct_ind] = self.get_particle(\
                             octet_singlet_dict[parts[oct_ind].get_pdg_code()])
            inter.set('particles', parts)
            # Set color string to T(3, 3bar)
            col_str = color.ColorString([color.T(trip_ind,
                                                                 antitrip_ind)])
            col_str.Nc_power = -1
            inter.set('color', [col_str])
            # Set orders, replacing QCD with singlet_QCD
            orders = copy.copy(inter.get('orders'))
            assert 'QCD' in orders
            del orders['QCD']
            orders['singlet_QCD'] = 1
            inter.set('orders', orders)

        interactions.extend(singlet_interactions)
        order_hierarchy = self.get('order_hierarchy')
        order_hierarchy['singlet_QCD'] = order_hierarchy['QCD']
        expansion_order = self.get('expansion_order')
        self.set('interactions', interactions)
        self.set('order_hierarchy', order_hierarchy)
        self.set('expansion_order', expansion_order)
        return 
