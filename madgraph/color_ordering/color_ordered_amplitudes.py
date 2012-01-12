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
                      model.get_interaction(self.vertex_id).get('particles') \
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

    def pass_restrictions(self, model, amp_link = False,
                          include_only_t = False,
                          allow_group_12 = False,
                          allow_12_sch = True):
        """Check if this subdiagram passes the following
        restrictions:

        - Maximum two external legs connecting to any s-channel
          propagator
        - If include_only_t, only allow pure t-channel diagrams
        - 
        """

        # First check if any daughter links fail
        if any([not l.pass_restrictions(model, amp_link,
                                        include_only_t,
                                        allow_group_12,
                                        allow_12_sch) \
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


        # If we only allow pure t-channels, fail if a level 1 vertex has non-t
        if include_only_t:
            if self.depth == 1:
                external_numbers = self.get_external_numbers()
                if 1 not in external_numbers and 2 not in external_numbers \
                   or set(external_numbers) == set([1,2]):
                    return False
            # Otherwise pass
            return True

        # If we allow non-pure t-channels, continue
        
        if allow_group_12:
            # If allow_group_12, allow grouping 1 and 2
            if self.depth == 1:
                return True
            elif not allow_12_sch:
                # Forbid any diagram which combines 1 or 2 with s-channel
                externals = [l.get_external_numbers()[0] for l in self.links \
                             if l.end_link] 
                if len(externals) == 1 and externals[0] in [1,2]:
                    # We are pairing 1 or 2 with an s-channel
                    return False
        else:
            end_links = [l.get_external_numbers()[0] for l in self.links \
                         if l.end_link] 

            if self.depth == 1:
                # If this vertex has depth 1, forbid the (1,2) combination
                if set(end_links) == set([1,2]):
                    return False
                # Otherwise pass
                return True

            # Forbid any diagram which combines 1 or 2 with s-channel
            if not allow_12_sch and \
                   len(end_links) == 1 and end_links[0] in [1,2]:
                return False
            elif set(end_links) == set([1,2]):
                # Forbid the combination (1,2)
                return False

        # The following applies whether or not allow_group_12
        # Note that self.depth > 1 below

        # Fail if this is pure s-channel
        external_numbers = sum([l.get_external_numbers() for l in \
                                self.links], [])

        if 1 not in external_numbers and 2 not in external_numbers:
            return False

        # If this is amplitude, fail if any single link contains both
        # 1 and 2 (meaning that this is an s-channel) unless allow_group_12
        if amp_link and not allow_group_12:
            return not any([sorted(l.get_external_numbers())[:2] == [1,2] for \
                            l in self.links])

        # Otherwise, return True
        return True

    def fill_comp_array(self, comp_array):
        """Fills a comparison array which allows fast comparison
        between tags, containing [depth] (if not external) and
        [depth, external number] (if external)."""

        comp_array.append(self.depth)
        if self.end_link:
            comp_array.append(self.links[0][1])
        else:
            for link in self.links:
                link.fill_comp_array(comp_array)

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

    def pass_restrictions(self, model,
                          include_only_t = False,
                          allow_group_12 = False,
                          allow_12_sch = True):
        """Check if this periferal diagram passes the following
        restrictions:

        - No s-channel propagators with more than 2 external particles
        - If include_only_t, only pure t-channels allowed (no s-channels)
        - If not allow_group_12, no (1,2) combinations allowed
        - If not allow_12_sch, initial state legs must couple directly
          to external particles (not to s-channel propagators).
        """
        
        return self.tag.pass_restrictions(model, True,
                                          include_only_t,
                                          allow_group_12,
                                          allow_12_sch)

    def get_comp_array(self):
        """Give a comparison array which allows fast comparison
        between tags"""

        comp_array = array.array('i')
        self.tag.fill_comp_array(comp_array)
        return comp_array

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

            # Create the color ordered flows for all combinations
            # numbers of octet and singlet gluon
            for iflow in range(max(1, ichain)):
                coprocess = copy.copy(process)
                coprocess.set('orders', copy.copy(process.get('orders')))
                coprocess.set('legs', colegs)
                if 'QCD' in process.get('orders'):
                    coprocess.get('orders')['QCD'] = \
                                        process.get('orders')['QCD'] - 2*iflow
                    if coprocess.get('orders')['QCD'] < 0:
                        continue
                coprocess.get('orders')['singlet_QCD'] = 2*iflow
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
                    
        self.set('color_flows', color_flows)
        return self.get('diagrams')

    def get_periferal_diagrams_from_flows(self,
                                          include_only_t = False,
                                          allow_group_12 = False,
                                          allow_12_sch = True):
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
        model = process.get('model')

        basic_diagrams = []
        basic_tags = []
        all_diagrams = []
        all_tags = []
        flow_permutations = []
        failed_tags = []

        # Pick out the periferal diagrams from the basic flows
        for iflow, flow in enumerate(self.get('color_flows')):
            for diag in flow.get('diagrams'):
                tag = PeriferalDiagramTag(diag)
                if not tag.check_periferal_diagram(model, order = 2):
                    # This diagram is not periferal
                    continue
                tag_array = tag.get_comp_array()
                # If the diagram is not already represented, add it to
                # basic_diagrams
                if tag_array not in basic_tags:
                    # Append tag to basic_tag
                    basic_tags.append(tag_array)
                    # Append diagram to basic_diagrams
                    basic_diagrams.append(diag)

        # Now go through all permutations to get the full set of diagrams
        permutations = self.get('color_flows')[0].get('permutations')
        for iperm, perm in enumerate(permutations):
            # Go through the basic diagrams and replace indices.
            perm_dict = dict(zip(permutations[0],perm))
            for basic_diag in basic_diagrams:
                diag = basic_diag.renumber_legs(perm_dict,
                                              process.get('legs'))
                tag = PeriferalDiagramTag(diag)
                # Use get_comp_array for very fast comparison between tags
                tag_array = tag.get_comp_array()
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
                    if tag.pass_restrictions(model,
                                             include_only_t,
                                             allow_group_12,
                                             allow_12_sch):
                        all_tags.append(tag_array)
                        all_diagrams.append(diag)
                        flow_permutations.append([(iflow,iperm)])
                    else:
                        failed_tags.append(tag_array)
                else:
                    flow_permutations[index].append((iflow, iperm))

        return base_objects.DiagramList(all_diagrams), flow_permutations

    @staticmethod
    def get_comp_list(process):
        """Return a list of numbers that are identical for particles
        that are permutated"""

        model = process.get('model')
        part_states = [(model.get_particle(l.get('id')), l.get('state')) \
                       for l in process.get('legs')]
        pdg_colors = [(p[0].get_pdg_code(), p[0].get_color()) if p[1] else \
                     (p[0].get_anti_pdg_code(), p[0].get_color()) \
                     for p in part_states]
        pdg_dict = {}
        comp_list = []
        comp_id = 0
        # For all-octet amplitudes, exclude the first particle
        octets_only = pdg_colors[0][1] == 8 and \
                      len(set([p[1] for p in pdg_colors])) == 1
        for ipart, (pdg, col) in enumerate(pdg_colors):
            if pdg in pdg_dict and col in [3, 8] and \
                   (not octets_only or ipart < len(pdg_colors) - 1):
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
        self.set('interactions', interactions)
        return 

