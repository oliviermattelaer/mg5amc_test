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
 
"""Classes for generation of color-ordered diagrams, and for
generating matrix elements using Behrends-Giele
currents.

ColorOrderedAmplitude keeps track of all color flows, and
ColorOrderedFlow performs the diagram generation. ColorOrderedLeg has
extra needed flags. See documentation for the functions
get_combined_legs and get_combined_vertices for the algorithm used to
generate only color-ordered diagrams.

ColorOrderedModel adds color singlets coupling only to triplets for
each color octet in the model.

BGHelasMatrixElement generates the matrix element for a given color
flow, by ensuring that all wavefunctions and amplitudes correspond to
the color flow defined by the process color ordering. If optimization
is turned on, wavefunctions combining the same external particles are
combined into BGHelasCurrents, to reduce the number of wavefunction
and amplitude calculations needed to the Behrends-Giele n^3 (or n^4 if
4-gluon vertices are present), or 3^n for color singlet processes.
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
import madgraph.core.helas_objects as helas_objects
from madgraph import MadGraph5Error

logger = logging.getLogger('madgraph.color_ordered_amplitudes')

# Dictionary from spin number to name used in the wavefunction sum functions
spin_dict = {1:'S', 2:'F', 3:'V', 4:'T'}

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
        legs = base_objects.LegList([copy.copy(l) for l in \
                                     process.get('legs')])

        if not isinstance(process.get('model'), ColorOrderedModel):
            process.set('model', ColorOrderedModel(process.get('model')))

        model = process.get('model')
        
        # Add color negative singlets to model corresponding to all
        # color octets, with only 3-3bar-8 interactions.
        # Color factor: -1/(2*Nc) * Id(1,2)

        # Set leg numbers for the process
        for i, leg in enumerate(legs):
            leg.set('number', i+1)
            # Reverse initial state legs to have all legs outgoing
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
        """Allow generating a HelasWavefunction from a Leg
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

#===============================================================================
# COHelasWavefunction
#===============================================================================
class COHelasWavefunction(helas_objects.HelasWavefunction):
    """COHelasWavefunction object, a HelasWavefunction with added
    fields for comparison
    """

    def default_setup(self):
        """Default values for all properties"""

        super(COHelasWavefunction, self).default_setup()
        # Comparison array, with mother external numbers, interaction
        # id and lorentz-color index
        self['compare_array'] = []
        # Comparison array to find wavefunctions that can be combined
        # into the same current, with external numbers, pdg code and
        # fermion flow state
        self['current_array'] = []
        # external wavefunction numbers
        self['external_numbers'] = array.array('I')
        # Information for color calculation
        self['color_string'] = color.ColorString()
        # Unique number for color substitution
        self['lastleg_number'] = 0
        # Factor for wavefunction in wf summation in form
        # (fraction, is_imaginary?)
        self['factor'] = (1, fractions.Fraction(1,1), False)
        # Fermion external numbers (for fermion exchange factor)
        self['external_fermion_numbers'] = []
        # Complete set of orders for the wavefunctions up to this
        self['all_orders'] = None

    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name in ['compare_array', 'current_array']:
            if not isinstance(value, list) and not \
                   isinstance(value, array.array):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list object" % str(value)
        elif name == 'external_numbers':
            if not isinstance(value, array.array):
                raise self.PhysicsObjectError, \
                        "%s is not a valid array object" % str(value)
        elif name == 'color_string':
            if not isinstance(value, color.ColorString):
                raise self.PhysicsObjectError, \
                        "%s is not a valid ColorString object" % str(value)
        elif name == 'factor':
            if not isinstance(value, tuple):
                raise self.PhysicsObjectError, \
                        "%s is not a valid tuple" % str(value)
        elif name == 'external_fermion_numbers':
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list" % str(value)
        elif name == 'all_orders':
            if not isinstance(value, dict) or value == None:
                raise self.PhysicsObjectError, \
                        "%s is not a valid dict" % str(value)
        else:
            super(COHelasWavefunction, self).filter(name, value)

        return True

    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""

        return super(COHelasWavefunction, self).get_sorted_keys() + \
               ['compare_array', 'current_array', 'external_numbers',
                'color_string', 'lastleg_number', 'factor']

    def get(self, name):
        """Enhanced get function to initialize compare_array,
        external_numbers and factor."""

        if name in ['compare_array', 'current_array', 'external_numbers'] and \
               not self[name]:
            self.create_arrays()
        if name == 'external_fermion_numbers' and not self[name]:
            self.set_external_fermion_numbers()
        if name == 'all_orders' and self[name] == None:
            self.set_all_orders()
            
        return super(COHelasWavefunction, self).get(name)        

    def create_arrays(self):
        """Create the comparison arrays compare_array, current_array and
        external_numbers for a COHelasWavefunction.

        external_numbers is the sorted set of external particle
        numbers that are included in this wavefunction.

        current_array is used to determine if two wavefunctions belong
        in the same current, i.e. if they have the same external
        numbers, pdg code, fermion flow state and interaction orders.
        
        compare_array is the tuple of [mothers current_arrays,
        interaction id, coupling key], used to find duplicate
        wavefunctions (after replacing mothers with currents).
        """

        if not self.get('mothers'):
            # This is an external wavefunction
            self.set('external_numbers',
                     array.array('I',[self.get('number_external')]))
            self.set('compare_array', self['external_numbers'])
            self.set('current_array', self['external_numbers'])
        else:
            self.set('external_numbers',
                     array.array('I',
                                 sorted(sum([list(m.get('external_numbers')) \
                                        for m in self.get('mothers')], []))))
            self.set('current_array', [self['external_numbers'],
                                       self.get('pdg_code'),
                                       self.get_with_flow('state'),
                                       sorted(self.get('all_orders').items())])
            self.set('compare_array', [[m.get('current_array') for \
                                        m in self.get('mothers')],
                                       self.get('interaction_id'),
                                       self.get('color_key')])

    def set_all_orders(self):
        """Set the dictionary of orders from this wf and all mothers"""

        all_orders = copy.copy(self.get('orders'))
        for mother in self.get('mothers'):
            for order in mother.get('all_orders'):
                if order in all_orders:
                    all_orders[order] += mother.get('all_orders')[order]
                else:
                    all_orders[order] = mother.get('all_orders')[order]
        self['all_orders'] = all_orders

    def set_color_and_fermion_factor(self):
        """Set the color and fermion factor for this wavefunction. The
        factor is of the type
        (fermion_factor, color_coeff, is_imaginary?)"""

        if self.get('color_string'):
            self['factor'] = (self.calculate_fermionfactor(),
                              self.get('color_string').coeff,
                              self.get('color_string').is_imaginary)        

    def set_external_fermion_numbers(self):
        """Set the fermion number used for fermionfactor calculation"""

        if self.is_fermion() and not self.get('mothers'):
            self['external_fermion_numbers'] = [self['number_external']]
        else:
            self['external_fermion_numbers'] = \
                            sum([sorted(wf.get('external_fermion_numbers')) for \
                                 wf in self.get('mothers')], [])

    def calculate_fermionfactor(self):
        """Calculate the fermion factor (needed sign flips for mother
        fermions, if this wavefunction has a pair of fermion mothers,
        times the product of fermion factors of the mothers).
        Sort external_fermion_numbers."""

        # Pick out fermion mothers
        
        return helas_objects.HelasAmplitude.sign_flips_to_order(\
            self.get('external_fermion_numbers')) * \
            reduce(lambda x1, x2: x1*x2, [m.get('factor')[0] for m in \
                                          self.get('mothers')])
        
    def get_base_vertices(self, wf_dict, vx_list = [], optimization = 1):
        """Recursive method to get a list of base_objects.VertexList,
        corresponding to this wavefunction and its mothers."""

        vertices = base_objects.VertexList()
        vertex_lists = [vertices]
        
        mothers = self.get('mothers')

        if not mothers:
            return vertex_lists

        # Add vertices for all mothers
        for mother in mothers:
            # This is where recursion happens
            mother_vertices = mother.get_base_vertices(\
                wf_dict, vx_list,optimization)
            new_vertex_lists = []
            for vertex in mother_vertices:
                new_vertex_list = [copy.copy(v) for v in vertex_lists]
                for v in new_vertex_list:
                    v.extend(vertex)
                new_vertex_lists.extend(new_vertex_list)
            vertex_lists = new_vertex_lists

        vertex = self.get_base_vertex(wf_dict, vx_list, optimization)

        try:
            index = vx_list.index(vertex)
            vertex = vx_list[index]
        except ValueError:
            pass
        
        for vertices in vertex_lists:
            vertices.append(vertex)

        return vertex_lists

    def get_wf_numbers(self):
        """Recursive function returning the set of wavefunction
        numbers that correspond to this wavefunction and its mothers."""

        numbers = sum([list(m.get_wf_numbers()) for m in self.get('mothers')], [])
        numbers.append(self.get('number'))
        return set(numbers)

#===============================================================================
# COHelasAmplitude
#===============================================================================
class COHelasAmplitude(helas_objects.HelasAmplitude):
    """COHelasAmplitude object, a HelasAmplitude with added
    fields for comparison and keeping track of color
    """

    def default_setup(self):
        """Default values for all properties"""

        super(COHelasAmplitude, self).default_setup()
        # external wavefunction numbers
        self['compare_array'] = []
        # Color string
        self['color_string'] = color.ColorString()
        # Factor for amplitude in JAMP summation in form
        # (fraction, is_imaginary?)
        self['factor'] = (1, fractions.Fraction(1,1), False)
        # Complete set of orders for this diagram
        self['all_orders'] = None

    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name == 'compare_array':
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list object" % str(value)
        elif name == 'color_string':
            if not isinstance(value, color.ColorString):
                raise self.PhysicsObjectError, \
                        "%s is not a valid ColorString object" % str(value)
        elif name == 'factor':
            if not isinstance(value, tuple):
                raise self.PhysicsObjectError, \
                        "%s is not a valid tuple" % str(value)
        elif name == 'all_orders':
            if not isinstance(value, dict) or value == None:
                raise self.PhysicsObjectError, \
                        "%s is not a valid dict" % str(value)
        else:
            super(COHelasAmplitude, self).filter(name, value)

        return True

    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""

        return super(COHelasAmplitude, self).get_sorted_keys() + \
               ['compare_array', 'color_string']

    def get(self, name):
        """Enhanced get function to initialize compare_array."""

        if name in ['compare_array'] and not self[name]:
            self.create_arrays()

        if name == 'all_orders' and self[name] == None:
            self.set_all_orders()
            
        return super(COHelasAmplitude, self).get(name)        

    def set(self, name, value):
        """Enhanced set function to correctly set color order and
        factor for singlet_QCD couplings - we need one factor in Nc
        for every two singlet_QCD couplings."""

        super(COHelasAmplitude, self).set(name, value)

        if name == 'color_string' and \
                      'singlet_QCD' in self.get('all_orders') and \
                      self.get('all_orders')['singlet_QCD'] > 0:
            self.set_singlet_color_string()
            
    def set_all_orders(self):
        """Set the dictionary of orders from this amp and all mothers"""

        all_orders = copy.copy(self.get('orders'))
        for mother in self.get('mothers'):
            for order in mother.get('all_orders'):
                if order in all_orders:
                    all_orders[order] += mother.get('all_orders')[order]
                else:
                    all_orders[order] = mother.get('all_orders')[order]
        self['all_orders'] = all_orders
        
    def set_singlet_color_string(self):
        """Modifiy the color factor due to singlet_QCD couplings - we
        need one factor in Nc and a coefficient -1/2 for every pair of
        singlet_QCD couplings."""

        color_string = self.get('color_string')
        singlet_orders = self.get('all_orders')['singlet_QCD'] / 2
        color_string.coeff *= fractions.Fraction(-1, 2)**singlet_orders
        color_string.Nc_power += singlet_orders

    def create_arrays(self):
        """Create the comparison array compare_array, to find
        identical amplitudes after replacing mothers with currents:
        [mother current arrays, interaction id, coupling key]."""

        self.set('compare_array', [sorted([list(m.get('current_array')) \
                                              for m in self.get('mothers')]),
                                      self.get('interaction_id'),
                                      self.get('color_key')])

    def set_color_and_fermion_factor(self):
        """Set the color and fermion factor for this wavefunction"""

        if self.get('color_string'):
            self['factor'] = (self.calculate_fermionfactor(),
                              self.get('color_string').coeff,
                              self.get('color_string').is_imaginary)
            
    def calculate_fermionfactor(self):
        """Calculate the fermion factor (needed sign flips for mother
        fermions, if this amplitude has multiple fermion mothers,
        times the product of fermion factors of the mothers)."""

        # Pick out fermion mothers
        fermion_numbers = sum([sorted(wf.get('external_fermion_numbers')) \
                                for wf in self.get('mothers')], [])

        return self.sign_flips_to_order(fermion_numbers) * \
            reduce(lambda x1, x2: x1*x2, [m.get('factor')[0] for m in \
                                          self.get('mothers')])
        
    def get_base_diagram(self, wf_dict, vx_list = [], optimization = 1):
        """Return the list of base_objects.Diagram which corresponds to this
        amplitude, using a recursive method for the wavefunctions."""

        vertices = base_objects.VertexList()
        vertex_lists = [vertices]

        # Add vertices for all mothers
        for mother in self.get('mothers'):
            # This is where recursion happens
            mother_vertices = mother.get_base_vertices(\
                wf_dict, vx_list,optimization)
            new_vertex_lists = []
            for vertex in mother_vertices:
                new_vertex_list = [copy.copy(v) for v in vertex_lists]
                for v in new_vertex_list:
                    v.extend(vertex)
                new_vertex_lists.extend(new_vertex_list)
            vertex_lists = new_vertex_lists

        # Generate last vertex
        vertex = self.get_base_vertex(wf_dict, vx_list, optimization)

        diagrams = base_objects.DiagramList()
        for vertices in vertex_lists:
            vertices.append(vertex)
            diagrams.append(base_objects.Diagram({'vertices': vertices}))

        return diagrams

    def get_wf_numbers(self):
        """Recursive function returning the set of wavefunction
        numbers that correspond to this wavefunction and its mothers."""

        numbers = sum([list(m.get_wf_numbers()) for m in self.get('mothers')], [])
        return set(numbers)

#===============================================================================
# BGHelasCurrent
#===============================================================================
class BGHelasCurrent(COHelasWavefunction):
    """BGHelasCurrent object, which combines HelasWavefunctions into
    Behrends-Giele currents. The color and fermion factor for each
    wavefunction is taken into account at the time of wavefunction
    combination, so the corresponding factors are reset for the
    current.
    """

    # Customized constructor
    def __init__(self, *arguments):
        """Correctly initialize a BGHelasCurrent with a COHelasWavefunction
        """

        if len(arguments) == 1:
            if isinstance(arguments[0], COHelasWavefunction) and not \
                   isinstance(arguments[0], BGHelasCurrent):
                super(BGHelasCurrent, self).__init__(arguments[0])
                self.set('mothers', helas_objects.HelasWavefunctionList())
                self.set('compare_array', [])
                self.set('external_numbers', array.array('I'))
                self.set('external_fermion_numbers',
                         sorted(arguments[0].get('external_fermion_numbers')))
                self.set('all_orders', arguments[0].get('all_orders'))
            else:
                super(BGHelasCurrent, self).__init__(*arguments)
        else:
            super(BGHelasCurrent, self).__init__(*arguments)

    def get(self, name):
        """Special get for lorentz"""

        if name == 'lorentz':
            self.set('lorentz', ['sum%s%d' % \
                     (spin_dict[self.get('mothers')[0].get('spin')],
                      len(self.get('mothers')))])
        return super(BGHelasCurrent, self).get(name)

    def set_color_string(self):
        """Set color string based on the first mother"""
        col_str = copy.copy(self.get('mothers')[0].get('color_string'))
        col_str.coeff = fractions.Fraction(1)
        col_str.is_imaginary = False
        self.set('color_string', col_str)
        # Since color and fermion factors are already included for
        # each wavefunction, we should set the factor to 1 for current
        self.set('factor', (1, fractions.Fraction(1,1), False))
        # Set Nc power for the combined wavefunction based on max Nc
        # power of wfs
        common_Nc_power = max([wf.get('color_string').Nc_power for wf \
                              in self.get('mothers')])
        self.get('color_string').Nc_power = common_Nc_power
        
    def get_call_key(self):
        """Generate the ('sum', spins) tuple used as key for
        the helas call dictionaries in HelasCallWriter"""

        res = [m.get('spin') for m in self.get('mothers')]

        return ('sum', tuple(res))

    def find_outgoing_number(self):
        """Return -1, to indicate that this is a special wavefunction."""

        return -1

    def get_base_vertices(self, wf_dict, vx_list = [], optimization = 1):
        """Returns a list of base_objects.VertexList,
        corresponding to the wavefunctions of the mothers."""

        mothers = self.get('mothers')

        vertex_lists = []

        # Add vertices for all mothers
        for mother in mothers:
            # This is where recursion happens
            vertices = mother.get_base_vertices(\
                wf_dict, vx_list,optimization)
            for vertex in vertices:
                if not vertex in vertex_lists:
                    vertex_lists.append(vertex)

        return vertex_lists

#===============================================================================
# COHelasFlow
#===============================================================================
class COHelasFlow(helas_objects.HelasMatrixElement):
    """COHelasFlow: Color ordered version of a
    HelasMatrixElement, starting from a ColorOrderedFlow.

    After regular helas diagram generation, performs the following:

    - If no BG optimization, go through all wavefunctions and
    amplitudes to select only those with the correct color structure

    If BG optimization:

    1. Go through all wavefunctions to combine the wavefunctions with
    the same external particle numbers using BGHelasCurrents, keeping
    only those wfs with the correct color structure

    2. Go through all amplitudes, and use only one amplitude per set
    of BGHelasCurrents (keeping only amplitudes with the correct color
    structure).
    """

    # Keep unique number for color substitution
    lastleg_number = 0

    def default_setup(self):
        """Add list of leg permutations and color string to HelasMatrixElement"""

        super(COHelasFlow, self).default_setup()

        # Permutations of external legs which give valid color flows.
        # The first permutation corresponds to the this color flow.
        self['permutations'] = []
        # Color string corresponding to this color flow
        self['color_string'] = color.ColorString()
        # Identifier for this color flow
        self['number'] = 0

    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name == 'permutations':
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list object" % str(value)
        elif name == 'color_string':
            if not isinstance(value, color.ColorString):
                raise self.PhysicsObjectError, \
                        "%s is not a valid color string" % str(value)
        elif name in ['number']:
            if not isinstance(value, int):
                raise self.PhysicsObjectError, \
                        "%s is not a valid integer" % str(value)
        else:
            super(COHelasFlow, self).filter(name, value)

        return True
    
    # Customized constructor
    def __init__(self, *arguments, **opts):
        """Add number to arguments of HelasMatrixElement
        """

        number = 0
        if 'number' in opts:
            number = opts['number']
            del opts['number']
        super(COHelasFlow, self).__init__(*arguments, **opts)        
        self.set('number', number)

    def generate_helas_diagrams(self, amplitude, optimization=3,
                                decay_ids=[]):
        """Generate Behrends-Giele diagrams for a color ordered amplitude
        """
        
        assert  isinstance(amplitude, ColorOrderedFlow), \
                    "Missing or erraneous arguments for generate_helas_diagrams"

        use_bg_currents = (optimization // 2 == 1)

        # Set permutations
        self.set('permutations', amplitude.get('permutations'))
        
        # Set color string
        self.set('color_string', self.get_color_string(range(1,
                                 len(amplitude.get('process').get('legs'))+1)))

        # First generate full set of wavefunctions and amplitudes
        super(COHelasFlow, self).generate_helas_diagrams(amplitude,
                                                         optimization,
                                                         decay_ids)

        # Collect all wavefunctions and turn into COHelasWavefunctions
        all_wavefunctions = self.get_all_wavefunctions()
        co_wavefunctions = helas_objects.HelasWavefunctionList(\
            [COHelasWavefunction(wf) for wf in all_wavefunctions])

        # Same thing for amplitudes
        for diag in self.get('diagrams'):
            diag.set('amplitudes', helas_objects.HelasAmplitudeList(\
                [COHelasAmplitude(amp) for amp in diag.get('amplitudes')]))

        # current_arrays for combined wfs
        combined_wf_arrays = []
        # lists of combined wfs with the same current_array
        combined_wavefunctions = []
        # list of compare_arrays for the wavefunctions in the currents
        compare_arrays = []
        # List of wf numbers that failed due to color
        removed_wfs = []
        # Dict from old wf number to new wf/BG current
        wf_co_dict = dict((i+1,wf) for i, wf in \
                               enumerate(co_wavefunctions))
        # Go through wavefunctions and combine into currents
        for wf in co_wavefunctions:
            # Replace mothers for this wf with CO wavefunction
            self.replace_mothers(wf, wf_co_dict)

            if use_bg_currents:
                # Combine wavefunctions with the same current_array
                wf_array = wf.get('current_array')                
                try:
                    index = combined_wf_arrays.index(wf_array)
                except ValueError:
                    combined_wf_arrays.append(wf_array)
                    combined_wavefunctions.append(\
                        helas_objects.HelasWavefunctionList([wf]))
                else:
                    combined_wavefunctions[index].append(wf)
            else:
                combined_wavefunctions.append(\
                    helas_objects.HelasWavefunctionList([wf]))

        # Sort currents according to length (number of wfs) to ensure
        # that all mothers come before their daughters.
        if use_bg_currents:
            combined_wavefunctions = sorted(combined_wavefunctions,
                                            lambda l1,l2:len(l1)-len(l2))

        # Dict from old wf number to new wf/BG current
        wf_current_dict = {}

        # Go through combined wavefunction and replace with currents
        combined_currents = helas_objects.HelasWavefunctionList()
        for wf_list in combined_wavefunctions:
            # The final mothers for this current
            bg_mothers = helas_objects.HelasWavefunctionList()
            # Numbers of wfs in this current
            wf_numbers = []
            while wf_list:
                wf = wf_list.pop(0)
                # Check if this wf has any removed mothers
                if any([m.get('number') in removed_wfs for m in \
                        wf.get('mothers')]):
                    removed_wfs.append(wf.get('number'))
                    continue

                if use_bg_currents:
                    # Replace the wavefunction mothers in this
                    # wavefunction with corresponding currents
                    self.replace_mothers(wf, wf_current_dict)

                # Check if color is incorrect
                if not self.check_color(wf):
                    removed_wfs.append(wf.get('number'))
                    continue

                # This is a valid wf, so append number to wfs
                wf_numbers.append(wf.get('number'))

                # Check if an identical wavefunction (after
                # replacing mothers with currents) is already
                # present in the current (if using BG currents)
                if use_bg_currents and wf.get('compare_array') in \
                       [m.get('compare_array') for m in bg_mothers]:
                    continue

                # Add the resulting wavefunction to
                # combined_currents and to bg_mothers
                combined_currents.append(wf)
                bg_mothers.append(wf)
            
            if len(bg_mothers) > 1:
                # Combine wavefunctions to a current
                combine_wf = BGHelasCurrent(bg_mothers[0])
                # Set the mothers
                combine_wf.set('mothers', bg_mothers)
                # Set color string for the new BG current
                combine_wf.set_color_string()
                # Add combine_wf to combined_currents
                combined_currents.append(combine_wf)
                # Add this BG current to the replacement dictionary
                for wf_number in wf_numbers:
                    wf_current_dict[wf_number] = combine_wf
            elif bg_mothers:
                for wf_number in wf_numbers:
                    wf_current_dict[wf_number] = bg_mothers[0]

        # left_diagrams is the diagrams that are left after BG
        # combinations
        left_diagrams = helas_objects.HelasDiagramList()
        idiag = 0
        amp_compare_arrays = []
        for diagram in self.get('diagrams'):
            valid_amps = helas_objects.HelasAmplitudeList()
            for amp in diagram.get('amplitudes'):
                # Replace wf mothers in this diagram
                self.replace_mothers(amp, wf_current_dict)
                # Check that this amp has any removed mothers
                if any([m.get('number') in removed_wfs for m in \
                        amp.get('mothers')]):
                    continue
                # Check if an amplitude with identical mothers,
                # interaction id and coupling key is already present -
                # in that case skip this amp
                if amp.get('compare_array') in amp_compare_arrays:
                    continue
                # Check that this amp passes color check
                if not self.check_color(amp):
                    continue
                # Update color and fermion factor, since mothers have changed
                amp.set_color_and_fermion_factor()
                # Otherwise add this diagram to the list
                valid_amps.append(amp)
                # If we use BG currents, add compare_array
                if use_bg_currents:
                    amp_compare_arrays.append(amp.get('compare_array'))
            # If no valid amplitudes, skip this diagram
            if not valid_amps:
                continue
            # Otherwise replace amplitudes and add to list
            left_diagrams.append(diagram)
            diagram.set('amplitudes', valid_amps)
            diagram.set('wavefunctions', helas_objects.HelasWavefunctionList())

        # Set diagram numbers
        for i,d in enumerate(left_diagrams):
            d.set('number', i+1)
        self.set('diagrams', left_diagrams)

        # Set wf number for all wavefunctions
        for i, wf in enumerate(combined_currents):
            wf.set('number', i+1)

        # Set amplitude number for all amplitudes
        for i, amp in enumerate(self.get_all_amplitudes()):
            amp.set('number', i+1)

        # Set wavefunctions in first diagram
        self.get('diagrams')[0].set('wavefunctions', combined_currents)

        # Set Nc power for the overall color string instead of every
        # amplitude, but only if it is negative
        common_Nc_power = 0
        if self.get_all_amplitudes():
            common_Nc_power = min(max([amp.get('color_string').Nc_power for \
                                       amp in self.get_all_amplitudes()]), 0)
        self.get('color_string').Nc_power = common_Nc_power
        for amp in self.get_all_amplitudes():
            amp.get('color_string').Nc_power -= common_Nc_power

    def get_color_amplitudes(self):
        """Return a list of (coefficient, amplitude number) lists,
        corresponding to the JAMPs for this matrix element. The
        coefficients are given in the format (fermion factor, color
        coeff (frac), imaginary, Nc power)."""

        col_amp = []
        for amplitude in self.get_all_amplitudes():
            col_amp.append((amplitude.get('factor') + 
                            (amplitude.get('color_string').Nc_power,),
                            amplitude.get('number')))
        return [col_amp]

    def replace_mothers(self, arg, wf_current_dict):
        """Find BG currents in combined_wavefunctions corresponding to
        the mothers of arg (wavefunction or amplitude), replace
        mothers with BG currents"""

        for i, m in enumerate(arg.get('mothers')):
            try:
                arg.get('mothers')[i] = wf_current_dict[m.get('number')]
            except:
                pass

    def check_color(self, arg):
        """Check that the wavefunction/amplitude color is consistent
        with the color ordering, and set the corresponding color of
        the wavefunction or amplitude."""

        assert isinstance(arg, COHelasWavefunction) or \
               isinstance(arg, COHelasAmplitude)

        if not arg.get('mothers'):
            return True

        process = self.get('processes')[0]
        model = process.get('model')
        # Create a Vertex that we can use to extract the color for
        # this wavefunction
        base_vertex = arg.get_base_vertex({})
        # Replace leg numbers in the base_vertex with the
        # lastleg_numbers of the mother wavefunctions to get correct
        # color strings
        for i, mother in enumerate(arg.get('mothers')):
            if mother.get('lastleg_number'):
                base_vertex.get('legs')[i].set('number',
                                               mother.get('lastleg_number'))
        lastleg = None
        if isinstance(arg, COHelasWavefunction):
            # We need to take care of the last leg, by giving it a
            # unique number
            self.lastleg_number -= 1
            lastleg = ColorOrderedLeg(base_vertex.get('legs').pop(-1))
            lastleg.set('id', model.get_particle(lastleg.get('id')).\
                        get_anti_pdg_code())
            lastleg.set('number', self.lastleg_number)
            base_vertex.get('legs').insert(0, lastleg)
            arg.set('lastleg_number', self.lastleg_number)
            color_indices = list(arg.get('external_numbers'))
            # Get the color string that we need to compare to
            comp_color_string = self.get_color_string(color_indices,
                                                      lastleg)
        else:
            comp_color_string = self.get('color_string')
            
        # Prepare for extracting the color dict for this vertex using
        # ColorBasis.add_vertex
        base_diagram = base_objects.Diagram({"vertices": \
                                    base_objects.VertexList([base_vertex])})
        col_basis = color_amp.ColorBasis()
        # Now extract the color dict for this vertex
        min_index, color_dict = col_basis.add_vertex(base_vertex,
                                                     base_diagram,
                                                     model,
                                                     {}, {},
                                                     self.lastleg_number - 1)
        # Pick out only the relevant color string for arg
        color_string = color_dict[(arg.get('color_key'),)]
        # Add the color strings of all mothers
        for mother in arg.get('mothers'):
            if mother.get('color_string'):
                color_string.product(mother.get('color_string'))

        # Now simplify color string, and check if we have a
        # contribution corresponding to comp_color_string
        col_fact = color.ColorFactor([color_string])
        col_fact = col_fact.full_simplify()
        similar_strings = [cs for cs in col_fact if \
                           cs.to_canonical() == comp_color_string.to_canonical()]
        if not similar_strings:
            return False
        assert(len(similar_strings) == 1)
        arg.set('color_string', similar_strings[0])
        # Set color and fermion factor once and for all
        arg.set_color_and_fermion_factor()
        return True

    def get_color_string(self, external_numbers, lastleg = None):
        """Get the color string corresponding to the given external
        numbers, using the color ordering flags in the process, and
        insert lastleg as appropriate (for wavefunctions)."""

        legs = self.get('processes')[0].get('legs')
        model = self.get('processes')[0].get('model')
        color_chains = {}
        chain_types = {}
        leg_number_dict = {}
        # Determine all participating color chains based on color
        # ordering flags in process
        for number in external_numbers:
            leg = [l for l in legs if l.get('number') == number][0]
            if not leg.get('color_ordering'):
                continue
            co = leg.get('color_ordering').keys()[0]
            co_val = leg.get('color_ordering').values()[0][0]
            try:
                color_chains[co].append(co_val)
            except:
                color_chains[co] = [co_val]
                chain_types[co] = 0
            leg_number_dict[(co, co_val)] = leg.get('number')
            # Set color (for final state particles) or anticolor (for IS)
            if leg.get('state'):
                leg_color = model.get_particle(leg.get('id')).get_color()
            else:
                leg_color = model.get_particle(leg.get('id')).get_anti_color()
            # Set chain_types to: 2 for chain with both triplet and
            # antitriplet, 1 for chain with only triplet, or -1 for
            # chain with only antitriplet
            if abs(leg_color) == 3:
                if chain_types[co]:
                    chain_types[co] = 2
                else:
                    chain_types[co] = leg_color//3

        color_string = color.ColorString()
        # Insert last leg in appropriate place
        if lastleg:
            lastleg_color = model.get_particle(lastleg.get('id')).get_color()
            if abs(lastleg_color) == 3:
                # Add to chain missing a 3 or 3bar
                key = [k for k in chain_types.keys() if \
                             abs(chain_types[k]) == 1][0]
                # Pick out all legs in the process contributing to
                # this chain
                leg_orderings = [l.get('color_ordering')[key][0] for l in \
                                 legs if key in l.get('color_ordering')]
                chain_types[key] = 2
                if lastleg_color == 3:
                    # A 3 should be put first in the chain
                    lastleg_order = min(leg_orderings)
                else:
                    # A 3bar should be put last in the chain
                    lastleg_order = max(leg_orderings)
                color_chains[key].append(lastleg_order)
                leg_number_dict[(key, lastleg_order)] = \
                                                  lastleg.get('number')
            elif abs(lastleg_color) == 8:
                # Add leg to all color chains that are not yet completed
                for key in color_chains.keys():
                    # Pick out all legs in the process contributing to
                    # this chain
                    leg_orderings = [l.get('color_ordering')[key][0] for l in \
                                     legs if key in l.get('color_ordering')]
                    if len(color_chains[key]) != len(leg_orderings):
                        color_chains[key].sort()
                        # If color chain is unfinished triplet
                        # (antitriplet) chain, place gluon at place of
                        # antitriplet (triplet)
                        if chain_types[key] == 1:
                            
                            lastleg_order = max(leg_orderings)
                        elif chain_types[key] == -1:
                            lastleg_order = min(leg_orderings)
                        # Otherwise, find gap in color chain
                        elif max(color_chains[key]) < len(leg_orderings):
                            lastleg_order = max(color_chains[key]) + 1
                        else:
                            lastleg_order = max([c for c in leg_orderings if \
                                                 c not in color_chains[key]])
                        color_chains[key].append(lastleg_order)
                        leg_number_dict[(key, lastleg_order)] = \
                                                         lastleg.get('number')

        for key in color_chains:
            # Order entries according to color chain order (3bar,8,...,3)
            color_chains[key].sort()
            # Replace with leg numbers
            color_chains[key] = [leg_number_dict[(key, co_val)] for \
                                                 co_val in color_chains[key]]
            
        # If we have triplet-antitriplet chains, combine them into a
        # common T chain
        triplet_chain_keys = [key for key in color_chains.keys() if \
                              abs(chain_types[key]) == 1]
        if triplet_chain_keys:
            assert len(triplet_chain_keys) == 2
            triplet_chain_key = [key for key in color_chains.keys() if \
                                 chain_types[key] == 1][0]
            antitriplet_chain_key = [key for key in color_chains.keys() if \
                                     chain_types[key] == -1][0]
            color_chains[triplet_chain_key].extend(\
                color_chains[antitriplet_chain_key][1:])
            chain_types[triplet_chain_key] = 2
            del color_chains[antitriplet_chain_key]
            del chain_types[antitriplet_chain_key]
            
        # Create color string based on color_chains
        for cckey in color_chains.keys():
            if chain_types[cckey] == 2:
                # Move triplet from first to next-to-last position
                color_chains[cckey].insert(-1,color_chains[cckey].pop(0))
                # A 88...33bar chain is a T
                color_string.append(color.T(*color_chains[cckey]))
            else:
                # Fix ordering to have lastleg first (since ordered
                # with minimum leg number first, and lastleg has
                # negative number)
                if lastleg:
                    chain = color_chains[cckey]
                    number = lastleg.get('number')
                    color_chains[cckey] = chain[chain.index(number):] + \
                                          chain[:chain.index(number)]
                # A 88..8 chain is a Tr
                color_string.append(color.Tr(*color_chains[cckey]))
                
        return color_string
    
    def get_base_amplitude(self):
        """Generate a diagram_generation.Amplitude from a
        HelasCOFlow. Note that this is not supposed to be used for
        color generation."""

        # Need to take care of diagram numbering for decay chains
        # before this can be used for those!

        optimization = 1
        if len(filter(lambda wf: wf.get('number') == 1,
                      self.get_all_wavefunctions())) > 1:
            optimization = 0

        model = self.get('processes')[0].get('model')

        wf_dict = {}
        vx_list = []
        diagrams = base_objects.DiagramList()
        for diag in self.get('diagrams'):
            if diag.get('amplitudes'):
                diagrams.extend(diag.get('amplitudes')[0].get_base_diagram(\
                                             wf_dict, vx_list, optimization))

        for diag in diagrams:
            diag.calculate_orders(self.get('processes')[0].get('model'))
            
        return diagram_generation.Amplitude({\
            'process': self.get('processes')[0],
            'diagrams': diagrams})

#===============================================================================
# COHelasFlowList
#===============================================================================
class COHelasFlowList(base_objects.PhysicsObjectList):
    """List of COHelasFlow objects
    """

    def is_valid_element(self, obj):
        """Test if object obj is a valid COHelasFlow for the list."""

        return isinstance(obj, COHelasFlow)

#===============================================================================
# COHelasMatrixElement
#===============================================================================
class COHelasMatrixElement(helas_objects.HelasMatrixElement):
    """COHelasMatrixElement: Holds a set of COHelasFlows and takes care
    of the color matrix calculation"""

    
    def default_setup(self):
        """COHelasMatrixElement has a list of COHelasFlows.

        Note that the color matrix is redefined so that it has rows
        corresponding to each flow, and columns corresponding to all
        permutations of all flows."""


        super(COHelasMatrixElement, self).default_setup()

        self['color_flows'] = COHelasFlowList()
        self['min_Nc_power'] = 0
        self['permutations'] = []

    def filter(self, name, value):
        """Filter for valid diagram property values."""

        if name == 'color_flows':
            if not isinstance(value, COHelasFlowList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid COHelasFlowList object" % str(value)
        
        elif name == 'min_Nc_power':
            if not isinstance(value, int):
                raise self.PhysicsObjectError, \
                        "%s is not a valid integer" % str(value)

        elif name == 'permutations':
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list" % str(value)
        
        else:
            super(COHelasMatrixElement, self).filter(name, value)

        return True
    
    def __init__(self, amplitude = None, optimization = 3, decay_ids = [],
                 gen_color = 2):
        """Initialize a COHelasMatrixElement with a ColorOrderedAmplitude"""
        
        if amplitude != None:
            if isinstance(amplitude, ColorOrderedAmplitude):
                super(COHelasMatrixElement, self).__init__()
                self.get('processes').append(amplitude.get('process'))
                self.set('has_mirror_process',
                         amplitude.get('has_mirror_process'))
                for iflow, flow in enumerate(amplitude.get('color_flows')):
                    self.get('color_flows').append(COHelasFlow(flow,
                                                   optimization, decay_ids,
                                                   gen_color = False,
                                                   number = iflow + 1))
                if self.get('color_flows'):
                    self.set('permutations',
                             self.get('color_flows')[0].get('permutations'))
                if gen_color and not self.get('color_matrix'):
                    self.build_color_matrix(gen_color)
            else:
                # In this case, try to use amplitude as a dictionary
                super(COHelasMatrixElement, self).__init__(amplitude)
        else:
            super(COHelasMatrixElement, self).__init__()
        

    def build_color_matrix(self, gen_color):
        """Build the relevant lines of the color matrix, to the order
        in color given in gen_color (0 = none, 1 = leading, 
        2 = next-to-leading, etc)"""
        
        if not gen_color:
            return

        col_matrix = self.get('color_matrix')

        if col_matrix:
            return

        # Build the color matrix based on the color_basis elements
        # corresponding to all permutations of all flows
        # but only include non-zero entries (to the wanted color order).
        # The resulting color matrix will only have for non-zero
        # entries, and only lines for the basic color flows (first
        # permutation)

        # First figure out maximum Nc power from the color strings
        Nc_powers = []
        for color_flow in self.get('color_flows'):
            col_str = color_flow.get('color_string').create_copy()
            col_str2 = col_str.create_copy()
            # Multiply color string with its complex conjugate
            col_str.product(col_str2.complex_conjugate())
            # Create a color factor to store the result and simplify it
            col_fact = color.ColorFactor([col_str])
            result = col_fact.full_simplify()
            # result now has all Nc powers
            Nc_powers.extend([cs.Nc_power for cs in result])

        # Set min Nc power to max(Nc_powers) - (gen_color - 1),
        # i.e., Nc^Nmax for leading, Nc^(Nmax-1) to subleading, etc.
        min_Nc_power = 0
        if Nc_powers:
            min_Nc_power = max(Nc_powers) - (gen_color - 1)
            self.set('min_Nc_power', min_Nc_power)
        # The col_basis will contain the basis strings for the basic flows
        basic_immutable_factors = []
        # canonical_dict to store previous calculation results
        canonical_dict = {}
        if gen_color == 1:
            permutations = self.get('permutations')[:1]
        else:
            permutations = self.get('permutations')

        for iperm, perm in enumerate(permutations):
            for iflow, color_flow in enumerate(self.get('color_flows')):
                # 1,2,3,4,5 -> 1,2,4,5,3 e.g.
                perm_replace_dict = \
                          dict(zip(self.get('permutations')[0], perm))
                col_str = color_flow.get('color_string').create_copy()
                col_str.replace_indices(perm_replace_dict)
                # Create the immutable string
                immutable_col_str = col_str.to_immutable()
                if iperm == 0:
                    basic_immutable_factors.append((immutable_col_str,
                                                    col_str.Nc_power))

                # Create the color matrix element between this entry
                # and all (previous) entries in the basic flows
                for ibasic, (basic_string, basic_Nc) in \
                                             enumerate(basic_immutable_factors):
                    # Get matrix entry
                    result, result_fixed_Nc = \
                           col_matrix.get_matrix_entry(canonical_dict,
                                                       immutable_col_str,
                                                       basic_string,
                                                       Nc_power_min = \
                                                       min_Nc_power - \
                                                       basic_Nc - \
                                                       col_str.Nc_power)
                    # Ignore null entries
                    if not result:
                        continue
                    # Set color matrix entry, if non-zero
                    irow = ibasic
                    icol = iperm * len(self.get('color_flows')) + iflow

                    col_matrix[(icol, irow)] = copy.copy(result)

                    # Account for combined Nc colors for results
                    for col in col_matrix[(icol, irow)]:
                        col.Nc_power += basic_Nc + col_str.Nc_power
                        
                    col_matrix.col_matrix_fixed_Nc[(icol, irow)] = \
                                      col_matrix[(icol, irow)].set_Nc(3)
                    if iperm == 0:
                        # Need also inverse entry
                        col_matrix[(irow, icol)] = col_matrix[(icol, irow)]
                        col_matrix.col_matrix_fixed_Nc[(irow, icol)] = \
                               col_matrix.col_matrix_fixed_Nc[(icol, irow)]
                
    def get_external_wavefunctions(self):
        """Redefine HelasMatrixElement.get_external_wavefunctions"""

        return self.get('color_flows')[0].get_external_wavefunctions()

    def get_used_lorentz(self):
        """Return a list of (lorentz_name, conjugate, outgoing) with
        all lorentz structures used by this HelasMultiProcess."""

        helas_list = []

        for me in self.get('color_flows'):
            helas_list.extend(me.get_used_lorentz())

        return list(set(helas_list))

    def get_nexternal_ninitial(self):
        """Gives (number or external particles, number of
        incoming particles)"""

        return self.get('color_flows')[0].get_nexternal_ninitial()

    def get_denominator_factor(self):
        """Calculate the denominator factor due to:
        Averaging initial state color and spin, and
        identical final state particles"""
        
        return self.get('color_flows')[0].get_denominator_factor()
        
#===============================================================================
# COHelasMultiProcess
#===============================================================================
class COHelasMultiProcess(helas_objects.HelasMultiProcess):
    """Version of HelasMultiProcess which generates COHelasMatrixElements"""

    matrix_element_class = COHelasMatrixElement

    def __init__(self, argument=None, gen_color = 3, optimization = 1):
        """Allow initialization with AmplitudeList"""

        if isinstance(argument, ColorOrderedMultiProcess):
            super(COHelasMultiProcess, self).__init__()
            self.set('matrix_elements',
                     self.generate_matrix_elements(argument.get('amplitudes'),
                                                   gen_color, optimization))
        if isinstance(argument, diagram_generation.AmplitudeList):
            super(COHelasMultiProcess, self).__init__()
            self.set('matrix_elements',
                     self.generate_matrix_elements(argument, gen_color,
                                                   optimization))
        elif argument:
            # call the mother routine
            super(COHelasMultiProcess, self).__init__(argument)
        else:
            # call the mother routine
            super(COHelasMultiProcess, self).__init__()

    #===========================================================================
    # generate_matrix_elements
    #===========================================================================
    @classmethod
    def generate_matrix_elements(cls, amplitudes, gen_color = 3,
                                 optimization = 1, decay_ids = []):
        """Generate the HelasMatrixElements for the amplitudes,
        identifying processes with identical matrix elements, as
        defined by HelasMatrixElement.__eq__. Returns a
        HelasMatrixElementList and an amplitude map (used by the
        SubprocessGroup functionality). decay_ids is a list of decayed
        particle ids, since those should not be combined even if
        matrix element is identical."""

        assert isinstance(amplitudes, diagram_generation.AmplitudeList), \
                  "%s is not valid AmplitudeList" % repr(amplitudes)

        # Keep track of already generated color objects, to reuse as
        # much as possible
        list_colorize = []
        list_color_basis = []
        list_color_matrices = []

        matrix_elements = helas_objects.HelasMatrixElementList()

        while amplitudes:
            # Pop the amplitude to save memory space
            amplitude = amplitudes.pop(0)
            if isinstance(amplitude, diagram_generation.DecayChainAmplitude):
                matrix_element_list = helas_objects.HelasDecayChainProcess(amplitude).\
                                      combine_decay_chain_processes()
            else:
                logger.info("Generating Helas calls for %s" % \
                         amplitude.get('process').nice_string().\
                                           replace('Process', 'process'))
                matrix_element_list = [cls.matrix_element_class(amplitude,
                                                   decay_ids=decay_ids,
                                                   gen_color=gen_color,
                                                   optimization = optimization)]
            for matrix_element in matrix_element_list:
                assert isinstance(matrix_element, helas_objects.HelasMatrixElement), \
                          "Not a HelasMatrixElement: %s" % matrix_element

                try:
                    # If an identical matrix element is already in the list,
                    # then simply add this process to the list of
                    # processes for that matrix element
                    me_index = matrix_elements.index(matrix_element)
                    other_processes = matrix_elements[me_index].get('processes')
                    logger.info("Combining process with %s" % \
                      other_processes[0].nice_string().replace('Process: ', ''))
                    other_processes.extend(matrix_element.get('processes'))
                except ValueError:
                    # Otherwise, if the matrix element has any diagrams,
                    # add this matrix element.
                    if matrix_element.get('processes') and \
                           matrix_element.get('color_flows'):
                        matrix_elements.append(matrix_element)

        return matrix_elements
