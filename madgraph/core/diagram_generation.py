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
"""Classes for diagram generation. Amplitude performs the diagram
generation, DecayChainAmplitude keeps track of processes with decay
chains, and MultiProcess allows generation of processes with
multiparticle definitions. DiagramTag allows to identify diagrams
based on relevant properties.
"""

import array
import copy
import itertools
import logging

import madgraph.core.base_objects as base_objects

from madgraph import MadGraph5Error, InvalidCmd
logger = logging.getLogger('madgraph.diagram_generation')

#===============================================================================
# DiagramTag mother class
#===============================================================================

class DiagramTag(object):
    """Class to tag diagrams based on objects with some __lt__ measure, e.g.
    PDG code/interaction id (for comparing diagrams from the same amplitude),
    or Lorentz/coupling/mass/width (for comparing AMPs from different MEs).
    Algorithm: Create chains starting from external particles:
    1 \        / 6
    2 /\______/\ 7
    3_ /  |   \_ 8
    4 /   5    \_ 9
                \ 10
    gives ((((9,10,id910),8,id9108),(6,7,id67),id910867)
           (((1,2,id12),(3,4,id34)),id1234),
           5,id91086712345)
    where idN is the id of the corresponding interaction. The ordering within
    chains is based on chain length (depth; here, 1234 has depth 2, 910867 has
    depth 3, 5 has depht 0), and if equal on the ordering of the chain elements.
    The determination of central vertex is based on minimizing the chain length
    for the longest subchain. 
    This gives a unique tag which can be used to identify diagrams
    (instead of symmetry), as well as identify identical matrix elements from
    different processes."""

    class DiagramTagError(Exception):
        """Exception for any problems in DiagramTags"""
        pass

    def __init__(self, diagram, model = None, ninitial = 2):
        """Initialize with a diagram. Create DiagramTagChainLinks according to
        the diagram, and figure out if we need to shift the central vertex."""

        # wf_dict keeps track of the intermediate particles
        leg_dict = {}
        # Create the chain which will be the diagram tag
        for vertex in diagram.get('vertices'):
            # Only add incoming legs
            legs = vertex.get('legs')[:-1]
            lastvx = vertex == diagram.get('vertices')[-1]
            if lastvx:
                # If last vertex, all legs are incoming
                legs = vertex.get('legs')
            # Add links corresponding to the relevant legs
            link = DiagramTagChainLink([leg_dict.setdefault(leg.get('number'),
                          DiagramTagChainLink(self.link_from_leg(leg, model))) \
                                        for leg in legs],
                                        self.vertex_id_from_vertex(vertex,
                                                                   lastvx,
                                                                   model,
                                                                   ninitial))
            # Add vertex to leg_dict if not last one
            if not lastvx:
                leg_dict[vertex.get('legs')[-1].get('number')] = link

        # The resulting link is the hypothetical result
        self.tag = link

        # Now make sure to find the central vertex in the diagram,
        # defined by the longest leg being as short as possible
        done = max([l.depth for l in self.tag.links]) == 0
        while not done:
            # Identify the longest chain in the tag
            longest_chain = self.tag.links[0]
            # Create a new link corresponding to moving one step
            new_link = DiagramTagChainLink(self.tag.links[1:],
                                           self.flip_vertex(\
                                               self.tag.vertex_id,
                                               longest_chain.vertex_id))
            # Create a new final vertex in the direction of the longest link
            other_link = DiagramTagChainLink(list(longest_chain.links) + \
                                             [new_link],
                                             self.flip_vertex(\
                                                 longest_chain.vertex_id,
                                                 self.tag.vertex_id))
            
            if other_link.links[0] < self.tag.links[0]:
                # Switch to new tag, continue search
                self.tag = other_link
            else:
                # We have found the central vertex
                done = True

    def get_external_numbers(self):
        """Get the order of external particles in this tag"""

        return self.tag.get_external_numbers()

    @staticmethod
    def reorder_permutation(perm, start_perm):
        """Reorder a permutation with respect to start_perm. Note that
        both need to start from 1."""
        if perm == start_perm:
            return range(len(perm))
        order = [i for (p,i) in \
                 sorted([(p,i) for (i,p) in enumerate(perm)])]
        return [start_perm[i]-1 for i in order]

    @staticmethod
    def link_from_leg(leg, model):
        """Returns the default end link for a leg: ((id, state), number).
        Note that the number is not taken into account if tag comparison,
        but is used only to extract leg permutations."""
        if leg.get('state'):
            # Identify identical final state particles
            return [((leg.get('id'), 0), leg.get('number'))]
        else:
            # Distinguish identical initial state particles
            return [((leg.get('id'), leg.get('number')), leg.get('number'))]

    @staticmethod
    def vertex_id_from_vertex(vertex, last_vertex, model, ninitial):
        """Returns the default vertex id: just the interaction id"""
        return vertex.get('id')

    @staticmethod
    def flip_vertex(new_vertex, old_vertex):
        """Returns the default vertex flip: just the new_vertex"""
        return new_vertex

    def __eq__(self, other):
        """Equal if same tag"""
        if type(self) != type(other):
            return False
        return self.tag == other.tag

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return str(self.tag)

    def __lt__(self, other):
        return self.tag < other.tag

    def __gt__(self, other):
        return self.tag > other.tag

    __repr__ = __str__

class DiagramTagChainLink(object):
    """Chain link for a DiagramTag. A link is a tuple + vertex id + depth,
    with a comparison operator defined"""

    def __init__(self, objects, vertex_id = None):
        """Initialize, either with a tuple of DiagramTagChainLinks and
        a vertex_id (defined by DiagramTag.vertex_id_from_vertex), or
        with an external leg object (end link) defined by
        DiagramTag.link_from_leg"""

        if vertex_id == None:
            # This is an end link, corresponding to an external leg
            self.links = tuple(objects)
            self.vertex_id = 0
            self.depth = 0
            self.end_link = True
            return
        # This is an internal link, corresponding to an internal line
        self.links = tuple(sorted(list(tuple(objects)), reverse=True))
        self.vertex_id = vertex_id
        # depth = sum(depth for links) + max(1, len(self.links)-1)
        # in order to get depth 2 for a 4-particle vertex
        self.depth = sum([l.depth for l in self.links],
                         max(1, len(self.links)-1))
        self.end_link = False

    def get_external_numbers(self):
        """Get the permutation of external numbers (assumed to be the
        second entry in the end link tuples)"""

        if self.end_link:
            return [self.links[0][1]]

        return sum([l.get_external_numbers() for l in self.links], [])

    def __lt__(self, other):
        """Compare self with other in the order:
        1. depth 2. len(links) 3. vertex id 4. measure of links"""

        if self == other:
            return False

        if self.depth != other.depth:
            return self.depth < other.depth

        if len(self.links) != len(other.links):
            return len(self.links) < len(other.links)

        if self.vertex_id != other.vertex_id:
            return self.vertex_id < other.vertex_id

        for i, link in enumerate(self.links):
            if i > len(other.links) - 1:
                return False
            if link != other.links[i]:
                return link < other.links[i]

    def __gt__(self, other):
        return self != other and not self.__lt__(other)

    def __eq__(self, other):
        """For end link,
        consider equal if self.links[0][0] == other.links[0][0],
        i.e., ignore the leg number (in links[0][1])."""

        if self.end_link and other.end_link and  \
               self.depth == other.depth and self.vertex_id == other.vertex_id:
            return self.links[0][0] == other.links[0][0]
        
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)


    def __str__(self):
        if self.end_link:
            return str(self.links)
        return "%s, %s; %d" % (str(self.links),
                               str(self.vertex_id),
                               self.depth)

    __repr__ = __str__

#===============================================================================
# Amplitude
#===============================================================================
class Amplitude(base_objects.PhysicsObject):
    """Amplitude: process + list of diagrams (ordered)
    Initialize with a process, then call generate_diagrams() to
    generate the diagrams for the amplitude
    """

    def default_setup(self):
        """Default values for all properties"""

        self['process'] = base_objects.Process()
        self['diagrams'] = None
        # has_mirror_process is True if the same process but with the
        # two incoming particles interchanged has been generated
        self['has_mirror_process'] = False

    def __init__(self, argument=None):
        """Allow initialization with Process"""

        if isinstance(argument, base_objects.Process):
            super(Amplitude, self).__init__()
            self.set('process', argument)
            self.generate_diagrams()
        elif argument != None:
            # call the mother routine
            super(Amplitude, self).__init__(argument)
        else:
            # call the mother routine
            super(Amplitude, self).__init__()

    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name == 'process':
            if not isinstance(value, base_objects.Process):
                raise self.PhysicsObjectError, \
                        "%s is not a valid Process object" % str(value)
        if name == 'diagrams':
            if not isinstance(value, base_objects.DiagramList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid DiagramList object" % str(value)
        if name == 'has_mirror_process':
            if not isinstance(value, bool):
                raise self.PhysicsObjectError, \
                        "%s is not a valid boolean" % str(value)
        return True

    def get(self, name):
        """Get the value of the property name."""

        if name == 'diagrams' and self[name] == None:
            # Have not yet generated diagrams for this process
            if self['process']:
                self.generate_diagrams()

        return Amplitude.__bases__[0].get(self, name)  #return the mother routine


    def get_sorted_keys(self):
        """Return diagram property names as a nicely sorted list."""

        return ['process', 'diagrams', 'has_mirror_process']

    def get_number_of_diagrams(self):
        """Returns number of diagrams for this amplitude"""
        return len(self.get('diagrams'))

    def get_amplitudes(self):
        """Return an AmplitudeList with just this amplitude.
        Needed for DecayChainAmplitude."""

        return AmplitudeList([self])

    def nice_string(self, indent=0):
        """Returns a nicely formatted string of the amplitude content."""
        return self.get('process').nice_string(indent) + "\n" + \
               self.get('diagrams').nice_string(indent)

    def nice_string_processes(self, indent=0):
        """Returns a nicely formatted string of the amplitude process."""
        return self.get('process').nice_string(indent)

    def generate_diagrams(self):
        """Generate diagrams. Algorithm:

        1. Define interaction dictionaries:
          * 2->0 (identity), 3->0, 4->0, ... , maxlegs->0
          * 2 -> 1, 3 -> 1, ..., maxlegs-1 -> 1 

        2. Set flag from_group=true for all external particles.
           Flip particle/anti particle for incoming particles.

        3. If there is a dictionary n->0 with n=number of external
           particles, create if possible the combination [(1,2,3,4,...)] 
           with *at least two* from_group==true. This will give a
           finished (set of) diagram(s) (done by reduce_leglist)

        4. Create all allowed groupings of particles with at least one
           from_group==true (according to dictionaries n->1):
           [(1,2),3,4...],[1,(2,3),4,...],...,
                          [(1,2),(3,4),...],...,[(1,2,3),4,...],... 
           (done by combine_legs)

        5. Replace each group with a (list of) new particle(s) with number 
           n = min(group numbers). Set from_group true for these
           particles and false for all other particles. Store vertex info.
           (done by merge_comb_legs)

        6. Stop algorithm when at most 2 particles remain.
           Return all diagrams (lists of vertices).

        7. Repeat from 3 (recursion done by reduce_leglist)

        8. Replace final p=p vertex
        
        Be aware that the resulting vertices have all particles outgoing,
        so need to flip for incoming particles when used.

        SPECIAL CASE: For A>BC... processes which are legs in decay
        chains, we need to ensure that BC... combine first, giving A=A
        as a final vertex. This case is defined by the Process
        property is_decay_chain = True.
        """

        process = self.get('process')
        model = process.get('model')
        legs = process.get('legs')

        # Make sure orders is the minimum of orders and overall_orders
        for key in process.get('overall_orders').keys():
            try:
                process.get('orders')[key] = \
                                 min(process.get('orders')[key],
                                     process.get('overall_orders')[key])
            except KeyError:
                process.get('orders')[key] = process.get('overall_orders')[key]

        assert model.get('particles'), \
           "particles are missing in model: %s" %  model.get('particles')

        assert model.get('interactions'), \
               "interactions are missing in model" 
                  

        res = base_objects.DiagramList()
        # First check that the number of fermions is even
        if len(filter(lambda leg: model.get('particle_dict')[\
                        leg.get('id')].is_fermion(), legs)) % 2 == 1:
            self['diagrams'] = res
            raise InvalidCmd, 'The number of fermion is odd' 
            return res

        # Then check same number of incoming and outgoing fermions (if
        # no Majorana particles in model)
        if not model.get('got_majoranas') and \
           len(filter(lambda leg: leg.is_incoming_fermion(model), legs)) != \
           len(filter(lambda leg: leg.is_outgoing_fermion(model), legs)):
            self['diagrams'] = res
            raise InvalidCmd, 'The number of of incoming/outcoming fermions are different' 
            return res
        
        # Finally check that charge (conserve by all interactions) of the process
        #is globally conserve for this process.
        for charge in model.get('conserved_charge'):
            total = 0
            for leg in legs:
                part = model.get('particle_dict')[leg.get('id')]
                try:
                    value = part.get(charge)
                except AttributeError, PhysicsObjectError:
                    value = 0
                    
                if (leg.get('id') != part['pdg_code']) != leg['state']:
                    total -= value
                else:
                    total += value

            if abs(total) > 1e-10:
                self['diagrams'] = res
                raise InvalidCmd, 'No %s conservation for this process ' % charge
                return res
                    
        logger.info("Trying %s " % process.nice_string().replace('Process', 'process'))

        # Give numbers to legs in process
        for i in range(0, len(process.get('legs'))):
            # Make sure legs are unique
            leg = copy.copy(process.get('legs')[i])
            process.get('legs')[i] = leg
            if leg.get('number') == 0:
                leg.set('number', i + 1)

        # Copy leglist from process, so we can flip leg identities
        # without affecting the original process
        leglist = base_objects.LegList(\
                       [ copy.copy(leg) for leg in process.get('legs') ])

        for leg in leglist:

            # For the first step, ensure the tag from_group 
            # is true for all legs
            leg.set('from_group', True)

            # Need to flip part-antipart for incoming particles, 
            # so they are all outgoing
            if leg.get('state') == False:
                part = model.get('particle_dict')[leg.get('id')]
                leg.set('id', part.get_anti_pdg_code())

        # Calculate the maximal multiplicity of n-1>1 configurations
        # to restrict possible leg combinations
        max_multi_to1 = max([len(key) for key in \
                             model.get('ref_dict_to1').keys()])


        # Reduce the leg list and return the corresponding
        # list of vertices

        # For decay processes, generate starting from final-state
        # particles and make sure the initial state particle is
        # combined only as the last particle. This allows to use these
        # in decay chains later on.
        is_decay_proc = process.get_ninitial() == 1
        if is_decay_proc:
            part = model.get('particle_dict')[leglist[0].get('id')]
            # For decay chain legs, we want everything to combine to
            # the initial leg. This is done by only allowing the
            # initial leg to combine as a final identity.
            ref_dict_to0 = {(part.get_pdg_code(),part.get_anti_pdg_code()):[0],
                            (part.get_anti_pdg_code(),part.get_pdg_code()):[0]}
            # Need to set initial leg from_group to None, to make sure
            # it can only be combined at the end.
            leglist[0].set('from_group', None)
            reduced_leglist = self.reduce_leglist(leglist,
                                                  max_multi_to1,
                                                  ref_dict_to0,
                                                  is_decay_proc,
                                                  process.get('orders'))
        else:
            reduced_leglist = self.reduce_leglist(leglist,
                                                  max_multi_to1,
                                                  model.get('ref_dict_to0'),
                                                  is_decay_proc,
                                                  process.get('orders'))

        if reduced_leglist:
            for vertex_list in reduced_leglist:
                res.append(base_objects.Diagram(
                            {'vertices':base_objects.VertexList(vertex_list)}))

        # Record whether or not we failed generation before required
        # s-channel propagators are taken into account
        failed_crossing = not res

        # Required s-channels is a list of id-lists. Select the
        # diagrams where all required s-channel propagators in any of
        # the lists are present (i.e., the different lists correspond
        # to "or", while the elements of the list correspond to
        # "and").
        if process.get('required_s_channels') and \
               process.get('required_s_channels')[0]:
            # We shouldn't look at the last vertex in each diagram,
            # since that is the n->0 vertex
            lastvx = -1
            # For decay chain processes, there is an "artificial"
            # extra vertex corresponding to particle 1=1, so we need
            # to exclude the two last vertexes.
            if is_decay_proc: lastvx = -2
            ninitial = len(filter(lambda leg: leg.get('state') == False,
                                  process.get('legs')))
            # Check required s-channels for each list in required_s_channels
            old_res = res
            res = base_objects.DiagramList()
            for id_list in process.get('required_s_channels'):
                res_diags = filter(lambda diagram: \
                          all([req_s_channel in \
                               [vertex.get_s_channel_id(\
                               process.get('model'), ninitial) \
                               for vertex in diagram.get('vertices')[:lastvx]] \
                               for req_s_channel in \
                               id_list]), old_res)
                # Add diagrams only if not already in res
                res.extend([diag for diag in res_diags if diag not in res])

        # Select the diagrams where no forbidden s-channel propagators
        # are present.
        # Note that we shouldn't look at the last vertex in each
        # diagram, since that is the n->0 vertex
        if process.get('forbidden_s_channels'):
            ninitial = len(filter(lambda leg: leg.get('state') == False,
                              process.get('legs')))
            verts = base_objects.VertexList(sum([[vertex for vertex \
                                                  in diagram.get('vertices')[:-1]
                                           if vertex.get_s_channel_id(\
                                               process.get('model'), ninitial) \
                                           in process.get('forbidden_s_channels')] \
                                               for diagram in res], []))
            for vert in verts:
                # Use onshell = False to indicate that this s-channel is forbidden
                newleg = copy.copy(vert.get('legs').pop(-1))
                newleg.set('onshell', False)
                vert.get('legs').append(newleg)
                
        # Set diagrams to res
        self['diagrams'] = res

        # Set actual coupling orders for each diagram
        for diagram in self['diagrams']:
            diagram.calculate_orders(model)

        # Set the coupling orders of the process
        if self['diagrams']:
            coupling_orders = {}        
            for key in sorted(list(model.get_coupling_orders())) + ['WEIGHTED']:
                coupling_orders[key] = max([d.get('orders')[key] for \
                                            d in self.get('diagrams')])
            process.set('orders', coupling_orders)

        # Replace final id=0 vertex if necessary
        if not process.get('is_decay_chain'):
            for diagram in self['diagrams']:
                vertices = diagram.get('vertices')
                if len(vertices) > 1 and vertices[-1].get('id') == 0:
                    # Need to "glue together" last and next-to-last
                    # vertex, by replacing the (incoming) last leg of the
                    # next-to-last vertex with the (outgoing) leg in the
                    # last vertex
                    vertices = copy.copy(vertices)
                    lastvx = vertices.pop()
                    nexttolastvertex = copy.copy(vertices.pop())
                    legs = copy.copy(nexttolastvertex.get('legs'))
                    ntlnumber = legs[-1].get('number')
                    lastleg = filter(lambda leg: leg.get('number') != ntlnumber,
                                     lastvx.get('legs'))[0]
                    # Reset onshell in case we have forbidden s-channels
                    if lastleg.get('onshell') == False:
                        lastleg.set('onshell', None)
                    # Replace the last leg of nexttolastvertex
                    legs[-1] = lastleg
                    nexttolastvertex.set('legs', legs)
                    vertices.append(nexttolastvertex)
                    diagram.set('vertices', vertices)

        if res:
            logger.info("Process has %d diagrams" % len(res))

        # Trim down number of legs and vertices used to save memory
        self.trim_diagrams()

        # Sort process legs according to leg number
        self.get('process').get('legs').sort()

        return not failed_crossing

    def reduce_leglist(self, curr_leglist, max_multi_to1, ref_dict_to0,
                       is_decay_proc = False, coupling_orders = None):
        """Recursive function to reduce N LegList to N-1
           For algorithm, see doc for generate_diagrams.
        """

        # Result variable which is a list of lists of vertices
        # to be added
        res = []

        # Stop condition. If LegList is None, that means that this
        # diagram must be discarded
        if curr_leglist is None:
            return None

        # Extract ref dict information
        model = self.get('process').get('model')
        ref_dict_to1 = self.get('process').get('model').get('ref_dict_to1')


        # If all legs can be combined in one single vertex, add this
        # vertex to res and continue.
        # Special treatment for decay chain legs

        if curr_leglist.can_combine_to_0(ref_dict_to0, is_decay_proc):
            # Extract the interaction id associated to the vertex 
            vertex_ids = self.get_combined_vertices(curr_leglist,
                       copy.copy(ref_dict_to0[tuple(sorted([leg.get('id') for \
                                                       leg in curr_leglist]))]))

            final_vertices = [base_objects.Vertex({'legs':curr_leglist,
                                                   'id':vertex_id}) for \
                              vertex_id in vertex_ids]
            # Check for coupling orders. If orders < 0, skip vertex
            for final_vertex in final_vertices:
                if self.reduce_orders(coupling_orders, model,
                                      [final_vertex.get('id')]) != False:
                    res.append([final_vertex])
        # Stop condition 2: if the leglist contained exactly two particles,
        # return the result, if any, and stop.
        if len(curr_leglist) == 2:
            if res:
                return res
            else:
                return None

        # Create a list of all valid combinations of legs
        comb_lists = self.combine_legs(curr_leglist,
                                       ref_dict_to1, max_multi_to1)

        # Create a list of leglists/vertices by merging combinations
        leg_vertex_list = self.merge_comb_legs(comb_lists, ref_dict_to1)

        # Consider all the pairs
        for leg_vertex_tuple in leg_vertex_list:

            # Remove forbidden particles
            if self.get('process').get('forbidden_particles') and \
                any([abs(vertex.get('legs')[-1].get('id')) in \
                self.get('process').get('forbidden_particles') \
                for vertex in leg_vertex_tuple[1]]):
                    continue

            # Check for coupling orders. If couplings < 0, skip recursion.
            new_coupling_orders = self.reduce_orders(coupling_orders,
                                                     model,
                                                     [vertex.get('id') for vertex in \
                                                      leg_vertex_tuple[1]])
            if new_coupling_orders == False:
                # Some coupling order < 0
                continue

            # This is where recursion happens
            # First, reduce again the leg part
            reduced_diagram = self.reduce_leglist(leg_vertex_tuple[0],
                                                  max_multi_to1,
                                                  ref_dict_to0,
                                                  is_decay_proc,
                                                  new_coupling_orders)
            # If there is a reduced diagram
            if reduced_diagram:
                vertex_list_list = [list(leg_vertex_tuple[1])]
                vertex_list_list.append(reduced_diagram)
                expanded_list = expand_list_list(vertex_list_list)
                res.extend(expanded_list)

        return res

    def reduce_orders(self, coupling_orders, model, vertex_id_list):
        """Return False if the coupling orders for any coupling is <
        0, otherwise return the new coupling orders with the vertex
        orders subtracted. If coupling_orders is not given, return
        None (which counts as success).
        WEIGHTED is a special order, which corresponds to the sum of
        order hierarchys for the couplings."""

        if not coupling_orders:
            return None

        present_couplings = copy.copy(coupling_orders)
        for id in vertex_id_list:
            # Don't check for identity vertex (id = 0)
            if not id:
                continue
            inter = model.get("interaction_dict")[id]
            for coupling in inter.get('orders').keys():
                # Note that we don't consider a missing coupling as a
                # constraint
                if coupling in present_couplings:
                    # Reduce the number of couplings that are left
                    present_couplings[coupling] -= \
                             inter.get('orders')[coupling]
                    if present_couplings[coupling] < 0:
                        # We have too many couplings of this type
                        return False
            # Now check for WEIGHTED, i.e. the sum of coupling hierarchy values
            if 'WEIGHTED' in present_couplings:
                weight = sum([model.get('order_hierarchy')[c]*n for \
                              (c,n) in inter.get('orders').items()])
                present_couplings['WEIGHTED'] -= weight
                if present_couplings['WEIGHTED'] < 0:
                        # Total coupling weight too large
                        return False
                
        return present_couplings

    def combine_legs(self, list_legs, ref_dict_to1, max_multi_to1):
        """Recursive function. Take a list of legs as an input, with
        the reference dictionary n-1->1, and output a list of list of
        tuples of Legs (allowed combinations) and Legs (rest). Algorithm:

        1. Get all n-combinations from list [123456]: [12],..,[23],..,[123],..

        2. For each combination, say [34]. Check if combination is valid.
           If so:

           a. Append [12[34]56] to result array

           b. Split [123456] at index(first element in combination+1),
              i.e. [12],[456] and subtract combination from second half,
              i.e.: [456]-[34]=[56]. Repeat from 1. with this array

        3. Take result array from call to 1. (here, [[56]]) and append
           (first half in step b - combination) + combination + (result
           from 1.) = [12[34][56]] to result array

        4. After appending results from all n-combinations, return
           resulting array. Example, if [13] and [45] are valid
           combinations:
            [[[13]2456],[[13]2[45]6],[123[45]6]] 
        """

        res = []

        # loop over possible combination lengths (+1 is for range convention!)
        for comb_length in range(2, max_multi_to1 + 1):

            # Check the considered length is not longer than the list length
            if comb_length > len(list_legs):
                return res

            # itertools.combinations returns all possible combinations
            # of comb_length elements from list_legs
            for comb in itertools.combinations(list_legs, comb_length):

                # Check if the combination is valid
                if base_objects.LegList(comb).can_combine_to_1(ref_dict_to1):

                    # Identify the rest, create a list [comb,rest] and
                    # add it to res
                    res_list = copy.copy(list_legs)
                    for leg in comb:
                        res_list.remove(leg)
                    res_list.insert(list_legs.index(comb[0]), comb)
                    res.append(res_list)

                    # Now, deal with cases with more than 1 combination

                    # First, split the list into two, according to the
                    # position of the first element in comb, and remove
                    # all elements form comb
                    res_list1 = list_legs[0:list_legs.index(comb[0])]
                    res_list2 = list_legs[list_legs.index(comb[0]) + 1:]
                    for leg in comb[1:]:
                        res_list2.remove(leg)

                    # Create a list of type [comb,rest1,rest2(combined)]
                    res_list = res_list1
                    res_list.append(comb)
                    # This is where recursion actually happens, 
                    # on the second part
                    for item in self.combine_legs(res_list2,
                                                  ref_dict_to1,
                                                  max_multi_to1):
                        final_res_list = copy.copy(res_list)
                        final_res_list.extend(item)
                        res.append(final_res_list)

        return res


    def merge_comb_legs(self, comb_lists, ref_dict_to1):
        """Takes a list of allowed leg combinations as an input and returns
        a set of lists where combinations have been properly replaced
        (one list per element in the ref_dict, so that all possible intermediate
        particles are included). For each list, give the list of vertices
        corresponding to the executed merging, group the two as a tuple.
        """

        res = []

        for comb_list in comb_lists:

            reduced_list = []
            vertex_list = []

            for entry in comb_list:

                # Act on all leg combinations
                if isinstance(entry, tuple):

                    # Build the leg object which will replace the combination:
                    # 1) leg ids is as given in the ref_dict
                    leg_vert_ids = copy.copy(ref_dict_to1[\
                        tuple(sorted([leg.get('id') for leg in entry]))])
                    # 2) number is the minimum of leg numbers involved in the
                    # combination
                    number = min([leg.get('number') for leg in entry])
                    # 3) state is final, unless there is exactly one initial 
                    # state particle involved in the combination -> t-channel
                    if len(filter(lambda leg: leg.get('state') == False,
                                  entry)) == 1:
                        state = False
                    else:
                        state = True
                    # 4) from_group is True, by definition

                    # Create and add the object. This is done by a
                    # separate routine, to allow overloading by
                    # daughter classes
                    new_leg_vert_ids = []
                    if leg_vert_ids:
                        new_leg_vert_ids = self.get_combined_legs(entry,
                                                                  leg_vert_ids,
                                                                  number,
                                                                  state)
                    
                    reduced_list.append([l[0] for l in new_leg_vert_ids])


                    # Create and add the corresponding vertex
                    # Extract vertex ids corresponding to the various legs
                    # in mylegs
                    vlist = base_objects.VertexList()
                    for (myleg, vert_id) in new_leg_vert_ids:
                        # Start with the considered combination...
                        myleglist = base_objects.LegList(list(entry))
                        # ... and complete with legs after reducing
                        myleglist.append(myleg)
                        # ... and consider the correct vertex id
                        vlist.append(base_objects.Vertex(
                                         {'legs':myleglist,
                                          'id':vert_id}))

                    vertex_list.append(vlist)

                # If entry is not a combination, switch the from_group flag
                # and add it
                else:
                    cp_entry = copy.copy(entry)
                    # Need special case for from_group == None; this
                    # is for initial state leg of decay chain process
                    # (see Leg.can_combine_to_0)
                    if cp_entry.get('from_group') != None:
                        cp_entry.set('from_group', False)
                    reduced_list.append(cp_entry)

            # Flatten the obtained leg and vertex lists
            flat_red_lists = expand_list(reduced_list)
            flat_vx_lists = expand_list(vertex_list)

            # Combine the two lists in a list of tuple
            for i in range(0, len(flat_vx_lists)):
                res.append((base_objects.LegList(flat_red_lists[i]), \
                            base_objects.VertexList(flat_vx_lists[i])))

        return res

    def get_combined_legs(self, legs, leg_vert_ids, number, state):
        """Create a set of new legs from the info given. This can be
        overloaded by daughter classes."""

        mylegs = [(base_objects.Leg({'id':leg_id,
                                    'number':number,
                                    'state':state,
                                    'from_group':True}),
                   vert_id)\
                  for leg_id, vert_id in leg_vert_ids]

        return mylegs
                          
    def get_combined_vertices(self, legs, vert_ids):
        """Allow for selection of vertex ids. This can be
        overloaded by daughter classes."""

        return vert_ids
                          
    def trim_diagrams(self, decay_ids = []):
        """Reduce the number of legs and vertices used in memory."""

        legs = []
        vertices = []

        for diagram in self.get('diagrams'):
            # Keep track of external legs (leg numbers already used)
            leg_external = set()
            for ivx, vertex in enumerate(diagram.get('vertices')):
                for ileg, leg in enumerate(vertex.get('legs')):
                    # Ensure that only external legs get decay flag
                    if leg.get('state') and leg.get('id') in decay_ids and \
                           leg.get('number') not in leg_external:
                        # Use onshell to indicate decaying legs,
                        # i.e. legs that have decay chains
                        leg = copy.copy(leg)
                        leg.set('onshell', True)
                    try:
                        index = legs.index(leg)
                    except ValueError:
                        vertex.get('legs')[ileg] = leg
                        legs.append(leg)
                    else: # Found a leg
                        vertex.get('legs')[ileg] = legs[index]
                    leg_external.add(leg.get('number'))
                try:
                    index = vertices.index(vertex)
                    diagram.get('vertices')[ivx] = vertices[index]
                except ValueError:
                    vertices.append(vertex)
        

#===============================================================================
# AmplitudeList
#===============================================================================
class AmplitudeList(base_objects.PhysicsObjectList):
    """List of Amplitude objects
    """

    def is_valid_element(self, obj):
        """Test if object obj is a valid Amplitude for the list."""

        return isinstance(obj, Amplitude)
    
#===============================================================================
# DecayChainAmplitude
#===============================================================================
class DecayChainAmplitude(Amplitude):
    """A list of amplitudes + a list of decay chain amplitude lists;
    corresponding to a ProcessDefinition with a list of decay chains
    """

    def default_setup(self):
        """Default values for all properties"""

        self['amplitudes'] = AmplitudeList()
        self['decay_chains'] = DecayChainAmplitudeList()

    def __init__(self, argument = None, collect_mirror_procs = False,
                 ignore_six_quark_processes = False):
        """Allow initialization with Process and with ProcessDefinition"""

        if isinstance(argument, base_objects.Process):
            super(DecayChainAmplitude, self).__init__()
            if isinstance(argument, base_objects.ProcessDefinition):
                self['amplitudes'].extend(\
                MultiProcess.generate_multi_amplitudes(argument,
                                                    collect_mirror_procs,
                                                    ignore_six_quark_processes))
            else:
                self['amplitudes'].append(Amplitude(argument))
                # Clean decay chains from process, since we haven't
                # combined processes with decay chains yet
                process = copy.copy(self.get('amplitudes')[0].get('process'))
                process.set('decay_chains', base_objects.ProcessList())
                self['amplitudes'][0].set('process', process)
            for process in argument.get('decay_chains'):
                process.set('overall_orders', argument.get('overall_orders'))
                if not process.get('is_decay_chain'):
                    process.set('is_decay_chain',True)
                if not process.get_ninitial() == 1:
                    raise MadGraph5Error,\
                          "Decay chain process must have exactly one" + \
                          " incoming particle"
                self['decay_chains'].append(\
                    DecayChainAmplitude(process, collect_mirror_procs,
                                        ignore_six_quark_processes))
            # Flag decaying legs in the core process by onshell = True
            decay_ids = sum([[a.get('process').get('legs')[0].get('id') \
                              for a in dec.get('amplitudes')] for dec in \
                             self['decay_chains']], [])
            decay_ids = set(decay_ids)
            for amp in self['amplitudes']:
                amp.trim_diagrams(decay_ids)                             
        elif argument != None:
            # call the mother routine
            super(DecayChainAmplitude, self).__init__(argument)
        else:
            # call the mother routine
            super(DecayChainAmplitude, self).__init__()

    def filter(self, name, value):
        """Filter for valid amplitude property values."""

        if name == 'amplitudes':
            if not isinstance(value, AmplitudeList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid AmplitudeList" % str(value)
        if name == 'decay_chains':
            if not isinstance(value, DecayChainAmplitudeList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid DecayChainAmplitudeList object" % \
                        str(value)
        return True

    def get_sorted_keys(self):
        """Return diagram property names as a nicely sorted list."""

        return ['amplitudes', 'decay_chains']

    # Helper functions

    def get_number_of_diagrams(self):
        """Returns number of diagrams for this amplitude"""
        return sum(len(a.get('diagrams')) for a in self.get('amplitudes')) \
               + sum(d.get_number_of_diagrams() for d in \
                                        self.get('decay_chains'))

    def nice_string(self, indent = 0):
        """Returns a nicely formatted string of the amplitude content."""
        mystr = ""
        for amplitude in self.get('amplitudes'):
            mystr = mystr + amplitude.nice_string(indent) + "\n"

        if self.get('decay_chains'):
            mystr = mystr + " " * indent + "Decays:\n"
        for dec in self.get('decay_chains'):
            mystr = mystr + dec.nice_string(indent + 2) + "\n"

        return  mystr[:-1]

    def nice_string_processes(self, indent = 0):
        """Returns a nicely formatted string of the amplitude processes."""
        mystr = ""
        for amplitude in self.get('amplitudes'):
            mystr = mystr + amplitude.nice_string_processes(indent) + "\n"

        if self.get('decay_chains'):
            mystr = mystr + " " * indent + "Decays:\n"
        for dec in self.get('decay_chains'):
            mystr = mystr + dec.nice_string_processes(indent + 2) + "\n"

        return  mystr[:-1]

    def get_decay_ids(self):
        """Returns a set of all particle ids for which a decay is defined"""

        decay_ids = []

        # Get all amplitudes for the decay processes
        for amp in sum([dc.get('amplitudes') for dc \
                        in self['decay_chains']], []):
            # For each amplitude, find the initial state leg
            decay_ids.append(amp.get('process').get_initial_ids()[0])
            
        # Return a list with unique ids
        return list(set(decay_ids))
    
    def get_amplitudes(self):
        """Recursive function to extract all amplitudes for this process"""

        amplitudes = AmplitudeList()

        amplitudes.extend(self.get('amplitudes'))
        for decay in self.get('decay_chains'):
            amplitudes.extend(decay.get_amplitudes())

        return amplitudes
            

#===============================================================================
# DecayChainAmplitudeList
#===============================================================================
class DecayChainAmplitudeList(base_objects.PhysicsObjectList):
    """List of DecayChainAmplitude objects
    """

    def is_valid_element(self, obj):
        """Test if object obj is a valid DecayChainAmplitude for the list."""

        return isinstance(obj, DecayChainAmplitude)

    
#===============================================================================
# MultiProcess
#===============================================================================
class MultiProcess(base_objects.PhysicsObject):
    """MultiProcess: list of process definitions
                     list of processes (after cleaning)
                     list of amplitudes (after generation)
    """

    def default_setup(self):
        """Default values for all properties"""

        self['process_definitions'] = base_objects.ProcessDefinitionList()
        # self['amplitudes'] can be an AmplitudeList or a
        # DecayChainAmplitudeList, depending on whether there are
        # decay chains in the process definitions or not.
        self['amplitudes'] = AmplitudeList()
        # Flag for whether to combine IS mirror processes together
        self['collect_mirror_procs'] = False
        # List of quark flavors where we ignore processes with at
        # least 6 quarks (three quark lines)
        self['ignore_six_quark_processes'] = []

    def __init__(self, argument=None, collect_mirror_procs = False,
                 ignore_six_quark_processes = []):
        """Allow initialization with ProcessDefinition or
        ProcessDefinitionList"""

        if isinstance(argument, base_objects.ProcessDefinition):
            super(MultiProcess, self).__init__()
            self['process_definitions'].append(argument)
        elif isinstance(argument, base_objects.ProcessDefinitionList):
            super(MultiProcess, self).__init__()
            self['process_definitions'] = argument
        elif argument != None:
            # call the mother routine
            super(MultiProcess, self).__init__(argument)
        else:
            # call the mother routine
            super(MultiProcess, self).__init__()

        self['collect_mirror_procs'] = collect_mirror_procs
        self['ignore_six_quark_processes'] = ignore_six_quark_processes
        
        if isinstance(argument, base_objects.ProcessDefinition) or \
               isinstance(argument, base_objects.ProcessDefinitionList):
            # Generate the diagrams
            self.get('amplitudes')


    def filter(self, name, value):
        """Filter for valid process property values."""

        if name == 'process_definitions':
            if not isinstance(value, base_objects.ProcessDefinitionList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid ProcessDefinitionList object" % str(value)

        if name == 'amplitudes':
            if not isinstance(value, diagram_generation.AmplitudeList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid AmplitudeList object" % str(value)

        if name in ['collect_mirror_procs']:
            if not isinstance(value, bool):
                raise self.PhysicsObjectError, \
                        "%s is not a valid boolean" % str(value)

        if name == 'ignore_six_quark_processes':
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list" % str(value)

        return True

    def get(self, name):
        """Get the value of the property name."""

        if (name == 'amplitudes') and not self[name]:
            for process_def in self.get('process_definitions'):
                if process_def.get('decay_chains'):
                    # This is a decay chain process
                    # Store amplitude(s) as DecayChainAmplitude
                    self['amplitudes'].append(\
                        DecayChainAmplitude(process_def,
                                       self.get('collect_mirror_procs'),
                                       self.get('ignore_six_quark_processes')))
                else:
                    self['amplitudes'].extend(\
                       self.generate_multi_amplitudes(process_def,
                                       self.get('collect_mirror_procs'),
                                       self.get('ignore_six_quark_processes')))

        return MultiProcess.__bases__[0].get(self, name) # call the mother routine

    def get_sorted_keys(self):
        """Return process property names as a nicely sorted list."""

        return ['process_definitions', 'amplitudes']

    amplitude_class = Amplitude

    @classmethod
    def generate_multi_amplitudes(cls, process_definition,
                                  collect_mirror_procs = False,
                                  ignore_six_quark_processes = []):
        """Generate amplitudes in a semi-efficient way.
        Make use of crossing symmetry for processes that fail diagram
        generation, but not for processes that succeed diagram
        generation.  Doing so will risk making it impossible to
        identify processes with identical amplitudes.
        """
        assert isinstance(process_definition, base_objects.ProcessDefinition), \
                                    "%s not valid ProcessDefinition object" % \
                                    repr(process_definition)

        # Set automatic coupling orders
        process_definition.set('orders', MultiProcess.\
                               find_optimal_process_orders(process_definition))
        # Check for maximum orders from the model
        process_definition.check_expansion_orders()

        processes = base_objects.ProcessList()
        amplitudes = AmplitudeList()

        # failed_procs and success_procs are sorted processes that have
        # already failed/succeeded based on crossing symmetry
        failed_procs = []
        success_procs = []
        # Complete processes, for identification of mirror processes
        non_permuted_procs = []
        # permutations keeps the permutations of the crossed processes
        permutations = []

        # Store the diagram tags for processes, to allow for
        # identifying identical matrix elements already at this stage.
        model = process_definition['model']
        
        isids = [leg['ids'] for leg in process_definition['legs'] \
                 if leg['state'] == False]
        fsids = [leg['ids'] for leg in process_definition['legs'] \
                 if leg['state'] == True]

        # Generate all combinations for the initial state
        
        for prod in apply(itertools.product, isids):
            islegs = [\
                    base_objects.Leg({'id':id, 'state': False}) \
                    for id in prod]

            # Generate all combinations for the final state, and make
            # sure to remove double counting

            red_fsidlist = []

            for prod in apply(itertools.product, fsids):

                # Remove double counting between final states
                if tuple(sorted(prod)) in red_fsidlist:
                    continue
                
                red_fsidlist.append(tuple(sorted(prod)));
                
                # Generate leg list for process
                leg_list = [copy.copy(leg) for leg in islegs]
                
                leg_list.extend([\
                        base_objects.Leg({'id':id, 'state': True}) \
                        for id in prod])
                
                legs = base_objects.LegList(leg_list)

                # Check for crossed processes
                sorted_legs = sorted([(l,i+1) for (i,l) in \
                                   enumerate(legs.get_outgoing_id_list(model))])
                permutation = [l[1] for l in sorted_legs]
                sorted_legs = array.array('i', [l[0] for l in sorted_legs])

                # Check for six-quark processes
                if ignore_six_quark_processes and \
                       len([i for i in sorted_legs if abs(i) in \
                            ignore_six_quark_processes]) >= 6:
                    continue
                    
                # Check if crossed process has already failed,
                # in that case don't check process
                if sorted_legs in failed_procs:
                    continue

                # Setup process
                process = base_objects.Process({\
                              'legs':legs,
                              'model':process_definition.get('model'),
                              'id': process_definition.get('id'),
                              'orders': process_definition.get('orders'),
                              'required_s_channels': \
                                 process_definition.get('required_s_channels'),
                              'forbidden_s_channels': \
                                 process_definition.get('forbidden_s_channels'),
                              'forbidden_particles': \
                                 process_definition.get('forbidden_particles'),
                              'is_decay_chain': \
                                 process_definition.get('is_decay_chain'),
                              'overall_orders': \
                                 process_definition.get('overall_orders')})
                
                fast_proc = \
                          array.array('i',[leg.get('id') for leg in legs])
                if collect_mirror_procs:
                    # Check if mirrored process is already generated
                    mirror_proc = \
                              array.array('i', [fast_proc[1], fast_proc[0]] + \
                                          list(fast_proc[2:]))
                    try:
                        mirror_amp = \
                               amplitudes[non_permuted_procs.index(mirror_proc)]
                    except:
                        # Didn't find any mirror process
                        pass
                    else:
                        # Mirror process found
                        mirror_amp.set('has_mirror_process', True)
                        logger.info("Process %s added to mirror process %s" % \
                                    (process.base_string(),
                                     mirror_amp.get('process').base_string()))
                        continue
                        
                # Check for successful crossings, unless we have specified
                # properties that break crossing symmetry
                if not process.get('required_s_channels') and \
                   not process.get('forbidden_s_channels') and \
                   not process.get('is_decay_chain'):
                    try:
                        crossed_index = success_procs.index(sorted_legs)
                    except ValueError:
                        # No crossing found, just continue
                        pass
                    else:
                        # Found crossing - reuse amplitude
                        amplitude = MultiProcess.cross_amplitude(\
                            amplitudes[crossed_index],
                            process,
                            permutations[crossed_index],
                            permutation)
                        amplitudes.append(amplitude)
                        success_procs.append(sorted_legs)
                        permutations.append(permutation)
                        non_permuted_procs.append(fast_proc)
                        logger.info("Crossed process found for %s, reuse diagrams." % \
                                    process.base_string())
                        continue
                    
                # Create new amplitude
                amplitude = cls.amplitude_class({"process": process})
                
                try:
                    result = amplitude.generate_diagrams()
                except InvalidCmd as error:
                    failed_procs.append(sorted_legs)
                else:
                    # Succeeded in generating diagrams
                    if amplitude.get('diagrams'):
                        amplitudes.append(amplitude)
                        success_procs.append(sorted_legs)
                        permutations.append(permutation)
                        non_permuted_procs.append(fast_proc)
                    elif not result:
                        # Diagram generation failed for all crossings
                        failed_procs.append(sorted_legs)
 
        # Raise exception if there are no amplitudes for this process
        if not amplitudes:
            if len(failed_procs) == 1 and 'error' in locals():
                raise error
            else:
                raise MadGraph5Error, \
            "No amplitudes generated from process %s. Please enter a valid process" % \
                  process_definition.nice_string()
        

        # Return the produced amplitudes
        return amplitudes
            
    @classmethod
    def find_optimal_process_orders(cls, process_definition):
        """Find the minimal WEIGHTED order for this set of processes.

        The algorithm:

        1) Check the coupling hierarchy of the model. Assign all
        particles to the different coupling hierarchies so that a
        particle is considered to be in the highest hierarchy (i.e.,
        with lowest value) where it has an interaction.
        
        2) Pick out the legs in the multiprocess according to the
        highest hierarchy represented (so don't mix particles from
        different hierarchy classes in the same multiparticles!)

        3) Find the starting maximum WEIGHTED order as the sum of the
        highest n-2 weighted orders

        4) Pick out required s-channel particle hierarchies, and use
        the highest of the maximum WEIGHTED order from the legs and
        the minimum WEIGHTED order extracted from 2*s-channel
        hierarchys plus the n-2-2*(number of s-channels) lowest
        leg weighted orders.

        5) Run process generation with the WEIGHTED order determined
        in 3)-4) - # final state gluons, with all gluons removed from
        the final state

        6) If no process is found, increase WEIGHTED order by 1 and go
        back to 5), until we find a process which passes. Return that
        order.

        7) Continue 5)-6) until we reach (n-2)*(highest hierarchy)-1.
        If still no process has passed, return
        WEIGHTED = (n-2)*(highest hierarchy)
        """

        assert isinstance(process_definition, base_objects.ProcessDefinition), \
                                    "%s not valid ProcessDefinition object" % \
                                    repr(process_definition)

        processes = base_objects.ProcessList()
        amplitudes = AmplitudeList()

        # If there are already couplings defined, return
        if process_definition.get('orders') or \
               process_definition.get('overall_orders'):
            return process_definition.get('orders')

        logger.info("Checking for minimal orders which gives processes.")
        logger.info("Please specify coupling orders to bypass this step.")

        # Calculate minimum starting guess for WEIGHTED order
        max_order_now, particles, hierarchy = \
                                       process_definition.get_minimum_WEIGHTED()
        coupling = 'WEIGHTED'

        model = process_definition.get('model')
        
        # Extract the initial and final leg ids
        isids = [leg['ids'] for leg in \
                 filter(lambda leg: leg['state'] == False, process_definition['legs'])]
        fsids = [leg['ids'] for leg in \
                 filter(lambda leg: leg['state'] == True, process_definition['legs'])]

        # Run diagram generation with increasing max_order_now until
        # we manage to get diagrams
        while max_order_now < len(fsids)*max(hierarchy):

            logger.info("Trying coupling order WEIGHTED=%d" % max_order_now)

            oldloglevel = logger.getEffectiveLevel()
            logger.setLevel(logging.WARNING)

            # failed_procs are processes that have already failed
            # based on crossing symmetry
            failed_procs = []
            
            # Generate all combinations for the initial state        
            for prod in apply(itertools.product, isids):
                islegs = [ base_objects.Leg({'id':id, 'state': False}) \
                        for id in prod]

                # Generate all combinations for the final state, and make
                # sure to remove double counting

                red_fsidlist = []

                for prod in apply(itertools.product, fsids):

                    # Remove double counting between final states
                    if tuple(sorted(prod)) in red_fsidlist:
                        continue

                    red_fsidlist.append(tuple(sorted(prod)));

                    # Remove gluons from final state if QCD is among
                    # the highest coupling hierarchy
                    nglue = 0
                    if 21 in particles[0]:
                        nglue = len([id for id in prod if id == 21])
                        prod = [id for id in prod if id != 21]

                    # Generate leg list for process
                    leg_list = [copy.copy(leg) for leg in islegs]

                    leg_list.extend([\
                            base_objects.Leg({'id':id, 'state': True}) \
                            for id in prod])

                    legs = base_objects.LegList(leg_list)

                    # Set summed coupling order according to max_order_now
                    # subtracting the removed gluons
                    coupling_orders_now = {coupling: max_order_now - \
                                           nglue}

                    # Setup process
                    process = base_objects.Process({\
                              'legs':legs,
                              'model':model,
                              'id': process_definition.get('id'),
                              'orders': coupling_orders_now,
                              'required_s_channels': \
                                 process_definition.get('required_s_channels'),
                              'forbidden_s_channels': \
                                 process_definition.get('forbidden_s_channels'),
                              'forbidden_particles': \
                                 process_definition.get('forbidden_particles'),
                              'is_decay_chain': \
                                 process_definition.get('is_decay_chain'),
                              'overall_orders': \
                                 process_definition.get('overall_orders')})

                    # Check for couplings with given expansion orders
                    process.check_expansion_orders()

                    # Check for crossed processes
                    sorted_legs = sorted(legs.get_outgoing_id_list(model))
                    # Check if crossed process has already failed
                    # In that case don't check process
                    if tuple(sorted_legs) in failed_procs:
                        continue

                    amplitude = cls.amplitude_class({"process": process})
                    try:
                        amplitude.generate_diagrams()
                    except InvalidCmd:
                        failed_procs.append(tuple(sorted_legs))
                    else:
                        if amplitude.get('diagrams'):
                            # We found a valid amplitude. Return this order number
                            logger.setLevel(oldloglevel)
                            return {coupling: max_order_now}
                        else:
                            failed_procs.append(tuple(sorted_legs))

            # No processes found, increase max_order_now
            max_order_now += 1
            logger.setLevel(oldloglevel)

        # If no valid processes found with nfinal-1 couplings, return maximal
        return {coupling: len(fsids)*max(hierarchy)}

    @staticmethod
    def cross_amplitude(amplitude, process, org_perm, new_perm):
        """Return the amplitude crossed with the permutation new_perm"""
        # Create dict from original leg numbers to new leg numbers
        perm_map = dict(zip(org_perm, new_perm))
        # Initiate new amplitude
        new_amp = copy.copy(amplitude)
        # Number legs
        for i, leg in enumerate(process.get('legs')):
            leg.set('number', i+1)
        # Set process
        new_amp.set('process', process)
        # Now replace the leg numbers in the diagrams
        diagrams = base_objects.DiagramList([d.renumber_legs(perm_map,
                                             process.get('legs'),) for \
                                             d in new_amp.get('diagrams')])
        new_amp.set('diagrams', diagrams)
        new_amp.trim_diagrams()

        # Make sure to reset mirror process
        new_amp.set('has_mirror_process', False)
        
        return new_amp
        
#===============================================================================
# Global helper methods
#===============================================================================

def expand_list(mylist):
    """Takes a list of lists and elements and returns a list of flat lists.
    Example: [[1,2], 3, [4,5]] -> [[1,3,4], [1,3,5], [2,3,4], [2,3,5]]
    """

    # Check that argument is a list
    assert isinstance(mylist, list), "Expand_list argument must be a list"

    res = []

    tmplist = []
    for item in mylist:
        if isinstance(item, list):
            tmplist.append(item)
        else:
            tmplist.append([item])

    for item in apply(itertools.product, tmplist):
        res.append(list(item))

    return res

def expand_list_list(mylist):
    """Recursive function. Takes a list of lists and lists of lists
    and returns a list of flat lists.
    Example: [[1,2],[[4,5],[6,7]]] -> [[1,2,4,5], [1,2,6,7]]
    """

    res = []

    if not mylist or len(mylist) == 1 and not mylist[0]:
        return [[]]

    # Check the first element is at least a list
    assert isinstance(mylist[0], list), \
              "Expand_list_list needs a list of lists and lists of lists"

    # Recursion stop condition, one single element
    if len(mylist) == 1:
        if isinstance(mylist[0][0], list):
            return mylist[0]
        else:
            return mylist

    if isinstance(mylist[0][0], list):
        for item in mylist[0]:
            # Here the recursion happens, create lists starting with
            # each element of the first item and completed with 
            # the rest expanded
            for rest in expand_list_list(mylist[1:]):
                reslist = copy.copy(item)
                reslist.extend(rest)
                res.append(reslist)
    else:
        for rest in expand_list_list(mylist[1:]):
            reslist = copy.copy(mylist[0])
            reslist.extend(rest)
            res.append(reslist)


    return res

