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
 
"""Classes for generation of color-ordered matrix elements using
Behrends-Giele currents.

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
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.color_ordering.color_ordered_amplitudes as \
       color_ordered_amplitudes
from madgraph import MadGraph5Error
#from madgraph import tracker

logger = logging.getLogger('madgraph.color_ordered_amplitudes')


# Dictionary from spin number to name used in the wavefunction sum functions
spin_dict = {1:'S', 2:'F', 3:'V', 4:'T'}

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
        
        assert  isinstance(amplitude,
                           color_ordered_amplitudes.ColorOrderedFlow), \
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
        co_wavefunctions = helas_objects.HelasWavefunctionList(\
            [COHelasWavefunction(wf) for wf in self.get_all_wavefunctions()])

        # Same thing for amplitudes
        for diag in self.get('diagrams'):
            diag.set('amplitudes', helas_objects.HelasAmplitudeList(\
                [COHelasAmplitude(amp) for amp in diag.get('amplitudes')]))
            # Remove old wavefunctions to free some memory
            diag.set('wavefunctions', helas_objects.HelasWavefunctionList())

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
            # Clean wfs to reduce memory
            wf.set('compare_array',[])
            wf.set('current_array', [])
            wf.set('external_numbers', array.array('I'))
            wf.set('color_string', color.ColorString())
            wf.set('lastleg_number', 0)
            wf.set('external_fermion_numbers', [])
            wf.set('all_orders', {})
            

        # Set amplitude number for all amplitudes
        for i, amp in enumerate(self.get_all_amplitudes()):
            amp.set('number', i+1)
            # Clean amps to reduce memory
            amp.set('compare_array', [])

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
            lastleg = color_ordered_amplitudes.ColorOrderedLeg(\
                                       base_vertex.get('legs').pop(-1))
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
        self['include_all_t'] = False
        self['tch_depth'] = 0
        self['identify_depth'] = 0

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

        elif name in ['permutations']:
            if not isinstance(value, list):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list" % str(value)
        
        elif name in ['include_all_t']:
            if not isinstance(value, bool):
                raise self.PhysicsObjectError, \
                        "%s is not a valid bool" % str(value)
        
        elif name in ['tch_depth', 'identify_depth']:
            if not isinstance(value, int):
                raise self.PhysicsObjectError, \
                        "%s is not a valid integer" % str(value)
        
        else:
            super(COHelasMatrixElement, self).filter(name, value)

        return True
    
    def __init__(self, amplitude = None, optimization = 3, decay_ids = [],
                 gen_color = 1, gen_periferal_diagrams = False,
                 include_all_t = True, tch_depth = 1, identify_depth = 1):
        """Initialize a COHelasMatrixElement with a ColorOrderedAmplitude"""
        
        if amplitude != None:
            if isinstance(amplitude,
                          color_ordered_amplitudes.ColorOrderedAmplitude):
                super(COHelasMatrixElement, self).__init__()
                self.get('processes').append(amplitude.get('process'))
                self.set('has_mirror_process',
                         amplitude.get('has_mirror_process'))
                for iflow, flow in enumerate(amplitude.get('color_flows')):
                    if gen_color == 1 and \
                            flow.get('process').get('orders')['singlet_QCD'] > 0:
                        continue
                    self.get('color_flows').append(COHelasFlow(flow,
                                                   optimization = optimization,
                                                   gen_color = False,
                                                   decay_ids = decay_ids,
                                                   number = iflow + 1))
                if self.get('color_flows'):
                    self.set('permutations',
                             self.get('color_flows')[0].get('permutations'))
                if gen_color and not self.get('color_basis'):
                    self.set_color_basis()
                if gen_color and not self.get('color_matrix'):
                    self.build_color_matrix(gen_color)
                if gen_periferal_diagrams:
                    self.set('include_all_t', include_all_t)
                    self.set('tch_depth', tch_depth)
                    self.set('identify_depth', identify_depth)
                    self.generate_periferal_diagrams(amplitude, decay_ids)
            else:
                # In this case, try to use amplitude as a dictionary
                super(COHelasMatrixElement, self).__init__(amplitude)
        else:
            super(COHelasMatrixElement, self).__init__()
        
    def set_color_basis(self):
        """Set the color basis using only the leading color flows."""

        max_Nc = max([flow.get('color_string').Nc_power \
                      for flow in self.get('color_flows')])

        list_color_dict = []
        # Build a color basis from the leading color flows
        for iflow, flow in enumerate(self.get('color_flows')):
            if flow.get('color_string').Nc_power < max_Nc: continue
            color_string = flow.get('color_string')
            list_color_dict.append(dict([((0,), copy.copy(color_string))]))

        self.get('color_basis')._list_color_dict = list_color_dict
        # Build the color basis
        self.get('color_basis').build()

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

        # Check if this is an all-octet amplitude
        process = self.get('processes')[0]
        model = process.get('model')
        colors = [model.get_particle(l.get('id')).get('color') for l in \
                  process.get('legs') if \
                  model.get_particle(l.get('id')).get('color') != 1]
        all_octets = (set(colors) == set([8]))
        
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
                    Nc_power = min_Nc_power - basic_Nc - col_str.Nc_power
                    # Get matrix entry
                    result, result_fixed_Nc = \
                           col_matrix.get_matrix_entry(canonical_dict,
                                                       immutable_col_str,
                                                       basic_string,
                                                       Nc_power_min = Nc_power)
                    # Ignore null entries
                    if not result:
                        continue
                    # Set color matrix entry, if non-zero
                    irow = ibasic
                    icol = iperm * len(self.get('color_flows')) + iflow

                    col_matrix[(icol, irow)] = copy.copy(result)

                    # For color octets, pick leading order contribution,
                    # and multiply by 1-1/Nc^2
                    if all_octets and gen_color == 1:
                        col_matrix[(icol, irow)][:] = \
                                      [copy.copy(result[0]),
                                       copy.copy(result[0])]
                        col_matrix[(icol, irow)][1].Nc_power -= 2
                        col_matrix[(icol, irow)][1].coeff *= -1

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
    
    def get_comp_list(self):
        """Return a list of numbers that are identical for particles
        that are permutated"""
        
        process = self.get('processes')[0]
        return color_ordered_amplitudes.ColorOrderedAmplitude.\
                                              get_comp_list(process)

    def generate_periferal_diagrams(self, amplitude, decay_ids):
        """Generate wavefunctions and amplitudes for all periferal
        diagrams for the color ordered amplitude, for use in phase
        space integration"""

        # Get all periferal diagrams as well as a list of color flows and
        # permutations that contribute to each diagram from the amplitude
        periferal_diagrams, tch_depth= \
                            amplitude.get_periferal_diagrams_from_flows(\
                             include_all_t = self.get('include_all_t'),
                             tch_depth = self.get('tch_depth'),
                             identify_depth = self.get('identify_depth'))
        self.set('tch_depth', tch_depth)
        periferal_amplitude = diagram_generation.Amplitude()
        periferal_amplitude.set('process', amplitude.get('process'))
        periferal_amplitude.set('diagrams', periferal_diagrams)

        self.generate_helas_diagrams(periferal_amplitude,
                                     optimization = 1,
                                     decay_ids = decay_ids)

    # Comparison between different amplitudes, to allow check for
    # identical processes. Note that we are then not interested in
    # interaction id, but in all other properties.
    def __eq__(self, other):
        """Comparison between different matrix elements, to allow check for
        identical processes.
        """

        if not isinstance(other, COHelasMatrixElement):
            return False

        # If no processes, this is an empty matrix element
        if not self['processes'] and not other['processes']:
            return True

        # Should only check if diagrams and process id are identical
        # Except in case of decay processes: then also initial state
        # must be the same
        if self['processes'] and not other['processes'] or \
               self['has_mirror_process'] != other['has_mirror_process'] or \
               self['processes'] and \
               self['processes'][0]['id'] != other['processes'][0]['id'] or \
               self['processes'][0]['is_decay_chain'] or \
               other['processes'][0]['is_decay_chain'] or \
               self['identical_particle_factor'] != \
                           other['identical_particle_factor'] or \
               self['diagrams'] != other['diagrams']:
            return False

        # Check if color flows are identical
        return self['color_flows'] == other['color_flows']

    def get_used_lorentz(self):
        """Return the used lorentz structures for this matrix element
        and all color flows"""

        return sum([f.get_used_lorentz() for f in self.get('color_flows')],
                   super(COHelasMatrixElement, self).get_used_lorentz())

    def get_used_couplings(self):
        """Return the used couplings structures for this matrix element
        and all color flows"""

        return sum([f.get_used_couplings() for f in self.get('color_flows')],
                   super(COHelasMatrixElement, self).get_used_couplings())

#===============================================================================
# DiagramTag class to compare color-ordered matrix elements
#===============================================================================

class IdentifyCOMETag(helas_objects.IdentifyMETag):
    """Combine IdentifyMETags for the matrix element (integration diagrams) 
       and all color flows."""

    @staticmethod
    def create_tag(amplitude, identical_particle_factor = 0):
        """Create a tag which identifies identical matrix elements"""

        assert(isinstance(amplitude, 
                          color_ordered_amplitudes.ColorOrderedAmplitude))

        process = amplitude.get('process')
        ninitial = process.get_ninitial()
        model = process.get('model')

        # Create tag for main amplitude
        tag = helas_objects.IdentifyMETag.create_tag(amplitude, 
                                                     identical_particle_factor)

        # Add tags for all color flows
        for flow in amplitude.get('color_flows'):
            tag.append(sorted([helas_objects.IdentifyMETag(d, model, ninitial) \
                                   for d in amplitude.get('diagrams')]))
        tag.sort()
        return tag

#===============================================================================
# COHelasMultiProcess
#===============================================================================
class COHelasMultiProcess(helas_objects.HelasMultiProcess):
    """Version of HelasMultiProcess which generates COHelasMatrixElements"""

    # Matrix element class
    matrix_element_class = COHelasMatrixElement
    # Tag class to compare matrix elements
    identify_tag_class = IdentifyCOMETag

    def __init__(self, argument=None, gen_color = 1, optimization = 3,
                 gen_periferal_diagrams = False,
                 include_all_t = True, tch_depth = 10, identify_depth = 1):
        """Allow initialization with AmplitudeList"""

        if isinstance(argument,
                      color_ordered_amplitudes.ColorOrderedMultiProcess):
            super(COHelasMultiProcess, self).__init__()
            self.set('matrix_elements',
                     self.generate_matrix_elements(argument.get('amplitudes'),
                                                   gen_color, optimization,
                                                   gen_periferal_diagrams,
                                                   include_all_t, tch_depth,
                                                   identify_depth))
        if isinstance(argument, diagram_generation.AmplitudeList):
            super(COHelasMultiProcess, self).__init__()
            self.set('matrix_elements',
                     self.generate_matrix_elements(argument,
                                                   gen_color = gen_color,
                                                   optimization = optimization,
                                                   gen_periferal_diagrams = \
                                                   gen_periferal_diagrams,
                                                   include_all_t = include_all_t,
                                                   tch_depth = tch_depth,
                                                   identify_depth = identify_depth))
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
                                 optimization = 1, decay_ids = [],
                                 gen_periferal_diagrams = False,
                                 include_all_t = True, tch_depth = 1,
                                 identify_depth = 1):
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

        # List of valid matrix elements
        matrix_elements = helas_objects.HelasMatrixElementList()
        # List of identified matrix_elements
        identified_matrix_elements = []
        # List of amplitude tags, synchronized with identified_matrix_elements
        amplitude_tags = []
        # List of the external leg permutations for the amplitude_tags,
        # which allows to reorder the final state particles in the right way
        # for maximal process combination
        permutations = []

        while amplitudes:
            # Pop the amplitude to save memory space
            amplitude = amplitudes.pop(0)
            if isinstance(amplitude, diagram_generation.DecayChainAmplitude):
                matrix_element_list = helas_objects.HelasDecayChainProcess(amplitude).\
                                      combine_decay_chain_processes()
            else:
                # Create tag identifying the matrix element using
                # IdentifyMETag. If two amplitudes have the same tag,
                # they have the same matrix element
                amplitude_tag = cls.identify_tag_class.create_tag(amplitude)
                try:
                    me_index = amplitude_tags.index(amplitude_tag)
                except ValueError:
                    logger.info("Generating Helas calls for color ordered %s" % \
                         amplitude.get('process').nice_string().\
                                           replace('Process', 'process'))
                    # Create matrix element for this amplitude
                    matrix_element_list = [cls.matrix_element_class(amplitude,
                                                   decay_ids=decay_ids,
                                                   gen_color=gen_color,
                                                   optimization = optimization,
                                                   gen_periferal_diagrams = \
                                                          gen_periferal_diagrams,
                                                   include_all_t = include_all_t,
                                                   tch_depth = tch_depth,
                                                   identify_depth = \
                                                                 identify_depth)]
                    me = matrix_element_list[0]
                    if me.get('processes') and me.get('color_flows'):
                        # Keep track of amplitude tags
                        amplitude_tags.append(amplitude_tag)
                        identified_matrix_elements.append(me)
                        permutations.append(amplitude_tag[-1][0].\
                                            get_external_numbers())
                else:
                    # Identical matrix element found
                    other_processes = identified_matrix_elements[me_index].\
                                      get('processes')
                    logger.info("Combining %s with %s" % \
                                (amplitude.get('process').nice_string().\
                                     replace('Process: ', ''),
                                 other_processes[0].nice_string().\
                                     replace('Process: ', '')))
                    other_processes.append(cls.reorder_process(\
                        amplitude.get('process'),
                        permutations[me_index],
                        amplitude_tag[-1][0].get_external_numbers()))
                    # Go on to next amplitude
                    continue
            # Deal with newly generated matrix element
            for matrix_element in matrix_element_list:
                assert isinstance(matrix_element, helas_objects.HelasMatrixElement), \
                          "Not a HelasMatrixElement: %s" % matrix_element

                # If the matrix element has no diagrams,
                # remove this matrix element.
                if not matrix_element.get('processes') or \
                       not matrix_element.get('color_flows'):
                    continue
                # Otherwise, add this matrix element to list
                matrix_elements.append(matrix_element)

        return matrix_elements

#===============================================================================
# ColorOrderedSubProcessGroup
#===============================================================================

class COSubProcessGroup(group_subprocs.SubProcessGroup):
    """Class to group a number of color ordered amplitudes with same
    initial states into a subprocess group"""

    helas_multi_process_class = COHelasMultiProcess

    @staticmethod
    def group_amplitudes(cls, amplitudes, gen_color, optimization,
                         gen_periferal_diagrams,
                         include_all_t = True, tch_depth = 10, identify_depth = 1):
        """Return a SubProcessGroupList with the color ordered amplitudes divided
        into subprocess groups, and store the additional options"""

        groups = group_subprocs.SubProcessGroup.group_amplitudes(cls,
                                                                 amplitudes)
        for group in groups:
            group.gen_color = gen_color
            group.optimization = optimization
            group.gen_periferal_diagrams = gen_periferal_diagrams
            group.include_all_t = include_all_t
            group.tch_depth = tch_depth
            group.identify_depth = identify_depth

        return groups

    #===========================================================================
    # generate_matrix_elements
    #===========================================================================
    def generate_matrix_elements(self):
        """Create a COHelasMultiProcess corresponding to the amplitudes
        in self"""

        if not self.get('amplitudes'):
            raise self.PhysicsObjectError, \
                  "Need amplitudes to generate matrix_elements"

        amplitudes = copy.copy(self.get('amplitudes'))

        self.set('matrix_elements',
                 self.helas_multi_process_class.\
                         generate_matrix_elements(amplitudes,
                         gen_color = self.gen_color,
                         optimization = self.optimization,
                         gen_periferal_diagrams = self.gen_periferal_diagrams,
                         include_all_t = self.include_all_t,
                         tch_depth = self.tch_depth,
                         identify_depth = self.identify_depth))

        self.set('amplitudes', diagram_generation.AmplitudeList())

