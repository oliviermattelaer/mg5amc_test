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

"""Definitions of the objects needed both for MadFKS from real 
and MadFKS from born"""

import madgraph.core.base_objects as MG
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.color_amp as color_amp
import madgraph.core.color_algebra as color_algebra
from operator import itemgetter
import copy
import logging
import array
import fractions
    
    
    
class FKSProcessError(Exception):
    """Exception for MadFKS"""
    pass


class FKSDiagramTag(diagram_generation.DiagramTag): #test written
    """Modified diagram tags to be used to link born and real configurations.
    """

    @staticmethod
    def link_from_leg(leg, model):
        """Returns the default end link for a leg: ((id, number), number).
        Note that the number is not taken into account if tag comparison,
        but is used only to extract leg permutations.
        """
        return [((leg.get('id'), leg.get('number')), leg.get('number'))]


def link_rb_conf(born_amp, real_amp, i, j, ij): #test written, sm and heft
    """finds the real configurations that match the born ones, i.e.
    for each born configuration, the real configuration that has 
    the ij -> i j splitting.
    i, j and ij are integers, and refer to the leg position in the real
    process (i, j) and in the born process (ij).
    """
    # find diagrams with 3 point functions and use them as configurations

    minvert = min([max([len(vert.get('legs')) \
                        for vert in diag.get('vertices')]) \
                        for diag in born_amp.get('diagrams')])

    born_confs = []
    real_confs = []

    for k, diag in enumerate(born_amp.get('diagrams')):
        if any([len(vert.get('legs')) > minvert for vert in
                diag.get('vertices')]):
            continue
        else:
            born_confs.append({'number' : k, 'diagram' : diag})

    for k, diag in enumerate(real_amp.get('diagrams')):
        if any([len(vert.get('legs')) > minvert \
                for vert in diag.get('vertices')]):
            continue
        else:
            real_confs.append({'number': k, 'diagram': diag})

    good_diags = []

    # find the real diagrams that have i and j attached to the same vertex
    real_confs_new = copy.deepcopy(real_confs)
    for diag in real_confs_new:
        for vert in diag['diagram'].get('vertices'):
            vert_legs = [ l.get('number') for l in vert.get('legs')]
            if (i in vert_legs and not j in vert_legs) or \
               (j in vert_legs and not i in vert_legs):
                   break

            if i in vert_legs and j in vert_legs:
                vert_legs.remove(i)
                vert_legs.remove(j)
                last_leg = vert_legs[0]
                diag['diagram']['vertices'].remove(vert)
                good_diags.append({'diagram': diag['diagram'], 
                                  'leg_ij': last_leg,
                                  'number': diag['number']})
                break #no need to continue once the vertex is found


    # now good_diags contains the real_confs which had the splitting, 
    #  with the vertex corresponding to the splitting removed

    # The legs in the born and real diags are ordered according to 
    #  the same criterion. Once we removed i-j, we have to re-label the
    #  real legs to match the born numbering.

    legs = []
    for d in good_diags: 
        for v in d['diagram'].get('vertices'):
            for l in v.get('legs'):
                if l not in legs:
                    legs.append(copy.copy(l))

    for good_diag in good_diags:
        replaced_ij = False
        for vert in good_diag['diagram'].get('vertices'):
            for l in vert.get('legs'):
                shift = 0
                #shift the legs
                if l.get('number') > 0 and l in legs:
                    if l.get('number') > j:
                        shift += 1
                    if l.get('number') > i:
                        shift += 1
                    if l.get('number') > ij and ij < max(i,j):
                        shift -= 1
                    if l.get('number') != good_diag['leg_ij'] or replaced_ij:
                        legs.remove(l)
                        l['number'] -= shift
                 #and relabel last_leg to ij
                    if l.get('number') == good_diag['leg_ij'] and not replaced_ij and l in legs:
                        legs.remove(l)
                        replaced_ij = True
                        l['number'] = ij

    # now create the tags
    born_tags = [FKSDiagramTag(d['diagram'], 
                               born_amp.get('process').get('model')) \
                               for d in born_confs]


    real_tags = [FKSDiagramTag(d['diagram'], 
                               real_amp.get('process').get('model')) \
                               for d in good_diags ]

    # and compare them
    if len(born_tags) != len(real_tags):
        raise FKSProcessError('Cannot map born/real configurations between \
                %s and %s: not same number of configurations: %d %d' % \
                (born_amp.get('process').nice_string().replace('Process:',''), 
                 real_amp.get('process').nice_string().replace('Process:',''),
                               len(born_tags),
                               len(real_tags)))
    
    links = []
    for ib, btag in enumerate(born_tags):
        try:
            ir = real_tags.index(btag)
            links.append({'real_conf': good_diags[ir]['number'],
                          'born_conf': born_confs[ib]['number']})
            real_tags.remove(btag)
            good_diags.pop(ir)
        except ValueError:
            print real_tags, i, j, ij
            print '\n'.join( d['diagram'].nice_string() for d in good_diags)
            raise FKSProcessError('Linking %s to %s: could not link born diagram %s' % \
                 (born_amp.get('process').nice_string().replace('Process:',''), 
                  real_amp.get('process').nice_string().replace('Process:',''),
                                  born_confs[ib]['diagram'].nice_string()) )

    return links


def find_orders(amp): #test_written
    """Takes an amplitude as input, and returns a dictionary with the
    orders of the couplings.
    """
    assert isinstance(amp, diagram_generation.Amplitude)
    orders = {}
    for diag in amp.get('diagrams'):
        for order, value in diag.get('orders').items():
            try:
                orders[order] = max(orders[order], value)
            except KeyError:
                orders[order] = value
    return orders


def find_splittings(leg, model, dict, pert='QCD'): #test written
    """Finds the possible splittings corresponding to leg
    """
    if dict == {}:
        dict = find_pert_particles_interactions(model, pert)
    splittings = []
#check that the leg is a qcd leg

    if leg.get('id') in dict['pert_particles']:
        part = model.get('particle_dict')[leg.get('id')]
        antipart = model.get('particle_dict')[part.get_anti_pdg_code()]
        for ii in dict['interactions']:
#check which interactions contain leg and at least one soft particles:
            parts = copy.deepcopy(ii['particles'])
            nsoft = 0
            if part in parts:
                #pops the ANTI-particle of part from the interaction
                parts.pop(parts.index(antipart))
                for p in parts:
                    if p.get_pdg_code() in dict['soft_particles']:
                        nsoft += 1
                if nsoft >= 1:
                    splittings.extend(split_leg(leg, parts, model))
    return splittings


def split_leg(leg, parts, model): #test written
    """Splits the leg into parts, and returns the two new legs.
    """
    #for an outgoing leg take the antiparticles
    split = []
    #for a final state particle one can have only a splitting
    if leg['state'] :
        split.append([])
        for part in parts:
            split[-1].append(to_fks_leg({'state': True, \
                                 'id': part.get_pdg_code()},model))
            ij_final(split[-1])
    #while for an initial state particle one can have two splittings 
    # if the two partons are different
    else:
        if parts[0] != parts[1]:
            for part in parts:
                cparts = copy.deepcopy(parts)
                split.append([\
                          to_fks_leg({'state': False, 
                                  'id': cparts.pop(cparts.index(part)).get_pdg_code(),
                                  'fks': 'j'}, model),
                          to_fks_leg({'state': True,
                                  'id': cparts[0].get_anti_pdg_code(),
                                  'fks': 'i'}, model)\
                          ])
        else:
            split.append([\
                            to_fks_leg({'state': False, 
                                  'id': parts[0].get_pdg_code(),
                                  'fks': 'j'}, model),
                            to_fks_leg({'state': True, 
                                  'id': parts[1].get_anti_pdg_code(),
                                  'fks': 'i'}, model)])
    return split


def ij_final(pair):
    """given a pair of legs in the final state, assigns the i/j fks id
    NOTE: the j partons is always put before the i one
    """
    #if a massless bosonic particle is in the pair, it is i
    #else by convention the anti-particle is labeled i
    #the order of the splitting is [ j, i]
    if len(pair) == 2:
        for i in range(len(pair)):
            set = 0
            if (pair[i]['massless'] and pair[i]['spin'] %2 == 1) or \
               (pair[i]['color'] == -3 and pair[1-i]['color'] == 3) and \
               not set:
                pair[i]['fks'] = 'i'
                pair[1-i]['fks'] = 'j'
                #check that first j then i
                if i < 1 - i:
                    pair.reverse()
                set = 1

def insert_legs(leglist_orig, leg, split):
    """Returns a new leglist with leg splitted into split.
    The convention is to remove leg ij, replace it with leg j, and put
    i at the end of the group of legs with the same color representation
    """
    # the deepcopy statement is crucial
    leglist = FKSLegList(copy.deepcopy(leglist_orig))         
    #find the position of the first final state leg
    for i in range(len(leglist)):
        if leglist[-i - 1].get('state'):
            firstfinal = len(leglist) - i - 1
    # replace leg with leg_j  (split[0])
    leglist[leglist.index(leg)] = split[0]
    # and find where to insert i  (split[1])
    col_maxindex = {}
    mass_col_maxindex = {}
    for col in set([l['color'] for l in leglist[firstfinal:] if l['massless']]):
        col_maxindex[col] = max([0] + [leglist.index(l) for l in leglist[firstfinal:] if l['color'] == col and l['massless']])
    for col in set([abs(l['color']) for l in leglist[firstfinal:] if not l['massless']]):
        mass_col_maxindex[col] = max([0] + [leglist.index(l) for l in leglist[firstfinal:] if abs(l['color']) == col and not l['massless']])
    #no need to keep info on particles with color > i
    for col in copy.copy(col_maxindex.keys()):
        if abs(col) > abs(split[1]['color']):
            del col_maxindex[col]
    for col in copy.copy(mass_col_maxindex.keys()):
        if abs(col) > abs(split[1]['color']):
            del mass_col_maxindex[col]
    #also remove antiquarks if i is a quark
    if split[1]['color'] > 0:
        try:
            del col_maxindex[-split[1]['color']]
        except KeyError:
            pass
    #so now the maximum of the max_col entries should be the position to insert leg i

    leglist.insert(max(col_maxindex.values() + mass_col_maxindex.values() + [firstfinal - 1] ) + 1, split[1])
#    for sleg in split:            
#        leglist.insert(i, sleg)
#        #keep track of the number for initial state legs
#        #if not sleg.get('state') and not leg.get('state'):
#        leglist[i]['number'] = leg['number']
#        i += 1
#        if i < firstfinal:
#            i = firstfinal
#            
#    leglist.sort()
    for i, leg in enumerate(leglist):
        leg['number'] = i + 1        
    return leglist 


def combine_ij( i, j, model, dict, pert='QCD'): #test written
    """checks whether FKSlegs i and j can be combined together in the given model
    and with given perturbation order and if so combines them into ij. 
    If dict is empty it is initialized with find_pert_particles_interactions
    """
    if dict == {}:
        dict = find_pert_particles_interactions(model, pert)
    ij = []
    num = copy.copy(min(i.get('number'), j.get('number')))
    
    # we do not want j being a massiless vector unless also i is or j is initial
    not_double_counting = (j.get('spin') == 3 and j.get('massless') and 
                           i.get('spin') == 3 and i.get('massless')) or \
                           j.get('spin') != 3 or not j.get('massless') or \
                           not j.get('state')

    #if i and j are a final state particle and antiparticle pair,
    # then we want i to be antipart and j to be 
    if j.get('state') and j.get('id') == - i.get('id'):  
        not_double_counting = not_double_counting and j.get('id') >0
                          
    if i.get('id') in dict['soft_particles'] and \
       j.get('id') in dict['pert_particles'] and \
       i.get('state') and not_double_counting:
        for int in dict['interactions']:
            parts= copy.copy(int['particles'])
                #remove i
            try:
                parts.remove(model.get('particle_dict')[i.get('id')])
            except ValueError:
                continue

            #remove j if final state, anti j if initial state
            if j.get('state'):
                j_id = j.get('id')
            else:
                j_id = model.get('particle_dict')[j.get('id')].get_anti_pdg_code()
            try:
                parts.remove(model.get('particle_dict')[j_id])
            except ValueError:
                continue
            
            #ij is what remains if j is initial, the anti of if j is final
            if j.get('state'):
                ij.append(MG.Leg({
                    'id': parts[0].get_anti_pdg_code(),
                    'state': True,
                    'number': num}))
            else:
                ij.append(MG.Leg({
                    'id': parts[0].get_pdg_code(),
                    'state': False,
                    'number': num}))
    return to_fks_legs(ij, model)       


def find_pert_particles_interactions(model, pert_order = 'QCD'): #test written
    """given a model and pert_order, returns a dictionary with as entries:
    --interactions : the interactions of order pert_order
    --pert_particles : pdgs of particles taking part to interactions
    --soft_particles : pdgs of massless particles in pert_particles
    """
    ghost_list = [-100033, 100033]
    ghost_list += [ p['get_pdg_code'] for p in model.get('particles') if p.get('spin') < 0]
    qcd_inter = MG.InteractionList()
    pert_parts = []
    soft_parts = []
    for i, ii in model.get('interaction_dict').items():
        # i want interections of pert_order: 1 (from LO to NLO), 
        # without any other orders
        if ii.get('orders') == {pert_order:1} and len(ii['particles']) == 3 :
            masslist = [p.get('mass').lower() for p in ii['particles']]
                # check that there is at least a massless particle, and that the 
                # remaining ones have the same mass 
                # (otherwise the real emission final state will not be degenerate
                # with the born one
            try:
                masslist.remove('zero')
            except:
                continue
            if len(set(masslist)) == 1 and not \
                    any( [ p['pdg_code'] in ghost_list for p in ii['particles']]) :
                qcd_inter.append(ii)
                for pp in ii['particles']:
                    pert_parts.append(pp.get_pdg_code())
                    if pp['mass'].lower() == 'zero':
                        soft_parts.append(pp.get_pdg_code())

    return {'interactions': sorted(qcd_inter), 
            'pert_particles': sorted(set(pert_parts)),
            'soft_particles': sorted(set(soft_parts))}    


def insert_color_links(col_basis, col_obj, links): #test written
    """insert the color links in col_obj: returns a list of dictionaries
    (one for each link) with the following entries:
    --link: the numbers of the linked legs
    --link_basis: the linked color basis
    --link_matrix: the color matrix created from the original basis and the linked one
    """   
    assert isinstance(col_basis, color_amp.ColorBasis)
    assert isinstance(col_obj, list)
    result =[]
    for link in links:
        this = {}
        #define the link
        l =[]
        for leg in link['legs']:
            l.append(leg.get('number'))
        this['link'] = l
            
        #replace the indices in col_obj of the linked legs according to
        #   link['replacements']
        # and extend-> product the color strings
            
        this_col_obj = []
        for old_dict in col_obj:
            dict = copy.copy(old_dict)
            for k, string in dict.items():
                dict[k] = string.create_copy()
                for col in dict[k]:
                    for ind in col:
                        for pair in link['replacements']:
                            if ind == pair[0]:
                                col[col.index(ind)] = pair[1]
                dict[k].product(link['string'])
            this_col_obj.append(dict)
        basis_link = color_amp.ColorBasis()
        for ind, dict in enumerate(this_col_obj):
            basis_link.update_color_basis(dict, ind)
   
        this['link_basis'] = basis_link
        this['link_matrix'] = color_amp.ColorMatrix(col_basis,basis_link)               
        result.append(this)
    basis_orig = color_amp.ColorBasis()
    for ind, dict in enumerate(col_obj):
            basis_orig.update_color_basis(dict, ind)
    
    for link in result:
        link['orig_basis'] = basis_orig
    return result



def find_color_links(leglist, symm = False): #test written
    """Finds all the possible color links between any 
    two legs of the born.
    If symm is true, only half of the color links are generated, those
    for which leg1['number'] <= leg2['number']
    """

    color_links = []
    for leg1 in leglist:
        for leg2 in leglist:
            #legs must be colored and different, unless massive
                if (leg1.get('color') != 1 and leg2.get('color') != 1) \
                  and (leg1 != leg2 or not leg1.get('massless')):
                    if not symm or leg1['number'] <= leg2['number']:
                        col_dict = legs_to_color_link_string(leg1,leg2)
                        color_links.append({
                            'legs': [leg1, leg2],
                            'string': col_dict['string'],
                            'replacements': col_dict['replacements']})
    return color_links
             

def legs_to_color_link_string(leg1, leg2): #test written, all cases
    """given two FKSlegs, returns a dictionary containing:
    --string: the color link between the two particles, to be appended to
        the old color string
        extra minus or 1/2 factor are included as it was done in MadDipole
    --replacements: a pair of lists containing the replacements of the color 
        indices in the old string to match the link 
    """
    #the second-to-last index of the t is the triplet,
    # the last is the anti-triplet

    legs = FKSLegList([leg1, leg2]) 
    dict = {}
    min_index = -3000
    iglu = min_index*2
    string = color_algebra.ColorString()
    replacements = []
    if leg1 != leg2:
        for leg in legs:
            min_index -= 1
            num = leg.get('number')
            replacements.append([num, min_index])
            icol = 1
            if not leg.get('state'):
                icol = - 1
            if leg.get('color') * icol == 3:
                string.product(color_algebra.ColorString([
                               color_algebra.T(iglu, num, min_index)]))
                string.coeff = string.coeff * (-1)
            elif leg.get('color') * icol == - 3:
                string.product(color_algebra.ColorString([
                               color_algebra.T(iglu, min_index, num)]))
            elif leg.get('color') == 8:
                string.product(color_algebra.ColorString(init_list = [
                               color_algebra.f(min_index,iglu,num)], 
                               is_imaginary =True))

    else:
        icol = 1
        if not leg1.get('state'):
            icol = - 1
        num = leg1.get('number')
        replacements.append([num, min_index -1])
        if leg1.get('color') * icol == 3:
            string = color_algebra.ColorString(
                     [color_algebra.T(iglu, iglu, num, min_index -1)])
        elif leg1.get('color') * icol == - 3:
            string = color_algebra.ColorString(
                     [color_algebra.T(iglu, iglu, min_index-1, num)])
        string.coeff = string.coeff * fractions.Fraction(1, 2) 
    dict['replacements'] = replacements
    dict['string'] = string
      
    return dict


def sort_proc(process):
    """Given a process, this function returns the same process 
    but with sorted FKSLegs.
    """
    leglist = to_fks_legs(process.get('legs'), process.get('model'))
    leglist.sort()
    for n, leg in enumerate(leglist):
        leg['number'] = n + 1
    process['legs'] = leglist

    return process


def to_leg(fksleg):
    """Given a FKSLeg, returns the original Leg.
    """
    leg = MG.Leg( \
        {'id': fksleg.get('id'),
         'number': fksleg.get('number'),
         'state': fksleg.get('state'),
         'from_group': fksleg.get('from_group'),
          })
    return leg


def to_legs(fkslegs):
    """Given a FKSLegList, returns the corresponding LegList.
    """
    leglist = MG.LegList()
    for leg in fkslegs:
        leglist.append(to_leg(leg))
    return leglist


def to_fks_leg(leg, model): #test written
    """Given a leg or a dict with Leg entries, 
    adds color, spin and massless entries, according to model"""
    fksleg = FKSLeg(leg)
    part = model.get('particle_dict')[leg['id']]
    fksleg['color'] = part.get_color()
    fksleg['massless'] = part['mass'].lower() == 'zero'
    fksleg['spin'] = part.get('spin')      
    return fksleg

    
def to_fks_legs(leglist, model): #test written
    """given leglist, sets color and massless entries according to the model 
    variable.
    return a FKSLeglist"""
    fkslegs = FKSLegList()     
    for leg in leglist:
        fkslegs.append(to_fks_leg(leg, model))
    return fkslegs     


class FKSLegList(MG.LegList):
    """list of FKSLegs"""
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid FKSLeg for the list."""
        return isinstance(obj, FKSLeg)

    def sort(self):
        """Sorting routine, sorting chosen to be optimal for madfks"""
        sorted_leglist = FKSLegList()
        #find initial state legs
        initial_legs = FKSLegList([l for l in copy.copy(self) if not l['state']])
        #find final state legs
        final_legs = FKSLegList([l for l in copy.copy(self) if l['state']])
        if len(initial_legs) == 1:
            sorted_leglist.extend(initial_legs)
        elif len(initial_legs) == 2:
            if initial_legs[0]['number'] > initial_legs[1]['number']:
                initial_legs.reverse()
            sorted_leglist.extend(initial_legs)
        else: 
            raise FKSProcessError('Too many initial legs')
        #find color representations
        colors = sorted(set([abs(l['color']) for l in final_legs]))
        for col in colors:
            col_legs = FKSLegList([l for l in final_legs if abs(l['color']) == col])
            #find massive and massless legs in this color repr
            massive_legs = [l for l in col_legs if not l['massless']]
            massless_legs = [l for l in col_legs if l['massless']]
            # sorting may be different for massive and massless particles
            # for color singlets, do not change order
            if col == 1:
                keys = [itemgetter('number'), itemgetter('number')]
                reversing = False
            else:
                keys = [itemgetter('id'), itemgetter('id')]
                reversing = True

            for i, list in enumerate([massive_legs, massless_legs]):
                init_pdg_legs = []
                if len(initial_legs) == 2:
                #put first legs which have the same abs(pdg) of the initial ones
                    for i in range(len(set([ abs(l['id']) for l in initial_legs]))):
                        pdg = abs(initial_legs[i]['id'])
                        init_pdg_legs = [l for l in list if abs(l['id']) == pdg]
                        if init_pdg_legs:
                            # sort in order to put first quarks then antiparticles,
                            #  and to put fks partons as n j i
                            init_pdg_legs.sort(key = keys[i], reverse=reversing)
                            sorted_leglist.extend(FKSLegList(init_pdg_legs))

                    init_pdgs = [ abs(l['id']) for l in initial_legs]
                    other_legs = [l for l in list if not abs(l['id']) in init_pdgs]
                    other_legs.sort(key = keys[i], reverse=reversing)
                    sorted_leglist.extend(FKSLegList(other_legs))
                else:
                    list.sort(key = keys[i], reverse=reversing)
                    sorted_leglist.extend(FKSLegList(list))

        for i, l in enumerate(sorted_leglist):
            self[i] = l




class FKSLeg(MG.Leg):
    """a class for FKS legs: it inherits from the ususal leg class, with two
    extra keys in the dictionary: 
    -'fks', whose value can be 'i', 'j' or 'n' (for "normal" particles) 
    -'color', which gives the color of the leg
    -'massless', boolean, true if leg is massless
    -'spin' which gives the spin of leg
    """

    def default_setup(self):
        """Default values for all properties"""
        super(FKSLeg, self).default_setup()

        self['fks'] = 'n'
        self['color'] = 0
        self['massless'] = True
        self['spin'] = 0
    
    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""
        keys = super(FKSLeg, self).get_sorted_keys()
        keys += ['fks', 'color', 'massless', 'spin']
        return keys

    
    def filter(self, name, value):
        """Filter for valid leg property values."""

        if name == 'fks':
            if not isinstance(value, str):
                raise self.PhysicsObjectError, \
                        "%s is not a valid string for leg fks flag" \
                                                        % str(value)
        if name in ['color', 'spin']:
            if not isinstance(value, int):
                raise self.PhysicsObjectError, \
                        "%s is not a valid leg color " % \
                                                 str(value)
        if name == 'massless':
            if not isinstance(value, bool):
                raise self.PhysicsObjectError, \
                        "%s is not a valid boolean for leg flag massless" % \
                                                                    str(value)                                                           
        return super(FKSLeg,self).filter(name, value)
    
     
