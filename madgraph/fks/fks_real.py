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

"""Definitions of the objects needed for the implementation of MadFKS"""

import madgraph.core.base_objects as MG
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.color_amp as color_amp
import madgraph.core.color_algebra as color_algebra
import madgraph.fks.fks_common as fks_common
import copy
import logging
import fractions

logger = logging.getLogger('madgraph.fks_real')
#from math import min

class FKSMultiProcessFromReals(diagram_generation.MultiProcess): #test written
    """a multi process class that contains informations on the born processes 
    and the reals"""
    
    def default_setup(self):
        """Default values for all properties"""
        super(FKSMultiProcessFromReals, self).default_setup()

        self['real_processes'] = FKSProcessFromRealsList()
    
    def get_sorted_keys(self):
        """Return particle property names as a nicely sorted list."""
        keys = super(FKSMultiProcessFromReals, self).get_sorted_keys()
        keys += ['real_processes']
        return keys

    def filter(self, name, value):
        """Filter for valid leg property values."""

        if name == 'real_processes':
            if not isinstance(value, FKSProcessFromRealsList):
                raise self.PhysicsObjectError, \
                        "%s is not a valid list for real_processes " % str(value)
                                                                
        return super(FKSMultiProcessFromReals,self).filter(name, value)
    
    def __init__(self, *arguments):
        """initialize the original multiprocess, then generates the amps for the 
        borns, then geneare the born processes and the reals"""
                
        super(FKSMultiProcessFromReals, self).__init__(*arguments)
        
        amps = self.get('amplitudes')
        born_amplist = []
        born_amp_id_list = []
        for amp in amps:
            real_proc = FKSProcessFromReals(amp)
            self['real_processes'].append(real_proc)
            
            
class FKSBornProcess(object):
    """contains informations about a born process
    -- born amplitude
    -- i/j_fks (in the real process leglist)
    -- ijglu -> 0 if ij is not a gluon, ij[number] otherwise
    -- need_color_link -> if color links are needed (ie i_fks is a gluon)
    -- color link list
    -- is_nbody_only
    -- is_to_integrate
    """
    
    def __init__(self, real_proc, leg_i, leg_j, leg_ij, perturbed_orders = ['QCD']):
        """initialize the born process starting from the real process and the
        combined leg"""

        self.i_fks = leg_i.get('number')
        self.j_fks = leg_j.get('number')
        if leg_ij.get('spin') == 3:
            self.ijglu = leg_ij.get('number')
        else:
            self.ijglu = 0
        
        self.need_color_links = leg_i.get('spin') == 3 and leg_i.get('massless')       
        self.process = copy.copy(real_proc)
#        orders = copy.copy(real_proc.get('orders'))
        
#        for order in perturbed_orders:
#            orders[order] += -1
#        self.process.set('orders', orders)
        
        born_legs = self.reduce_real_leglist(leg_i, leg_j, leg_ij)
        self.process.set('legs', MG.LegList(born_legs))
        self.amplitude = diagram_generation.Amplitude(self.process)
        self.is_to_integrate = True
        self.is_nbody_only = False
        self.color_links = []
    
    
    def reduce_real_leglist(self, leg_i, leg_j, leg_ij):
        """removes from the leglist of self.process leg_i, leg_j 
        and inserts leg_ij (fkslegs)"""
        red_leglist = fks_common.to_fks_legs(
                        copy.deepcopy(self.process.get('legs')),
                                      self.process.get('model'))

        red_leglist.remove(leg_i)
        red_leglist.remove(leg_j)
        red_leglist.insert(leg_ij.get('number')-1, leg_ij)
        for n, leg in enumerate(red_leglist):
            red_leglist[n].set('number', n+1)
        return red_leglist

    def find_color_links(self): #test written
        """finds all the possible color links between two legs of the born.
        Uses the find_color_links function in fks_common"""
        if self.need_color_links:
            self.color_links = fks_common.find_color_links(\
                              self.process.get('legs'))
        return self.color_links
    
    def get_born_fks_inc_string(self, iconf):#test written
        """returns the part of the fks.inc file relative to this born configuration
        the configuration number is given"""
        replace = {'conf_num' : iconf, 
                   'i_fks' : self.i_fks, 
                   'j_fks' : self.j_fks}
        string =" \n\
c     FKS configuration number  %(conf_num)d \n\
data fks_i(  %(conf_num)d  ) /  %(i_fks)d  / \n\
data fks_j(  %(conf_num)d  ) /  %(j_fks)d  / \n\
      " % replace
        return string
 
#===============================================================================
# FKS Process
#===============================================================================

class FKSProcessFromRealsList(MG.PhysicsObjectList):
    """class to handle lists of FKSProcesses"""
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid FKSProcessFromReals for the list."""
        return isinstance(obj, FKSProcessFromReals)
   
   
class FKSProcessFromReals(object):
    """the class for an FKS process, starting from reals """ 
    
    def __init__(self, start_proc = None, remove_borns = True): #test written
        """initialization: starts either from an amplitude or a process,
        then init the needed variables.
        remove_borns tells if the borns not needed for integration will be removed
        from the born list (mainly used for testing)
        --real_proc/real_amp
        --model
        --leglist, nlegs
        --pdg codes, colors ###, ipos, j_from_i to be written in fks.inc
        --borns"""
        
        self.borns = []
        self.leglist = []
        self.pdg_codes = []
        self.colors = []
        self.nlegs = 0
        self.fks_j_from_i = {}
        self.real_proc = None
        self.real_amp = None
        self.model = None
        self.nincoming = 0
 
        if start_proc:
            if isinstance(start_proc, MG.Process):
                self.real_proc = start_proc 
                self.real_amp = diagram_generation.Amplitude(self.real_proc)
            elif isinstance(start_proc, diagram_generation.Amplitude):
                self.real_amp = start_proc
                self.real_proc = start_proc.get('process') 
            self.model = self.real_proc.get('model')   

            self.leglist = fks_common.to_fks_legs(self.real_proc.get('legs'), self.model)
            self.nlegs = len(self.leglist)
            for leg in self.leglist:
                self.pdg_codes.append(leg['id'])
                self.colors.append(leg['color'])
                if not leg['state']:
                    self.nincoming += 1
            
            self.find_borns()
            self.find_borns_to_integrate(remove_borns)
            self.find_born_nbodyonly()
    
            
    def get_fks_inc_string(self): #test written
        """returns the list of configurations corrresponding to the various 
        underlying borns"""
        string = "integer fks_configs, ipos, jpos \n\
data fks_configs / %(nconfs)d / \n\
integer fks_i(%(nconfs)d), fks_j(%(nconfs)d) \n\
INTEGER FKS_IPOS(0:NEXTERNAL) \n\
INTEGER FKS_J_FROM_I(NEXTERNAL, 0:NEXTERNAL) \n\
INTEGER PARTICLE_TYPE(NEXTERNAL), PDG_TYPE(NEXTERNAL) \n" %{'nconfs' : 
                                                        len(self.borns)}
        n = 0
        for born in self.borns:
            if born.is_to_integrate:
                n +=1
                string += born.get_born_fks_inc_string(n)
        
        ipos = []
        for ii, js in self.fks_j_from_i.items():
            if js:
                ipos.append(ii)
                string += "\n\
data (fks_j_from_i(%d, jpos), jpos = 0, %d)  / %d, %s /" \
        % (ii, len(js), len(js), ', '.join(["%d" % j for j in js]))
        string += "\n\
data (fks_ipos(ipos), ipos = 0, %d)  / %d, %s /\n" \
        % (len(ipos), len(ipos), ', '.join(["%d" % i for i in ipos]))
        
        string += "\n\
C\n\
C     Particle type:\n\
C     octet = 8, triplet = 3, singlet = 1\n\
DATA (PARTICLE_TYPE(IPOS), IPOS=1, NEXTERNAL) / %s / \n" % \
                                ', '.join(["%d" % i for i in self.colors])
        string += "\n\
C\n\
C     Particle type according to PDG:\n\
C\n\
      DATA (PDG_TYPE(IPOS), IPOS=1, NEXTERNAL) / %s / \n " % \
                                ', '.join(["%d" % i for i in self.pdg_codes])
        return string
             
                 
    def find_borns(self): #test written
        """finds the underlying borns for a given fks real process"""
        dict ={}
        
        for i in self.leglist:
            self.fks_j_from_i[i.get('number')] = []
            if i.get('state'):
                for j in self.leglist:
                    if j.get('number') != i.get('number') :
                        ijlist = fks_common.combine_ij(i, j, self.model, dict)
                        for ij in ijlist:
                            born = FKSBornProcess(self.real_proc, i, j, ij)
                            if born.amplitude.get('diagrams'):
                                self.borns.append(born)
                                self.fks_j_from_i[i.get('number')].append(\
                                                        j.get('number'))                                
                                
                                
    def find_borns_to_integrate(self, remove): #test written
        """finds double countings in the born configurations, sets the 
        is_to_integrate variable and if "remove" is True removes the 
        not needed ones from the born list"""
        #find the initial number of born configurations
        ninit = len(self.borns)
        
        for m in range(ninit):
            for n in range(m+1, ninit):
                born_m = self.borns[m]
                born_n = self.borns[n]
                # j j outgoing
                if born_m.j_fks > self.nincoming and \
                   born_n.j_fks > self.nincoming:
                    if (self.pdg_codes[born_m.i_fks - 1] == \
                        self.pdg_codes[born_n.i_fks - 1] and \
                        self.pdg_codes[born_m.j_fks - 1] == \
                        self.pdg_codes[born_n.j_fks - 1]) \
                        or \
                       (self.pdg_codes[born_m.i_fks - 1] == \
                        self.pdg_codes[born_n.j_fks - 1] and \
                        self.pdg_codes[born_m.j_fks - 1] == \
                        self.pdg_codes[born_n.i_fks - 1]):
                        
                        if born_m.i_fks > born_n.i_fks:
                            self.borns[m].is_to_integrate = False
                        elif born_m.i_fks == born_n.i_fks and \
                             born_m.j_fks > born_n.j_fks:
                            self.borns[n].is_to_integrate = False
                        else:
                            self.borns[m].is_to_integrate = False
                elif born_m.j_fks <= self.nincoming and \
                     born_n.j_fks == born_m.j_fks:
                    if self.pdg_codes[born_m.i_fks - 1] == \
                       self.pdg_codes[born_n.i_fks - 1]:
                        if born_m.i_fks > born_n.i_fks:
                            self.borns[n].is_to_integrate = False
                        else:
                            self.borns[m].is_to_integrate = False
        if remove:
            newborns = []
            for born in self.borns:
                if born.is_to_integrate:
                    newborns.append(born)
            self.borns = newborns
    
    
    def find_born_nbodyonly(self): #test written
        """finds the born configuration that includes the n-body contribution 
        in the virt0 and born0 running mode. By convention it is the  one 
        with soft singularities for which i/j are the largest"""
        imax = 0
        jmax = 0
        chosen = -1
        for n, born in enumerate (self.borns):
            if born.need_color_links:
                if born.i_fks >= imax and born.j_fks >= jmax and \
                        born.is_to_integrate:
                    chosen = n
                
        if chosen >=0:
            self.borns[chosen].is_nbody_only = True