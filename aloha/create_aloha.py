################################################################################
#
# Copyright (c) 2010 The MadGraph Development team and Contributors
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
from __future__ import division
import cmath
import copy
import cPickle
import glob
import logging
import numbers
import os
import re
import shutil
import sys
import time

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(root_path)
from aloha.aloha_object import *
import aloha
import aloha.aloha_writers as aloha_writers
import aloha.aloha_lib as aloha_lib
import aloha.aloha_object as aloha_object
import aloha.aloha_parsers as aloha_parsers
import aloha.aloha_fct as aloha_fct
try:
    import madgraph.iolibs.files as files
except Exception:
    import aloha.files as files
    
aloha_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger('ALOHA')

_conjugate_gap = 50
_spin2_mult = 1000


class ALOHAERROR(Exception): pass

class AbstractRoutine(object):
    """ store the result of the computation of Helicity Routine
    this is use for storing and passing to writer """
    
    def __init__(self, expr, outgoing, spins, name, infostr):
        """ store the information """

        self.spins = spins
        self.expr = expr
        self.name = name
        self.outgoing = outgoing
        self.infostr = infostr
        self.symmetries = []
        self.combined = []
        self.tag = []
        self.contracted = {}
        

        
    def add_symmetry(self, outgoing):
        """ add an outgoing """
        
        if not outgoing in self.symmetries:
            self.symmetries.append(outgoing)
    
    def add_combine(self, lor_list):
        """add a combine rule """
        
        if lor_list not in self.combined:
            self.combined.append(lor_list)
        
    def write(self, output_dir, language='Fortran', mode='self', combine=True,**opt):
        """ write the content of the object """
        writer = aloha_writers.WriterFactory(self, language, output_dir, self.tag)
        text = writer.write(mode=mode, **opt)
        if combine:
            for grouped in self.combined:
                if isinstance(text, tuple):
                    text = tuple([old.__add__(new)  for old, new in zip(text, 
                             writer.write_combined(grouped, mode=mode+'no_include', **opt))])
                else:
                    text += writer.write_combined(grouped, mode=mode+'no_include', **opt)        
        if aloha.mp_precision and 'MP' not in self.tag:
            self.tag.append('MP')
            text += self.write(output_dir, language, mode, **opt)
        return text

class AbstractRoutineBuilder(object):
    """ Launch the creation of the Helicity Routine"""
    
    prop_lib = {} # Store computation for the propagator
    counter = 0   # counter for statistic only
    
    class AbstractALOHAError(Exception):
        """ An error class for ALOHA"""
    
    def __init__(self, lorentz):
        """ initialize the run
        lorentz: the lorentz information analyzed (UFO format)
        language: define in which language we write the output
        modes: 0 for  all incoming particles 
              >0 defines the outgoing part (start to count at 1)
        """

        self.spins = self.spins = [abs(s) for s  in lorentz.spins]
        self.name = lorentz.name
        self.conjg = []
        self.tag = []
        self.outgoing = None
        self.lorentz_expr = lorentz.structure        
        self.routine_kernel = None
        self.spin2_massless = False
        self.contracted = {}
        self.fct = {}
        
    
    def compute_routine(self, mode, tag=[], factorize=True):
        """compute the expression and return it"""
        self.outgoing = mode
        self.tag = tag
        self.expr = self.compute_aloha_high_kernel(mode, factorize)
        return self.define_simple_output()
    
    def define_all_conjugate_builder(self, pair_list):
        """ return the full set of AbstractRoutineBuilder linked to fermion 
        clash"""
        
        solution = []

        for i, pair in enumerate(pair_list):
            new_builder = self.define_conjugate_builder(pair)
            solution.append(new_builder)
            solution += new_builder.define_all_conjugate_builder(pair_list[i+1:])
        return solution
                   
    def define_conjugate_builder(self, pairs=1):
        """ return a AbstractRoutineBuilder for the conjugate operation.
        If they are more than one pair of fermion. Then use pair to claim which 
        one is conjugated"""
        
        new_builder = copy.copy(self)
        new_builder.conjg = self.conjg[:]
        try:
            for index in pairs:
                new_builder.apply_conjugation(index) 
        except TypeError:
            new_builder.apply_conjugation(pairs) 
        return new_builder
    
    def apply_conjugation(self, pair=1):
        """ apply conjugation on self object"""
        
        nb_fermion = len([1 for s in self.spins if s % 2 == 0])        
        if (pair > 1 or nb_fermion >2) and not self.conjg:
            # self.conjg avoif multiple check
            data = aloha_fct.get_fermion_flow(self.lorentz_expr, nb_fermion)
            target = dict([(2*i+1,2*i+2) for i in range(nb_fermion//2)])
            if not data == target:
                text = """Unable to deal with 4(or more) point interactions
in presence of majorana particle/flow violation"""
                raise ALOHAERROR, text
        
        old_id = 2 * pair - 1
        new_id = _conjugate_gap + old_id
        
        self.kernel_tag = set()
        if not self.routine_kernel or isinstance(self.routine_kernel, str):
            self.routine_kernel = eval(self.lorentz_expr)
        
        # We need to compute C Gamma^T C^-1 = C_ab G_cb (-1) C_cd 
        #                  = C_ac G_bc (-1) C_bd = C_ac G_bc C_db
        self.routine_kernel = \
             C(new_id, old_id + 1) * self.routine_kernel * C(new_id + 1, old_id)
             
        self.lorentz_expr = '('+self.lorentz_expr+') * C(%s,%s) * C(%s,%s)' % \
                        (new_id, old_id + 1, new_id + 1, old_id ) 

        self.conjg.append(pair)

    
    def define_simple_output(self):
        """ define a simple output for this AbstractRoutine """
    
        infostr = str(self.lorentz_expr)        
        output = AbstractRoutine(self.expr, self.outgoing, self.spins, self.name, \
                                                                        infostr)
        output.contracted = dict([(name, aloha_lib.KERNEL.reduced_expr2[name])
                                          for name in aloha_lib.KERNEL.use_tag
                                          if name.startswith('TMP')])
        output.fct = dict([(name, aloha_lib.KERNEL.reduced_expr2[name])
                                          for name in aloha_lib.KERNEL.use_tag
                                          if name.startswith('FCT')])

        output.tag = [t for t in self.tag if not t.startswith('C')]
        output.tag += ['C%s' % pair for pair in self.conjg]
        return output

    def change_sign_for_outcoming_fermion(self):
        """change the sign of P for outcoming fermion in order to 
        correct the mismatch convention between HELAS and FR"""
        #flip_sign = []
        #for i in range(1,len(self.spins),2):
        #    if self.spins[i] == 2:
        #        flip_sign.append(str(i))
        
        #momentum_pattern = re.compile(r'\bP\(([\+\-\d]+),(%s)\)' % '|'.join(flip_sign))
        #lorentz_expr = momentum_pattern.sub(r'-P(\1,\2)', self.lorentz_expr)
        
        lorentz_expr = self.lorentz_expr
        calc = aloha_parsers.ALOHAExpressionParser()
        lorentz_expr = calc.parse(lorentz_expr)
        return lorentz_expr
                
    def compute_aloha_high_kernel(self, mode, factorize=True):
        """compute the abstract routine associate to this mode """

        # reset tag for particles
        aloha_lib.KERNEL.use_tag=set()
        #multiply by the wave functions
        nb_spinor = 0
        outgoing = self.outgoing
        if (outgoing + 1) // 2 in self.conjg:
            #flip the outgoing tag if in conjugate
            outgoing = outgoing + outgoing % 2 - (outgoing +1) % 2
        
        if not self.routine_kernel:
            AbstractRoutineBuilder.counter += 1
            logger.info('aloha creates %s routines' % self.name)
            try:
                lorentz = self.change_sign_for_outcoming_fermion()  
                self.routine_kernel = lorentz
                lorentz = eval(lorentz)
            except NameError as error:
                logger.error('unknow type in Lorentz Evaluation:%s'%str(error))
                raise ALOHAERROR, 'unknow type in Lorentz Evaluation: %s ' % str(error) 
            else:
                self.kernel_tag = set(aloha_lib.KERNEL.use_tag)
        elif isinstance(self.routine_kernel,str):
            lorentz = eval(self.routine_kernel)
            aloha_lib.KERNEL.use_tag = set(self.kernel_tag) 
        else:
            lorentz = copy.copy(self.routine_kernel)
            aloha_lib.KERNEL.use_tag = set(self.kernel_tag)
        for (i, spin ) in enumerate(self.spins):   
            id = i + 1
            #Check if this is the outgoing particle
            if id == outgoing:
                if spin == 1: 
                    lorentz *= complex(0,1)
                elif spin == 2:
                    # shift and flip the tag if we multiply by C matrices
                    if (id + 1) // 2 in self.conjg:
                        id += _conjugate_gap + id % 2 - (id +1) % 2
                    if (id % 2):
                        #propagator outcoming
                        lorentz *= SpinorPropagatorout(id, 'I2', outgoing)
                    else:
                    #    #propagator incoming
                        lorentz *= SpinorPropagatorin('I2', id, outgoing)
                elif spin == 3 :
                    lorentz *= VectorPropagator(id, 'I2', id)
                elif spin == 4:
                    # shift and flip the tag if we multiply by C matrices
                    if (id + 1) // 2 in self.conjg:
                        spin_id = id + _conjugate_gap + id % 2 - (id +1) % 2
                    else:
                        spin_id = id
                    nb_spinor += 1
                    if id %2:
                        lorentz *= Spin3halfPropagatorout(id, 'I2', spin_id,'I3', outgoing)
                    else:
                        lorentz *= Spin3halfPropagatorin('I2', id, 'I3', spin_id, outgoing)                      
                elif spin == 5 :
                    #lorentz *= 1 # delayed evaluation (fastenize the code)
                    if self.spin2_massless:
                        lorentz *= Spin2masslessPropagator(_spin2_mult + id, \
                                             2 * _spin2_mult + id,'I2','I3')
                    else:
                        lorentz *= Spin2Propagator(_spin2_mult + id, \
                                             2 * _spin2_mult + id,'I2','I3', id)
                else:
                    raise self.AbstractALOHAError(
                                'The spin value %s is not supported yet' % spin)
            else:
                # This is an incoming particle
                if spin == 1:
                    lorentz *= Scalar(id)
                elif spin == 2:
                    # shift the tag if we multiply by C matrices
                    if (id+1) // 2 in self.conjg:
                        spin_id = id + _conjugate_gap + id % 2 - (id +1) % 2
                    else:
                        spin_id = id
                    lorentz *= Spinor(spin_id, id)
                elif spin == 3:        
                    lorentz *= Vector(id, id)
                elif spin == 4:
                    # shift the tag if we multiply by C matrices
                    if (id+1) // 2 in self.conjg:
                        spin_id = id + _conjugate_gap + id % 2 - (id +1) % 2
                    else:
                        spin_id = id
                    nb_spinor += 1
                    lorentz *= Spin3Half(id, spin_id, id)
                elif spin == 5:
                    lorentz *= Spin2(1 * _spin2_mult + id, 2 * _spin2_mult + id, id)
                else:
                    raise self.AbstractALOHAError(
                                'The spin value %s is not supported yet' % spin)                    

        # If no particle OffShell
        if not outgoing:
            lorentz *= complex(0,-1)
            # Propagator are taken care separately
        
        lorentz = lorentz.simplify()
        
        # Modify the expression in case of loop-pozzorini
        if any((tag.startswith('L') for tag in self.tag if len(tag)>1)):
            return self.compute_loop_coefficient(lorentz, outgoing)
            
        lorentz = lorentz.expand()
        lorentz = lorentz.simplify()
        
        if factorize:
            lorentz = lorentz.factorize()
            
        lorentz.tag = set(aloha_lib.KERNEL.use_tag)
        return lorentz     
    
    def compute_loop_coefficient(self, lorentz, outgoing):
        

        l_in = [int(tag[1:]) for tag in self.tag if tag.startswith('L')][0]
        if (l_in + 1) // 2 in self.conjg:
            #flip the outgoing tag if in conjugate                                                                                                                                         
            l_in = l_in + l_in % 2 - (l_in +1) % 2                    
        assert l_in != outgoing, 'incoming Open Loop can not be the outcoming one'
        
        # modify the expression for the momenta
        # P_i -> P_i + P_L and P_o -> -P_o - P_L
        Pdep = [aloha_lib.KERNEL.get(P) for P in lorentz.get_all_var_names()
                                                      if P.startswith('_P')]

        Pdep = set([P for P in Pdep if P.particle in [outgoing, l_in]])
        for P in Pdep:
            if P.particle == l_in:
                sign = 1
            else:
                sign = -1
            id = P.id
            lorentz_ind = P.lorentz_ind[0]
            P_Lid = aloha_object.P(lorentz_ind, 'L')
            P_obj = aloha_object.P(lorentz_ind, P.particle)
            new_expr = sign*(P_Lid + P_obj)
            lorentz = lorentz.replace(id, new_expr)

        # Compute the variable from which we need to split the expression
        var_veto =  ['PL_0', 'PL_1', 'PL_2', 'PL_3']
        spin = aloha_writers.WriteALOHA.type_to_variable[self.spins[l_in-1]]
        size = aloha_writers.WriteALOHA.type_to_size[spin]-1
        var_veto += ['%s%s_%s' % (spin,l_in,i) for i in range(1,size)]
        # compute their unique identifiant
        veto_ids = aloha_lib.KERNEL.get_ids(var_veto)
        
        lorentz = lorentz.expand(veto = veto_ids)
        lorentz = lorentz.simplify()
        coeff_expr = lorentz.split(veto_ids)
        
        for key, expr in coeff_expr.items():
            expr = expr.simplify()
            coeff_expr[key] = expr.factorize()
        coeff_expr.tag = set(aloha_lib.KERNEL.use_tag)

        return coeff_expr
                        
    def define_lorentz_expr(self, lorentz_expr):
        """Define the expression"""
        
        self.expr = lorentz_expr
    
    def define_routine_kernel(self, lorentz=None):
        """Define the kernel at low level"""
        
        if not lorentz:
            logger.info('compute kernel %s' % self.counter)
            AbstractRoutineBuilder.counter += 1  
            lorentz = eval(self.lorentz_expr)
                 
            if isinstance(lorentz, numbers.Number):
                self.routine_kernel = lorentz
                return lorentz
            lorentz = lorentz.simplify()
            lorentz = lorentz.expand()
            lorentz = lorentz.simplify()        
        
        self.routine_kernel = lorentz
        return lorentz

    
    @staticmethod
    def get_routine_name(name, outgoing):
        """return the name of the """
        
        name = '%s_%s' % (name, outgoing) 
        return name
            
    @classmethod
    def load_library(cls, tag):
        # load the library
        if tag in cls.prop_lib:
            return
        else:
            cls.prop_lib = create_prop_library(tag, cls.aloha_lib)
        

class CombineRoutineBuilder(AbstractRoutineBuilder):
    """A special builder for combine routine if needed to write those
        explicitely.
    """
    def __init__(self, l_lorentz):
        """ initialize the run
        l_lorentz: list  of lorentz information analyzed (UFO format)
        language: define in which language we write the output
        modes: 0 for  all incoming particles 
              >0 defines the outgoing part (start to count at 1)
        """

        lorentz = l_lorentz[0]
        self.spins = lorentz.spins
        l_name = [l.name for l in l_lorentz]
        self.name = aloha_writers.combine_name(l_name[0], l_name[1:], None)
        self.conjg = []
        self.tag = []
        self.outgoing = None
        self.lorentz_expr = []
        for i, lor in enumerate(l_lorentz):
            self.lorentz_expr.append( 'Coup(%s) * (%s)' % (i+1, lor.structure))
        self.lorentz_expr = ' + '.join(self.lorentz_expr)
        self.routine_kernel = None
        self.spin2_massless = False
        self.contracted = {}
        self.fct = {}

class AbstractALOHAModel(dict):
    """ A class to build and store the full set of Abstract ALOHA Routine"""


    def __init__(self, model_name, write_dir=None, format='Fortran', 
                 explicit_combine=False):
        """ load the UFO model and init the dictionary """
        
        aloha_lib.KERNEL.clean()
        # Option
        self.explicit_combine = explicit_combine
        
        # Extract the model name if combined with restriction
        model_name_pattern = re.compile("^(?P<name>.+)-(?P<rest>[\w\d_]+)$")
        model_name_re = model_name_pattern.match(model_name)
        if model_name_re:
            name = model_name_re.group('name')
            rest = model_name_re.group("rest")
            if rest == 'full' or \
               os.path.isfile(os.path.join(root_path, "models", name,
                                           "restrict_%s.dat" % rest)):
                model_name = model_name_re.group("name")

        # load the UFO model
        try:
            python_pos = model_name 
            __import__(python_pos)
        except Exception:
            python_pos = 'models.%s' % model_name 
            __import__(python_pos)
        self.model = sys.modules[python_pos]
        # find the position on the disk
        self.model_pos = os.path.dirname(self.model.__file__)

        # list the external routine
        self.external_routines = [] 

        # init the dictionary
        dict.__init__(self)
        self.symmetries = {}
        self.multiple_lor = {}
        
        # check the mass of spin2 (if any)
        self.massless_spin2 = self.has_massless_spin2()
        
        if write_dir:
            self.main(write_dir,format=format)
            
    def main(self, output_dir, format='Fortran'):
        """ Compute if not already compute. 
            Write file in models/MY_MODEL/MY_FORMAT.
            copy the file to output_dir
        """
        ext = {'Fortran':'f','Python':'py','CPP':'h'}
        
        
        # Check if a pickle file exists
        if not self.load():
            self.compute_all()
        logger.info(' %s aloha routine' % len(self))
            
        # Check that output directory exists
        if not output_dir:
            output_dir = os.path.join(self.model_pos, format.lower())
            logger.debug('aloha output dir is %s' % output_dir) 
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        # Check that all routine are generated at default places:
        for (name, outgoing), abstract in self.items():
            routine_name = AbstractRoutineBuilder.get_routine_name(name, outgoing)
            if not glob.glob(os.path.join(output_dir, routine_name) + '.' + ext[format]):
                abstract.write(output_dir, format) 
            else:
                logger.info('File for %s already present, skip the writing of this file' % routine_name)
                   
        
    def save(self, filepos=None):
        """ save the current model in a pkl file """
        
        logger.info('save the aloha abstract routine in a pickle file')
        if not filepos:
            filepos = os.path.join(self.model_pos,'aloha.pkl') 
        
        fsock = open(filepos, 'w')
        cPickle.dump(dict(self), fsock)
        
    def load(self, filepos=None):
        """ reload the pickle file """
        return False
        if not filepos:
            filepos = os.path.join(self.model_pos,'aloha.pkl') 
        if os.path.exists(filepos):
            fsock = open(filepos, 'r')
            self.update(cPickle.load(fsock))        
            return True
        else:
            return False
        
    def get(self, lorentzname, outgoing):
        """ return the AbstractRoutine with a given lorentz name, and for a given
        outgoing particle """
        
        try:
            return self[(lorentzname, outgoing)]
        except Exception:
            logger.warning('(%s, %s) is not a valid key' % 
                                                       (lorentzname, outgoing) )
            return None
    
    def set(self, lorentzname, outgoing, abstract_routine):
        """ add in the dictionary """
    
        self[(lorentzname, outgoing)] = abstract_routine
    
    def compute_all(self, save=True, wanted_lorentz = []):
        """ define all the AbstractRoutine linked to a model """


        # Search identical particles in the vertices in order to avoid
        #to compute identical contribution
        self.look_for_symmetries()
        conjugate_list = self.look_for_conjugate()
        self.look_for_multiple_lorentz_interactions()
        
        if not wanted_lorentz:
            wanted_lorentz = [l.name for l in self.model.all_lorentz]
        for lorentz in self.model.all_lorentz:
            if not lorentz.name in wanted_lorentz:
                # Only include the routines we ask for
                continue
            
            if -1 in lorentz.spins:
                # No Ghost in ALOHA
                continue 
            
            if lorentz.structure == 'external':
                for i in range(len(lorentz.spins)):
                    self.external_routines.append('%s_%s' % (lorentz.name, i))
                continue
            
            builder = AbstractRoutineBuilder(lorentz)
            # add information for spin2mass
            if 5 in lorentz.spins and self.massless_spin2 is not None:
                builder.spin2_massless = self.massless_spin2
            self.compute_aloha(builder)
            
            if lorentz.name in self.multiple_lor:
                for m in self.multiple_lor[lorentz.name]:
                    for outgoing in range(len(lorentz.spins)+1):
                        try:
                            self[(lorentz.name, outgoing)].add_combine(m)
                        except Exception:
                            pass # this routine is a symmetric one, so it 
                                 # already has the combination.
                    
            if lorentz.name in conjugate_list:
                conjg_builder_list= builder.define_all_conjugate_builder(\
                                                   conjugate_list[lorentz.name])
                for conjg_builder in conjg_builder_list:
                    # No duplication of conjugation:
                    assert conjg_builder_list.count(conjg_builder) == 1
                    self.compute_aloha(conjg_builder, lorentz.name)
                    if lorentz.name in self.multiple_lor:
                        for m in self.multiple_lor[lorentz.name]:
                            for outgoing in range(len(lorentz.spins)+1):
                                realname = conjg_builder.name + ''.join(['C%s' % pair for pair in conjg_builder.conjg])
                                try:
                                    self[(realname, outgoing)].add_combine(m)
                                except Exception,error:
                                    self[(realname, self.symmetries[lorentz.name][outgoing])].add_combine(m)          
                       
        if save:
            self.save()
    
    def add_Lorentz_object(self, lorentzlist):
        """add a series of Lorentz structure created dynamically"""
        
        for lor in lorentzlist:
            if not hasattr(self.model.lorentz, lor.name):
                setattr(self.model.lorentz, lor.name, lor)
    
    # Notice that when removing the quick fix, one should also remove the two
    # additional arguments byPassFix and forceLoop
    # To emulate the behavior without the quick fix, simply replace the default
    # value of byPassFix by True.
    # == START AD-HOC QUICK FIX ==
    def compute_subset(self, data, byPassFix=False, forceLoop=False):
    # == END AD-HOC QUICK FIX ==
        """ create the requested ALOHA routine. 
        data should be a list of tuple (lorentz, tag, outgoing)
        tag should be the list of special tag (like conjugation on pair)
        to apply on the object """

        # == START AD-HOC QUICK FIX WARNING ==        

        # In order to be able to use open loops in unitary gauge, one must make
        # sure ALOHA does not include the longitudinal part of the loop
        # propagator for the gluon.
        # In further versions, this will be insured by having separate routines
        # for the massive and massless vector propagators. For now, I use a 
        # quick fix which consists in calling twice aloha (rather compute_subset)
        # with only the loop routines (then in feynman gauge not matter what)
        # and a second time with the rest as usual.

        if not byPassFix: 
            
            data_loop = [d for d in data if any((t.startswith('L') for t in d[1]))]
            data_tree = [d for d in data if not any((t.startswith('L') for t in d[1]))]
            
            if data_loop == []:
                self.compute_subset(data,byPassFix=True)
                return
            
            self.compute_subset(data_tree, byPassFix=True, forceLoop=True)
            old_aloha_gauge = aloha.unitary_gauge
            aloha.unitary_gauge = False
            self.compute_subset(data_loop, byPassFix=True, forceLoop=True)
            aloha.unitary_gauge = old_aloha_gauge
            return
        # == END AD-HOC QUICK FIX ==

        # Search identical particles in the vertices in order to avoid
        #to compute identical contribution
        self.look_for_symmetries()
        # reorganize the data (in order to use optimization for a given lorentz
        #structure
	    # HSS,27/11/2012
        aloha.loop_mode = False
	    # self.explicit_combine = False
	    # HSS
        request = {}
        for list_l_name, tag, outgoing in data:
            #allow tag to have integer for retro-compatibility
            conjugate = [i for i in tag if isinstance(i, int)]
            tag =  [i for i in tag if isinstance(i, str)]
            tag = tag + ['C%s'%i for i in conjugate] 

            conjugate = tuple([int(c[1:]) for c in tag if c.startswith('C')])
            loop = any((t.startswith('L') for t in tag))
            # When removing the QUICK FIX, also remove 'forceLoop'
            # == START AD-HOC QUICK FIX ==
            if loop or forceLoop:
            # == END AD-HOC QUICK FIX ==
                aloha.loop_mode = True
                self.explicit_combine = True
                

                
                       
            for l_name in list_l_name:
                try:
                    request[l_name][conjugate].append((outgoing,tag))
                except Exception:
                    try:
                        request[l_name][conjugate] = [(outgoing,tag)]
                    except Exception:
                        request[l_name] = {conjugate: [(outgoing,tag)]}
                           
        # Loop on the structure to build exactly what is request
        for l_name in request:
            lorentz = eval('self.model.lorentz.%s' % l_name)
            if lorentz.structure == 'external':
                for tmp in request[l_name]:
                    for outgoing, tag in request[l_name][tmp]:
                        name = aloha_writers.get_routine_name(lorentz.name,outgoing=outgoing,tag=tag)
                        if name not in self.external_routines:
                            print name
                            self.external_routines.append(name)
                continue
            
            builder = AbstractRoutineBuilder(lorentz)
            # add information for spin2mass
            if 5 in lorentz.spins and self.massless_spin2 is not None:
                builder.spin2_massless = self.massless_spin2 
            
            for conjg in request[l_name]:
                #ensure that routines are in rising order (for symetries)
                def sorting(a,b):
                    if a[0] < b[0]: return -1
                    else: return 1
                routines = request[l_name][conjg]
                routines.sort(sorting)
                if not conjg:
                    # No need to conjugate -> compute directly
                    self.compute_aloha(builder, routines=routines)
                else:
                    # Define the high level conjugate routine
                    conjg_builder = builder.define_conjugate_builder(conjg)
                    # Compute routines
                    self.compute_aloha(conjg_builder, symmetry=lorentz.name,
                                         routines=routines)
            
        
        # Build mutiple lorentz call
        for list_l_name, tag, outgoing in data:
            if len(list_l_name) ==1:
                continue
            #allow tag to have integer for retrocompatibility
            conjugate = [i for i in tag if isinstance(i, int)]
            tag =  [i for i in tag if isinstance(i, str)]
            tag = tag + ['C%s'%i for i in conjugate] 
            
            if not self.explicit_combine:
                lorentzname = list_l_name[0]
                lorentzname += ''.join(tag)
                if self.has_key((lorentzname, outgoing)):
                    self[(lorentzname, outgoing)].add_combine(list_l_name[1:])
                else:
                    lorentz = eval('self.model.lorentz.%s' % lorentzname)
                    assert lorentz.structure == 'external'
            else:
                l_lorentz = []
                for l_name in list_l_name: 
                    l_lorentz.append(eval('self.model.lorentz.%s' % l_name))
                builder = CombineRoutineBuilder(l_lorentz)
                # add information for spin2mass
                if 5 in l_lorentz[0].spins and self.massless_spin2 is not None:
                    builder.spin2_massless = self.massless_spin2 
            
                for conjg in request[list_l_name[0]]:
                    #ensure that routines are in rising order (for symetries)
                    def sorting(a,b):
                        if a[0] < b[0]: return -1
                        else: return 1
                    routines = request[list_l_name[0]][conjg]
                    routines.sort(sorting)
                    if not conjg:
                        # No need to conjugate -> compute directly
                        self.compute_aloha(builder, routines=routines)
                    else:
                        # Define the high level conjugate routine
                        conjg_builder = builder.define_conjugate_builder(conjg)
                        # Compute routines
                        self.compute_aloha(conjg_builder, symmetry=lorentz.name,
                                        routines=routines)
                      
  
                            
    def compute_aloha(self, builder, symmetry=None, routines=None, tag=[]):
        """ define all the AbstractRoutine linked to a given lorentz structure
        symmetry authorizes to use the symmetry of anoter lorentz structure.
        routines to define only a subset of the routines."""

        name = builder.name
        if not symmetry:
            symmetry = name
        if not routines:
            tag = ['C%s' % i for i in builder.conjg]
            routines = [ tuple([i,tag]) for i in range(len(builder.spins) + 1 )]

        # Create the routines
        for outgoing, tag in routines:
            symmetric = self.has_symmetries(symmetry, outgoing, valid_output=routines)
            realname = name + ''.join(tag)
            if (realname, outgoing) in self:
                continue # already computed
            
            if symmetric:
                self.get(realname, symmetric).add_symmetry(outgoing)
            else:
                wavefunction = builder.compute_routine(outgoing, tag)
                #Store the information
                self.set(realname, outgoing, wavefunction)
          

    def compute_aloha_without_kernel(self, builder, symmetry=None, routines=None):
        """define all the AbstractRoutine linked to a given lorentz structure
        symmetry authorizes to use the symmetry of anoter lorentz structure.
        routines to define only a subset of the routines. 
        Compare to compute_aloha, each routines are computed independently.
        """

        name = builder.name
        if not routines:
            routines = [ tuple([i,[]]) for i in range(len(builder.spins) + 1 )]         
        
        for outgoing, tag in routines:
            builder.routine_kernel = None
            wavefunction = builder.compute_routine(outgoing, tag)
            self.set(name, outgoing, wavefunction)


    def write(self, output_dir, language):
        """ write the full set of Helicity Routine in output_dir"""
        for abstract_routine in self.values():
            abstract_routine.write(output_dir, language)

        for routine in self.external_routines:
            self.locate_external(routine, language, output_dir)
        
        #self.write_aloha_file_inc(output_dir)
    
    def locate_external(self, name, language, output_dir=None):
        """search a valid external file and copy it to output_dir directory"""
        
        language_to_ext = {'Python': 'py',
                           'Fortran' : 'f',
                           'CPP': 'C'}
        ext = language_to_ext[language]
        paths = [os.path.join(self.model_pos, language), self.model_pos, 
                           os.path.join(root_path, 'aloha', 'template_files', )]

        ext_files  = []
        for path in paths:
            ext_files = glob.glob(os.path.join(path, '%s.%s' % (name, ext)))
            if ext_files:
                break
        else: 

            raise ALOHAERROR, 'No external routine \"%s.%s\" in directories\n %s' % \
                        (name, ext, '\n'.join(paths))
       
        if output_dir:
            for filepath in ext_files:
                
                files.cp(filepath, output_dir)
        return ext_files
                    
        

    def look_for_symmetries(self):
        """Search some symmetries in the vertices.
        We search if some identical particles are in a vertices in order
        to avoid to compute symmetrical contributions"""
        
        for vertex in self.model.all_vertices:
            for i, part1 in enumerate(vertex.particles):
                for j in range(i-1,-1,-1):
                    part2 = vertex.particles[j]
                    if part1.pdg_code == part2.pdg_code and part1.color == 1:
                        if part1.spin == 2 and (i % 2 != j % 2 ):
                            continue 
                        for lorentz in vertex.lorentz:
                            if self.symmetries.has_key(lorentz.name):
                                if self.symmetries[lorentz.name].has_key(i+1):
                                    self.symmetries[lorentz.name][i+1] = max(self.symmetries[lorentz.name][i+1], j+1)
                                else:
                                    self.symmetries[lorentz.name][i+1] = j+1
                            else:
                                self.symmetries[lorentz.name] = {i+1:j+1}
                        break
    
    def look_for_multiple_lorentz_interactions(self):
        """Search the interaction associate with more than one lorentz structure.
        If those lorentz structure have the same order and the same color then
        associate a multiple lorentz routines to ALOHA """
        
        orders = {}
        for coup in self.model.all_couplings:
            orders[coup.name] = str(coup.order)
        
        for vertex in self.model.all_vertices:
            if len(vertex.lorentz) == 1:
                continue
            #remove ghost
            if -1 in vertex.lorentz[0].spins:
                continue
            
            # assign each order/color to a set of lorentz routine
            combine = {}
            for (id_col, id_lor), coups in vertex.couplings.items():
                if not isinstance(coups, list):
                    coups = [coups]
                for coup in coups:
                    order = orders[coup.name]
                    key = (id_col, order)
                    if key in combine:
                        combine[key].append(id_lor)
                    else:
                        combine[key] = [id_lor]
                    
            # Check if more than one routine are associated
            for list_lor in combine.values():
                if len(list_lor) == 1:
                    continue
                list_lor.sort() 
                main = vertex.lorentz[list_lor[0]].name 
                if main not in self.multiple_lor:
                    self.multiple_lor[main] = []
                
                info = tuple([vertex.lorentz[id].name for id in list_lor[1:]])
                if info not in self.multiple_lor[main]:
                    self.multiple_lor[main].append(info)
                
    def has_massless_spin2(self):
        """Search if the spin2 particles are massless or not"""
        
        massless = None
        for particle in self.model.all_particles:
            if particle.spin == 5:
                if massless is None:
                    massless = (particle.mass == 'Zero')
                elif massless != (particle.mass == 'Zero'):
                    raise ALOHAERROR, 'All spin 2 should be massive or massless'
        return massless     
                    
    def has_symmetries(self, l_name, outgoing, out=None, valid_output=None):
        """ This returns out if no symmetries are available, otherwise it finds 
        the lowest equivalent outgoing by recursivally calling this function.
        auth is a list of authorize output, if define"""

        try:
            equiv = self.symmetries[l_name][outgoing]
        except Exception:
            return out
        else:
            if not valid_output or equiv in valid_output:
                return self.has_symmetries(l_name, equiv, out=equiv, 
                                                      valid_output=valid_output)
            else:
                return self.has_symmetries(l_name, equiv, out=out,              
                                                      valid_output=valid_output)
        
    def look_for_conjugate(self):
        """ create a list for the routine needing to be conjugate """

        # Check if they are majorana in the model.
        need = False
        for particle in self.model.all_particles:
            if particle.spin == 2 and particle.selfconjugate:
                need = True
                break

        if not need:
            for interaction in self.model.all_vertices:
                fermions = [p for p in interaction.particles if p.spin == 2]
                for i in range(0, len(fermions), 2):
                    if fermions[i].pdg_code * fermions[i+1].pdg_code > 0:
                        # This is a fermion flow violating interaction
                        need = True
                        break

        # No majorana particles    
        if not need:
            return {}
        
        conjugate_request = {}
        # Check each vertex if they are fermion and/or majorana
        for vertex in self.model.all_vertices:
            for i in range(0, len(vertex.particles), 2):
                part1 = vertex.particles[i]
                if part1.spin !=2:
                    # deal only with fermion
                    break
                # check if this pair contains a majorana
                if part1.selfconjugate:
                    continue
                part2 = vertex.particles[i + 1]
                if part2.selfconjugate:
                    continue
                
                # No majorana => add the associate lorentz structure
                for lorentz in vertex.lorentz:
                    try:
                        conjugate_request[lorentz.name].add(i//2+1)
                    except Exception:
                        conjugate_request[lorentz.name] = set([i//2+1])
        
        for elem in conjugate_request:
            conjugate_request[elem] = list(conjugate_request[elem])
        
        return conjugate_request
            
        
            
def write_aloha_file_inc(aloha_dir,file_ext, comp_ext):
    """find the list of Helicity routine in the directory and create a list 
    of those files (but with compile extension)"""

    aloha_files = []
    
    # Identify the valid files
    alohafile_pattern = re.compile(r'''_\d%s''' % file_ext)
    for filename in os.listdir(aloha_dir):
        if os.path.isfile(os.path.join(aloha_dir, filename)):
            if alohafile_pattern.search(filename):
                aloha_files.append(filename.replace(file_ext, comp_ext))

    text="ALOHARoutine = "
    text += ' '.join(aloha_files)
    text +='\n'
    file(os.path.join(aloha_dir, 'aloha_file.inc'), 'w').write(text) 


            
def create_prop_library(tag, lib={}):
    
    def create(obj):
        """ """
        obj= obj.simplify()
        obj = obj.expand()
        obj = obj.simplify()
        return obj        
    
    # avoid to add tag in global
    old_tag = set(aloha_lib.KERNEL.use_tag)
    print 'create lib',tag
    name, i = tag
    if name == "Spin2Prop":
        lib[('Spin2Prop',i)] = create( Spin2Propagator(_spin2_mult + i, \
                                             2 * _spin2_mult + i,'I2','I3', i) )
    elif name == "Spin2PropMassless":
        lib[('Spin2PropMassless',i)] = create( Spin2masslessPropagator(
                             _spin2_mult + i, 2 * _spin2_mult + i,'I2','I3'))
    
    aloha_lib.KERNEL.use_tag = old_tag
    return lib


if '__main__' == __name__:       
    logging.basicConfig(level=0)
    #create_library()
    import profile       
    #model 
      
    start = time.time()
    def main():
        alohagenerator = AbstractALOHAModel('sm') 
        alohagenerator.compute_all(save=False)
        return alohagenerator
    def write(alohagenerator):
        alohagenerator.write('/tmp/', 'Python')
    alohagenerator = main()
    logger.info('done in %s s' % (time.time()-start))
    write(alohagenerator)
    #profile.run('main()')
    #profile.run('write(alohagenerator)')
    stop = time.time()
    logger.info('done in %s s' % (stop-start))
  







