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
import os
import logging
import sys
import cPickle
import numbers
import time


from helas.helasamp_object import *
import helas.WriteHelas as WriteHelas


helas_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger('HelasGenerator')

class AbstractHelas(object):
    """ Launch the creation of the Helas Routine"""
    
    spin_to_tag = {1:'S', 2:'F', 3:'V', 5:'T'}
    helas_lib = None
    counter = 0
    
    class AbstractHelasError(Exception):
        """ An error class for the AbstractHelas class"""
    
    def __init__(self, lorentz, mode, auto_run=True, **opt):
        """ initialize the run
        lorentz: the lorentz information analyzed (UFO format)
        language: define in which language we write the output
        modes: 0 for  all incoming particles 
               >0 defines the outgoing part (start to count at 1)
        """

        self.spins = lorentz.spins
        self.name = lorentz.name
        self.outgoing = mode
        self.lorentz_expr = lorentz.structure        
        self.helas_kernel = None
        
        if auto_run:
            self.expr = self.compute_helas(mode)
    
    def compute_helas(self, mode):
        """ """
        self.compute_helas_high_kernel(mode)
    
    def compute_helas_high_kernel(self, mode):
        """compute the abstract routine associate to this mode """
        
        #multiply by the wave functions
        nb_spinor = 0
        if not self.helas_kernel:
            AbstractHelas.counter +=1
            logger.info('new kernel %s' % self.counter)
            try:       
                lorentz = eval(self.lorentz_expr)
            except NameError, why:
                print 'unknow type in Lorentz Evaluation'
                raise
            else:
                self.helas_kernel = lorentz
        else:
            lorentz = self.helas_kernel
            
        for (i, spin ) in enumerate(self.spins):
            id = i + 1
            
            #Check if this is the outgoing particle
            if id == self.outgoing:
                if spin == 1: 
                    lorentz *= complex(0,1)
                elif spin == 2:
                    nb_spinor += 1
                    if nb_spinor %2:
                        lorentz *= SpinorPropagator(i, 'I2', i)
                    else:
                        lorentz *= SpinorPropagator('I2', i, i) 
                elif spin == 3 :
                    lorentz *= VectorPropagator(i, 'I2', i)
                elif spin == 5 :
                    lorentz *= 1 # delayed evaluation (fastenize the code)
                else:
                    raise self.AbstractHelasError(
                                'The spin value %s is not supported yet' % spin)
            else:
                # This is an incoming particle
                if spin == 1:
                    lorentz *= Scalar(id)
                elif spin == 2:
                    nb_spinor += 1
                    lorentz *= Spinor(id, id)
                elif spin == 3:        
                    lorentz *= Vector(id, id)
                elif spin == 5:
                    lorentz *= Spin2(10*id+1, 10*id+2, 'I2', 'I3', id)
                else:
                    raise self.AbstractHelasError(
                                'The spin value %s is not supported yet' % spin)                    

        # If no particle OffShell
        if self.outgoing:
            lorentz /= DenominatorPropagator(self.outgoing)
            lorentz.tag.add('OM%s' % self.outgoing )  
            lorentz.tag.add('P%s' % self.outgoing)  
            lorentz.tag.add('W%s' % self.outgoing)  
        else:
            lorentz *= complex(0,-1)


            
        lorentz = lorentz.simplify()
        lorentz = lorentz.expand()
        lorentz = lorentz.simplify()
        if self.spins[self.outgoing-1] == 5:
            if not self.helas_lib:
                AbstractHelas.load_library()
            lorentz *= self.helas_lib[('Spin2Prop', id)]
        lorentz = lorentz.factorize()

        return lorentz         
        
    def compute_helas_low_kernel(self, mode):
        """ compute the abstract routine associate to this mode """
        
        if not self.helas_kernel:
            self.helas_kernel = self.define_helas_kernel()
        
        if not AbstractHelas.helas_lib:
            AbstractHelas.load_library()
        
        #multiply by the wave functions
        nb_spinor = 0
        lorentz = self.helas_kernel
        for (i, spin ) in enumerate(self.spins):
            id = i + 1
            
            #Check if this is the outgoing particle
            if id == self.outgoing:
                if spin == 1 : 
                    lorentz *= self.helas_lib[('ScalarProp', id)]
                elif spin == 2 :
                    nb_spinor += 1
                    lorentz *= self.helas_lib[('SpinorProp', id, nb_spinor %2)]
                elif spin == 3 :
                    lorentz *= self.helas_lib[('VectorProp', id)]
                elif spin == 5 :
                    lorentz *= self.helas_lib[('Spin2Prop', id)]
                else:
                    raise self.AbstractHelasError(
                                'The spin value %s is not supported yet' % spin)
            else:
                # This is an incoming particle
                if spin == 1:
                    lorentz *= self.helas_lib[('Scalar', id)]
                elif spin == 2:
                    nb_spinor += 1
                    lorentz *= self.helas_lib[('Spinor', id)]
                elif spin == 3:        
                    lorentz *= self.helas_lib[('Vector', id)]
                elif spin == 5:
                    lorentz *= self.helas_lib[('Spin2', id)]
                else:
                    raise self.AbstractHelasError(
                                'The spin value %s is not supported yet' % spin)                    
            
        # If no particle OffShell
        if self.outgoing:
            lorentz /= self.helas_lib[('Denom', id)]
            lorentz.tag.add('OM%s' % self.outgoing )  
            lorentz.tag.add('P%s' % self.outgoing)  
            lorentz.tag.add('W%s' % self.outgoing)  
        else:
            lorentz *= complex(0,-1)
                        
        lorentz = lorentz.simplify()
        lorentz = lorentz.factorize()
        return lorentz 
                
    def define_helas(self, lorentz_expr):
        """Define the expression"""
        
        self.lorentz_expr = lorentz_expr
    
    def define_helas_kernel(self, lorentz=None):
        """Define the kernel at low level"""
        
        if not lorentz:
            logger.info('compute kernel %s' % self.counter)
            AbstractHelas.counter += 1  
            try:        
                lorentz = eval(self.lorentz_expr)
            except NameError, why:
                print dir(why)
                print why.args
                raise
                
            if isinstance(lorentz, numbers.Number):
                self.helas_kernel = lorentz
                return lorentz
            lorentz = lorentz.simplify()
            lorentz = lorentz.expand()
            lorentz = lorentz.simplify()        
        
        self.helas_kernel = lorentz
        return lorentz

        
 
   
    def gethelasname(self):
        """return the name of the """
        
        helasname = self.name + '_'
        for i in len(self.spins):
            if i + 1 == self.outgoing:
                helasname += '0'
            else:
                helasname += '1'
        return helasname
        

    def write(self, output_dir, language):
        """write the lorentz structure in the Helas Routine file filepos"""

        name = os.path.join(output_dir, self.gethelasname())
        infostr = str(self.expr)
        isoutgoing = lambda i: i + 1 == self.outgoing
        info = [( self.spin_to_tag[spin], isoutgoing(i) ) \
                                           for i, spin in enumerate(self.spins)]

        getattr(WriteHelas, 'HelasWriterFor%s' % language)\
                            (self.expr, info, name, infostr).write()
    
    @classmethod
    def load_library(cls):
    # load the library
        fsock = open(os.path.join(helas_path, 'HelasLib.pkl'), 'r')
        cls.helas_lib = cPickle.load(fsock)
        

class AbstractHelasModel(dict):
    """ A class to buid and store the full set of Abstract Helas Routine"""

    def __init__(self, model_name):
        """ load the UFO model and init the dictionary """
        
        # load the UFO model
        python_pos = 'models.%s' % model_name 
        __import__(python_pos)
        self.model = sys.modules[python_pos]
        
        # find the position on the disk
        self.model_pos = os.path.dirname(self.model.__file__)

        #init the dictionary
        dict.__init__(self)
        self.symmetries = {}
        
    def save(self, filepos=None):
        """ save the current model in a pkl file """
        
        ff = open(os.path.join(self.model_pos,'helas.pkl'), 'w')
        cPickle.dump(dict(self), ff)
        
    def get(self, lorentzname, outgoing):
        """ return the AbstractHelas with a given lorentz name, and for a given
        outgoing particle """
        
        try:
            return self[(lorentzname, outgoing)]
        except:
            print self.keys()
            logger.warning('(%s, %s) is not a valid key' % 
                                                       (lorentzname, outgoing) )
            return None
    
    def set(self, lorentzname, outgoing, abstracthelas):
        """ add in the dictionary """
    
        self[(lorentzname, outgoing)] = abstracthelas
    
    def compute_all(self):
        """ define all the AbstractHelas linked to a model """

        # Search identical particles in the vertices in order to avoid
        #to compute identical contribution
        #self.look_for_symmetries()
        
        for lorentz in self.model.all_lorentz:
            self.compute_for_lorentz(lorentz)
        
        self.save()
        
    def compute_for_lorentz(self, lorentz):
        """ """
        self.compute_lorentz_with_kernel(lorentz)
        
    def compute_lorentz_whithout_kernel(self, lorentz):
        """define all the AbstractHelas"""
        
        name = lorentz.name
        for outgoing in range(len(lorentz.spins) + 1 ):
            wavefunctions = AbstractHelas(lorentz, outgoing)
            self.set(name, outgoing, wavefunctions)
        
        
    def compute_lorentz_with_kernel(self, lorentz):
        """ define all the AbstractHelas linked to a given lorentz structure"""
        
        name = lorentz.name
        # first compute the amplitude contribution
        wavefunctions = AbstractHelas(lorentz, 0)
        self.set(name, 0, wavefunctions)
        
        # Take the kernel for future operation
        kernel = wavefunctions.helas_kernel
            
        # Create the routine associate to an externel particles
        for outgoing in range(1, len(lorentz.spins) + 1 ):
            # Make the computation
            wavefunctions = AbstractHelas(lorentz, outgoing, auto_run=False)
            if self.symmetries.has_key(lorentz.name) and \
                                   outgoing -1 in self.symmetries[lorentz.name]:
                equiv_out = self.symmetries[lorentz.name][outgoing-1]
                lorentz_expr = self.get(lorentz.name, equiv_out).lorentz_expr
                wavefunctions.define_helas( lorentz_expr )
            else:
                wavefunctions.define_helas_kernel(kernel)
                wavefunctions.compute_helas(outgoing)
            
            #Store the information
            self.set(name, outgoing, wavefunctions)

    def write(self, output_dir, language):
        """ write the full set of Helas Routine in output_dir"""
        
        for abstract_helas in self.values():
            abstract_helas.write(output_dir, language)
        

    def look_for_symmetries(self):
        """Search some symmetries in the vertices.
        We search if some identical particles are in a vertices in order
        to avoid to compute symmetrical contributions"""
        
        for vertex in self.model.all_vertices:
            for i, part1 in enumerate(vertex.particles):
                for j in range(i):
                    part2 = vertex.particles[j]
                    if part1.name == part2.name:
                        for lorentz in vertex.lorentz:
                            if self.symmetries.has_key(lorentz.name):
                                self.symmetries[lorentz.name][i] = j
                            else:
                                self.symmetries[lorentz.name] = {i:j}
                        break

def create_library():
    
    def create(obj):
        """ """
        obj= obj.simplify()
        obj = obj.expand()
        obj = obj.simplify()
        return obj        
    
    
    lib = {} # key: (name, part_nb, special) -> object
    for i in range(1, 10):
        logger.info('step %s/9' % i)
        lib[('Scalar', i)] = create( Scalar(i) )
        lib[('ScalarProp', i)] = complex(0,1)
        lib[('Denom', i )] = create( DenominatorPropagator(i) )
        lib[('Spinor', i )] = create( Spinor(i, i) )
        lib[('SpinorProp', i, 0)] = create( SpinorPropagator(i, 'I2', i) )
        lib[('SpinorProp', i, 1)] = create( SpinorPropagator('I2', i, i) )
        lib[('Vector', i)] = create( Vector(i+1, i+1) )
        lib[('VectorProp', i)] = create( VectorPropagator(i,'I2', i) )
        lib[('Spin2', i )] = create( Spin2(10*i+1, 10*i+2, i) )
        lib[('Spin2Prop',i)] = create( Spin2Propagator(10*i+1, \
                                            10*i+2,'I2','I3', i) )
    logger.info('writing')         
    fsock = open('./HelasLib.pkl','wb')
    cPickle.dump(lib, fsock, -1)
    logger.info('done')
    


       
if '__main__' == __name__:
    logging.basicConfig(level=0)
    #create_library()
    import cProfile
    start = time.time()
    helasgenerator = AbstractHelasModel('sm')
    helasgenerator.compute_all()
    logger.info('done')
    stop = time.time()
    logger.info('done in %s s' % (stop-start))








