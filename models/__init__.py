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
"""All models for MG5, in particular UFO models (by FeynRules)"""

import os
import sys

def load_model(name, decay=False):
    
    # avoid final '/' in the path
    if name.endswith('/'):
        name = name[:-1]
    
    path_split = name.split(os.sep)
    if len(path_split) == 1:
        model_pos = 'models.%s' % name
        __import__(model_pos)
        return sys.modules[model_pos]
    elif path_split[-1] in sys.modules:
        if not sys.modules[path_split[-1]].__file__.startswith(os.sep.join(path_split[:-1])):
            raise Exception, 'name %s already consider as a python library cann\'t be reassigned' % \
                path_split[-1] 

    sys.path.insert(0, os.sep.join(path_split[:-1]))
    __import__(path_split[-1])
    output = sys.modules[path_split[-1]]
    if decay:
        dec_name = '%s.decays' % path_split[-1]
        __import__(dec_name)
        output.all_decays = sys.modules[dec_name].all_decays
    
    sys.path.pop(0)
    
    
    
    return sys.modules[path_split[-1]]
