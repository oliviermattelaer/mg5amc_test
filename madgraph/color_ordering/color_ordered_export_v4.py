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
"""Methods and classes to export matrix elements to v4 format."""

from __future__ import absolute_import
import copy
import fractions
import glob
import logging
import os
import re
import shutil
import subprocess

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.color_amp as color_amp
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.drawing_eps as draw
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.files as files
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.gen_infohtml as gen_infohtml
import madgraph.iolibs.template_files as iolibs_template_files
import madgraph.color_ordering.color_ordered_amplitudes as \
       color_ordered_amplitudes
import madgraph.color_ordering.color_ordered_helas_objects as \
       color_ordered_helas_objects
import madgraph.color_ordering.template_files as template_files
import madgraph.iolibs.ufo_expression_parsers as parsers
import madgraph.various.diagram_symmetry as diagram_symmetry
import madgraph.various.misc as misc

import madgraph.color_ordering.List_py.list_fundamental_1qq as list_NLC

import aloha.create_aloha as create_aloha
import models.write_param_card as param_writer
from madgraph import MadGraph5Error, MG5DIR
from madgraph.iolibs.files import cp, ln, mv
from six.moves import range
_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0] + '/'
logger = logging.getLogger('madgraph.export_v4')

pjoin = os.path.join
#===============================================================================
# ProcessExporterFortranCO
#===============================================================================
class ProcessExporterFortranCO(export_v4.ProcessExporterFortran):
    """Base class to take care of exporting a set of color ordered matrix
    elements to MadGraph v4 format."""

    co_flow_file = "co_flow_v4.inc"

    #===========================================================================
    # write_co_flow_v4
    #===========================================================================
    def write_co_flow_v4(self, writer, flow, helas_call_writer, matrix_element, me_number = '',):
        """Export a matrix element to a flow.f file in MG4 standalone format"""

        if not flow.get('processes') or \
               not flow.get('diagrams'):
            return 0

        if not isinstance(writer, writers.FortranWriter):
            raise writers.FortranWriter.FortranWriterError(\
                "writer not FortranWriter")

    #    misc.sprint(type(flow))
    #    misc.sprint(flow)
        # Set lowercase/uppercase Fortran code
        writers.FortranWriter.downcase = False

        replace_dict = {}
        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(flow)
        replace_dict['process_lines'] = process_lines

        # Set number
        replace_dict['number'] = "%s%d" % (me_number, flow.get('number'))

        # Extract number of external particles
        (nexternal, ninitial) = flow.get_nexternal_ninitial()
        replace_dict['nexternal'] = nexternal

        # Extract ngraphs
        ngraphs = flow.get_number_of_amplitudes()
        replace_dict['ngraphs'] = ngraphs

        # Extract nperms, perms
        replace_dict['nperms'] = len(flow.get('permutations')) 
        perms = matrix_element.get('permutations')

        reduced_rows = True


        if reduced_rows:
            perms = list_NLC.list_jamps('gluons',nexternal,ninitial)
        else:
            perms = list_NLC.list_jamps_all_rows('gluons',nexternal,ninitial)

        # Extract IC data line
        ic_data_line = self.get_ic_data_line(flow)
        replace_dict['ic_data_line'] = ic_data_line

        # Extract helas calls
        helas_calls = helas_call_writer.get_matrix_element_calls(flow)
        replace_dict['helas_calls'] = "\n".join(helas_calls)
        # misc.sprint(helas_calls)
        
        # Extract nwavefuncs
        nwavefuncs = flow.get_number_of_wavefunctions()
        replace_dict['nwavefuncs'] = nwavefuncs

        # Extract fermion permutation factors
        iferm_lines = self.get_all_iferm_lines(matrix_element,nexternal,ninitial)
        replace_dict['iferm_lines_all'] = iferm_lines
        # misc.sprint(iferm_lines)

        # Extract new stuff that AL put in after perms
        helas_calls_2, calls_dict, nwavefuncs, namps = self.get_permuted_helas_calls(helas_calls,flow, perms)
        replace_dict['helas_calls'] = "\n".join(helas_calls_2)
        # misc.sprint(helas_calls_2)
        replace_dict['nwavefuncs'] = nwavefuncs
        replace_dict['ngraphs'] = namps

        # Extract number of jamps
        replace_dict['njampsAL'] = len(matrix_element.get('permutations'))\
                                   *len(matrix_element.get('color_flows'))

        
        if reduced_rows:
            njamps = len(list_NLC.list_jamps('gluon',nexternal,ninitial))
        else:
            njamps = len(list_NLC.list_jamps_all_rows('gluon',nexternal,ninitial))

        replace_dict['njampsAL'] = njamps


        # Extract JAMP lines

        jamp_lines, nb_tmp_jamp = self.get_JAMP_lines(flow)
        jamp_lines_temp = self.get_permuted_jamp_lines(jamp_lines, calls_dict)

        # Test jamp permutation dict
    #    self.get_jamp_dict(flow)

        replace_dict['jamp_lines'] = '\n'.join(jamp_lines_temp)



        replace_dict['nb_temp_jamp'] = nb_tmp_jamp
        #adding the support for the fake width (forbidding too small width)
        mass_width = flow.get_all_mass_widths()
        width_list = set([e[1] for e in mass_width])
        
        replace_dict['fake_width_declaration'] = \
            ('  double precision fk_%s \n' * len(width_list)) % tuple(width_list)
        replace_dict['fake_width_declaration'] += \
            ('  save fk_%s \n' * len(width_list)) % tuple(width_list)
        fk_w_defs = []
        one_def = ' IF(%(w)s.ne.0d0) fk_%(w)s = SIGN(MAX(ABS(%(w)s), ABS(%(m)s*small_width_treatment)), %(w)s)'    
        for m, w in mass_width:
            if w == 'zero':
                if ' fk_zero = 0d0' not in fk_w_defs: 
                    fk_w_defs.append(' fk_zero = 0d0')
                continue    
            fk_w_defs.append(one_def %{'m':m, 'w':w})
        replace_dict['fake_width_definitions'] = '\n'.join(fk_w_defs)


        file = open(os.path.join(_file_path, \
                                 'color_ordering/template_files/%s' % \
                                 self.co_flow_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)
        # misc.sprint([call for call in helas_calls if call.find('#') != 0])

        return len([call for call in helas_calls if call.find('#') != 0])


    def get_permuted_helas_calls(self, helas_calls, flow, perms):
        """Function to go over all permutations and output the relevant
           helas_calls to flow.f"""

        nperms = len(perms)
        (nexternal, ninitial) = flow.get_nexternal_ninitial()
        flow_num = flow.get('number')

   
    #    misc.sprint(nperms, perms, flow_num, type(flow))

        # replace IP(num) with num
        # TODO: Learn how to REGEX to get rid of everything except the number
        helas_calls = [call.replace('IP','') for call in helas_calls]
        # helas_calls_noIP2 = [re.sub(, call) for call in helas_calls]

        # Come up with a list of dictionaries of permutations of wfs, amps, and jamps
        perm_dicts = {}

        # first go through initial permutation
        jamp_num = 1 + nperms*(flow_num-1)
        wf_dict = {}
        amp_dict = {}

        # get wavefunctions from first permutation
        for iwf, wf in enumerate(flow.get('diagrams')[0].get('wavefunctions')):
            wf_dict['W(1,'+ str(wf.get('number')) + ')'] = \
                    'Wn(1,' + str(wf.get('number')) + ')'

        # keep track of number of wavefunctions
        nwfs = len(wf_dict)

        # get helas calls for wfs
        helas_calls_wfs = copy.copy(helas_calls[:nwfs])
        nwfs_in_perm = len(helas_calls_wfs)
        
        # replace W -> Wn for first permutation in helas_calls_wfs
        for iwf in range(len(flow.get('diagrams')[0].get('wavefunctions'))):
            for wf_n in wf_dict:
                if wf_n in helas_calls_wfs[iwf]:
                    helas_calls_wfs[iwf] = helas_calls_wfs[iwf].replace(wf_n, wf_dict[wf_n])

        # get amplitudes from first permutation
        for diag in flow.get('diagrams'):
            for amp in diag.get('amplitudes'):
                amp_dict['AMP(' + str(amp.get('number')) + ')'] = \
                    'AMPN(' + str(amp.get('number')) + ')'

        # keep track of number of amps
        namps = len(amp_dict)
    #    misc.sprint(namps)

        # get helas calls for amps
        helas_calls_amps = copy.copy(helas_calls[nwfs:])
        helas_calls_amps_only = [call for call in helas_calls_amps if "AMP" in call]
        namps_in_perm = len(helas_calls_amps_only)
        # misc.sprint(namps_in_perm, helas_calls_amps_only)
        
        # # replace W -> Wn and AMP -> AMPN for first permutation in helas_calls_amps
        # for icall, call in enumerate(helas_calls_amps):
        #     if 'CALL' in call:
        #         for wf_n in wf_dict:
        #             helas_calls_amps[icall] = helas_calls_amps[icall].\
        #                 replace(wf_n, wf_dict[wf_n])
        #         for amp_n in amp_dict:
        #             helas_calls_amps[icall] = helas_calls_amps[icall].\
        #                 replace(amp_n, amp_dict[amp_n])
                    

        # replace W -> Wn and AMP -> AMPN for first permutation in helas_calls_amps
        for icall, call in enumerate(helas_calls_amps_only):
            for wf_n in wf_dict:
                helas_calls_amps_only[icall] = helas_calls_amps_only[icall].\
                                          replace(wf_n, wf_dict[wf_n])
            for amp_n in amp_dict:
                helas_calls_amps_only[icall] = helas_calls_amps_only[icall].\
                                          replace(amp_n, amp_dict[amp_n])



        # add information to dictionary of jamps to wfs, amps
        perm_dicts[jamp_num] = wf_dict, amp_dict


        # loop through rest of perms
        for iperm, perm in enumerate(perms):
            # first perm already done
            if iperm == 0: continue

            # reset wf and amp dicts
            wf_dict = {}
            amp_dict = {}

            # permute wavefunctions by first permuting external particles and adding 
            # these permutations to the dictionary, then for each internal particle
            # replace external particles in helas_call, check if this helas_call 
            # already exists, and add the permutation of the internal particle wavefunction
            # to the dictionary.

            # first permute external wavefunctions and store in dictionary
            for iwf in range(nexternal):
                wf_dict['W(1,'+ str(iwf+1) + ')'] = \
                       'Wn(1,' + str(perm[iwf]) + ')'

            # make a copy of the helas_calls 
            helas_calls_copy = copy.copy(helas_calls)

            # now go through internal particles
            for iwf in range(nexternal, nwfs_in_perm):
                # loop over wfs already in dictionary and replace with new wf
                for wf_n in wf_dict:
                    if wf_n in helas_calls_copy[iwf]:
                        helas_calls_copy[iwf] = helas_calls_copy[iwf].\
                            replace(wf_n, wf_dict[wf_n])
                
                # now check if internal wavefunction iwf already exists by checking 
                # its call. If it doesn't exist, give it a new number and add it to helas_calls_wfs. 
                # In both cases update the wf dictionary
                
                # get current call minus the last wf
                last_wf = helas_calls_copy[iwf].split(',')[-2:]
                last_wf = ','.join(last_wf)
                curr_call_less_wf = helas_calls_copy[iwf].replace(last_wf,'')
                # loop over already used calls to see if this call exists
                for call in helas_calls_wfs[nexternal:]:
                    # misc.sprint(call)
                    # if call exists, update dictionary to map to that call
                    if curr_call_less_wf in call:
                        last_wf_call = call.split(',')[-2:]
                        last_wf_call = ','.join(last_wf_call)
                        # remove extra ) from wf
                        last_wf_call = last_wf_call[:-1]
                        wf_dict['W(1,' + str(iwf+1) + ')'] = last_wf_call
                        break
                    # else call doesn't exist yet, create it and give wf new number
                else:
                    nwfs += 1
                    iW = 'W(1,' + str(iwf+1) + ')'
                    iWN = 'Wn(1,' + str(nwfs) + ')'
                    wf_dict[iW] = iWN
                    helas_calls_copy[iwf] = helas_calls_copy[iwf].replace(iW,iWN)
                    helas_calls_wfs.append(helas_calls_copy[iwf])

            # now go through amps and add new amps to dictionary/helas_calls

            # make a copy of the helas_calls for amps
            helas_calls_copy = copy.copy(helas_calls[nwfs_in_perm:])
            helas_calls_copy = [call for call in helas_calls_copy if "AMP" in call]

            for iamp in range(namps_in_perm):
                # first go through wavefunctions and get permuted ones
                for wf_n in wf_dict:
                    helas_calls_copy[iamp] = helas_calls_copy[iamp].\
                                             replace(wf_n,wf_dict[wf_n])
                
                # now check if amp iamp already exists by checking its call
                # If it doesn't exist, give it a new number and add it to helas_calls_amps. 
                # In both cases update the amp dictionary
               
                # get current call minus amp
                amp = helas_calls_copy[iamp].split(',')[-1]
                curr_call_less_amp = helas_calls_copy[iamp].replace(amp,'')
                
                # loop over already used calls to see if this call exists
                for call in helas_calls_amps_only:
                    # misc.sprint(call)
                    # if call exists, update dictionary to map to that call
                    if curr_call_less_amp in call:
                        amp = call.split(',')[-1]
                        iAMP = 'AMP(' + str(iamp+1) + ')'
                        iAMPN = amp[:-1]
                        # amp_dict[iamp] = amp
                        amp_dict[iAMP] = iAMPN
                        break
                    # else call doesn't exist yet, create it and give amp new number
                else:
                    namps += 1
                    iAMP = 'AMP(' + str(iamp+1) + ')'
                    iAMPN = 'AMPN(' + str(namps) + ')'
                    amp_dict[iAMP] = iAMPN
                    helas_calls_copy[iamp] = helas_calls_copy[iamp].replace(iAMP,iAMPN)
                    # not sure yet if want to add description so just append to both helas
                    # calls of amps for now
                    helas_calls_amps_only.append(helas_calls_copy[iamp])
                    # helas_calls_amps.append(helas_calls_copy[iamp])
        
                    
            # add dictionaries to jamp dictionary
            perm_dicts[iperm + 1 + nperms*(flow_num-1)] = wf_dict, amp_dict
        
        # misc.sprint(helas_calls_wfs)
        # misc.sprint(helas_calls_amps_only)


        # return variable = list of permuted helas_calls (W->Wn, AMP->AMPN)
        helas_calls_permuted = []
        helas_calls_permuted.extend([ wf_call for \
                         wf_call in helas_calls_wfs ])
        helas_calls_permuted.extend([ amp_call for \
                         amp_call in helas_calls_amps_only ])
    #    misc.sprint(helas_calls_permuted)



        return helas_calls_permuted, perm_dicts, nwfs, namps

    def get_permuted_jamp_lines(self, jamp_lines, jamp_dict):
        """return jamp lines after W->WN, AMP -> AMPN permutations"""
        
        # get list of jamps
        jamp_lines_ret = ['JAMP(' + str(i) + ')' for i in jamp_dict]
    #    misc.sprint(jamp_lines_ret)
        
        # loop through jamps and use dictionary to replace amps in amp_sum with 
        # those from dictionary
    #    misc.sprint(jamp_dict)
        # for ijamp in range(len(jamp_dict)):
        # counter for jamps in given flow file
        jcounter = 0
        for ijamp in jamp_dict:
            # get unpermuted amps
            amp_sum = jamp_lines[0].replace('JAMP(1)','')
            # multiply amp sum by +-1 depending on fermion perms
            iferm_line = 'IFERM(' + str(ijamp) + ')*('
            amp_sum = amp_sum[:1] + iferm_line + amp_sum[1:] + ')'
            # misc.sprint(amp_sum)
            # get dictionary with amps for this jamp
            curr_dict = jamp_dict[ijamp]
    #        misc.sprint(ijamp, len(jamp_dict), curr_dict)
            # replace amps with permuted amps
            for key in curr_dict[1]:
                # misc.sprint(type(key), key, type(amp_sum), amp_sum)
                # misc.sprint(curr_dict)
                val = curr_dict[1][key]
                if key in amp_sum:
                    amp_sum = amp_sum.replace(key, val)
            
            # now add amp_sum to jamp
            jamp_lines_ret[jcounter] += amp_sum
            jcounter +=1

        return jamp_lines_ret
    
    def get_jamp_dict(self, flow):
        """Get dictionary of jamp permutations for a given row of the 
        colour matrix. I.e., the standard trace basis colour matrix """

        (nexternal, ninitial) = flow.get_nexternal_ninitial()
        perms = flow.get('permutations')

        nperms = len(flow.get('permutations'))

    #    misc.sprint(nexternal,ninitial)
    #    misc.sprint(perms)
    #    misc.sprint(nperms)

        if reduced_rows:
            perms = list_NLC.list_jamps('gluons',nexternal,ninitial)
        else:
            perms = list_NLC.list_jamps_all_rows('gluons',nexternal,ninitial)

        nperms = len(perms)

    #    misc.sprint(perms)
    #    misc.sprint(nperms)

        # make dict of permuted jamps
        jamp_perm_dict = {}
        
        # loop over permutations/jamps
        for ijamp in range(1,nperms+1):
            # make dictionary of external wavefunction permutations
            w_dict_tmp = {}
            for i in range(nexternal):
                w_dict_tmp[i+1] = perms[ijamp-1][i]
            # for ikey, key in enumerate(jamp_dict[ijamp][0]):
            #     # only want to track permutations of external particles
            #     if ikey >= nexternal: 
            #         continue
            #     val = jamp_dict[ijamp][0][key]
            #     # get just numbers, not W(1,num) etc.
            #     # first get num)
            #     key = key.split(',')[1]
            #     val = val.split(',')[1]
            #     # now get num
            #     key = int(key.split(')')[0])
            #     val = int(val.split(')')[0])
            #     w_dict_tmp[key] = val
            #     # misc.sprint(key,val)
            # misc.sprint(w_dict_tmp)

            # loop over perms for this ijamp to find which corresponds to which jamp
    #        misc.sprint(w_dict_tmp)

            jamp_perms_temp = {}
            # get a list of external particle orderings for each jamp
            ext_nums = []
            for i in range(nperms):
                ext_nums_tmp = []
                for inum in perms[i]:
                    ext_nums_tmp.append(w_dict_tmp[inum])
                ext_nums.append(ext_nums_tmp)
        #    misc.sprint(ext_nums)
        #    misc.sprint(len(ext_nums))
        #    misc.sprint('passed') 
            
            # put the list of particle orderings into the jamp dictionary
            ind_list = []
            for iperm, perm in enumerate(perms):
                ind_list.append(ext_nums.index(perm)+1)

            #    misc.sprint(iperm, ind_list)
                jamp_perms_temp[iperm+1] = ind_list[iperm]

            jamp_perm_dict[ijamp] = jamp_perms_temp
        #    misc.sprint(jamp_perm_dict)
    #    misc.sprint(jamp_perm_dict, type(flow))       

    #    misc.sprint(jamp_perm_dict)
        return jamp_perm_dict

    def get_jamp_data_lines(self, matrix_element):
        """Get the permutations needed to permute JAMPs"""
        
        # get jamp dict to be printed as DATA
        permuted_jamp_dict = self.get_jamp_dict(matrix_element)
        nflows = len(matrix_element.get('color_flows'))
    #    misc.sprint(type(matrix_element), nflows)
    #    misc.sprint(permuted_jamp_dict)
        nperms = len(permuted_jamp_dict)
        
        # return data a list
        jamp_data_list = []

        for iflow in range(0,nflows):
          for irow in permuted_jamp_dict:
              # get list of permuted jamps for this row
              jperm_list = [int(val + iflow*nperms) for val in permuted_jamp_dict[irow].values()]
            #   jperm_list = iflow
              jamp_data_list.append(\
                  "DATA (JPERM(I, %d),I=%d,%d)"  % (irow,1+nperms*iflow,nperms*(iflow+1))
                  + " /" + ",".join(['%d' % iperm for iperm in jperm_list]) + "/" )

        return jamp_data_list

    def get_ic_data_line(self, flow):
        """Get the IC line, giving sign for external HELAS wavefunctions"""

        ret_line = "DATA IC/"
        for wf in flow.get_external_wavefunctions():
            if wf.is_fermion():
                # For fermions, need particle/antiparticle
                ret_line += "%d" % (- (-1) ** wf.get('is_part'))
            else:
                # For boson, need initial/final
                ret_line += "%d" % ((-1) ** (wf.get('state') == "initial"))
            ret_line += ','

        ret_line = ret_line[:-1] + '/'
        
        return ret_line

    def get_flow_info(self, matrix_element):
        """Return the lines for defining all needed
        permutations, fermion factors, and the color sum lines. Note
        that this is all only for one row in the color matrix (per
        basic color flow)."""

        color_matrix = matrix_element.get('color_matrix')


        flows = matrix_element.get('color_flows')
        
        nflows = len(flows)

        # The permutations with non-zero color matrix elements with
        # the basic color flows


        needed_perms = sorted(list(set([icol // nflows for (icol, irow) in \
                                        color_matrix.keys()])))

        nperms = len(needed_perms)

        n_all_perms = len(matrix_element.get('permutations'))

        # Get maximum color factor Nc
        max_Nc = max([max([c.Nc_power for c in color_matrix[(icol, irow)][0]]) \
                      for (icol, irow) in color_matrix.keys()])

        # Get maximum color factor for each flow in each perm
        # (for comments see below)
        perm_flow_factors = {}

    #    for (icol, irow) in sorted(color_matrix.keys()):
        for (icol, irow) in color_matrix.keys():

            iperm = icol // nflows
            iflow = icol % nflows


            flow_Nc = max([c.Nc_power for c in color_matrix[(icol, irow)][0]]) - max_Nc

            perm_flow_factors[(iperm, irow, iflow)] = \
              max(flow_Nc, perm_flow_factors.setdefault((iperm, irow, iflow), flow_Nc))

        # Each row in the color matrix corresponds to one of the basic
        # flows, while each column corresponds to a flow. Keep track
        # of the needed flows for each permutation

        perm_needed_flows = {}
        row_flow_factors = {}
        # Mapping from basic flows (in 1st perm) to jamp number
        flow_jamp_dict = {}
        # nflows_needed keeps track of the present JAMP number
        jamp = 0

#        for (icol, irow) in sorted(color_matrix.keys()):
        for (icol, irow) in color_matrix.keys():

            # irow is the number of the basic flow (from first permutation)
            # iperm is the permutation (among the full set, all_perms)

            iperm = icol // nflows
            # iflow is the flow number (for this permutation)
            iflow = icol % nflows

            # Calculate Nc for this flow in this row
            row_Nc = max([c.Nc_power for c in color_matrix[(icol, irow)][0]])
            flow_Nc = row_Nc - max_Nc

            # Add this flow to the needed flows for this permutation
            # (used for the flow call lines generated below)


            if not iflow in [i for (i,n,c) in \
                             perm_needed_flows.setdefault(iperm, [])]:
                jamp += 1
                # perm_needed_flows[iperm].append((iflow, jamp,
                #                          perm_flow_factors[(iperm, iflow)]))
                # AL: Keep perm number rather than starting from 1 again

                perm_needed_flows[iperm].append((iflow, iperm,
                                         perm_flow_factors[(iperm, irow, iflow)]))

                #if iperm == 0: 
                flow_jamp_dict[irow] = jamp

            # Make sure that also the basic flow is included
            #if not irow in [i for (i,n,c) in perm_needed_flows[0]]:
            #    jamp += 1
            #    # perm_needed_flows[0].append((irow, jamp, 
            #    #                              perm_flow_factors[(0, irow)]))
            #    # AL: Keep perm number rather than starting from 1 again

            #    perm_needed_flows[0].append((irow, iperm, 
            #                                 perm_flow_factors[(irow, 0)]))
            #    flow_jamp_dict[irow] = jamp

            # Add the factor needed for this JAMP
        #    misc.sprint(row_flow_factors.setdefault(irow, []))
            row_flow_factors.setdefault(irow, []).append(\
                            (row_Nc,
                            # AL: Keep perm number rather than starting from 1 again
                            #  jamp if iperm > 0 else flow_jamp_dict[iflow],

                             (iperm+1)+ n_all_perms*iflow,
                            # (iperm+1) + nperms*iflow if iperm > 0 else flow_jamp_dict[iflow],
                             color_matrix.col_matrix_fixed_Nc[(icol, irow)][0],
                             flow_Nc))


        return jamp, needed_perms, perm_needed_flows, row_flow_factors, \
               flow_jamp_dict

    def get_perm_lines(self, matrix_element, needed_perms):
        """Get the permutations needed to calculate JAMPs for this
        color order"""

        all_perms = matrix_element.get('permutations')
        nexternal = len(all_perms[0])
    #    misc.sprint(needed_perms)

        # The data lines giving the needed permutations
        iperm_line_list = []
        for iperm, perm in enumerate(needed_perms):
            int_list = [iperm+1, nexternal]
            # int_list = [needed_perms[iperm]+1, nexternal]
            int_list.extend(all_perms[perm])
            iperm_line_list.append(\
                ("DATA (PERMS(I,%4r),I=1,%d) /" + \
                 ",".join(['%2r'] * nexternal) + "/") % tuple(int_list))

        return iperm_line_list

    def get_all_iferm_lines(self, matrix_element,nexternal,ninitial):
        """Get the fermion factors for the needed permutations"""

        # The data line for iferm
        iferm_list = []


        reduced_rows = True

        if reduced_rows:
            all_perms = list_NLC.list_jamps('gluons',nexternal,ninitial)
        else:
            all_perms = matrix_element.get('permutations')

        external_fermions = [i for (i,w) in enumerate(matrix_element.\
                             get_external_wavefunctions()) if w.is_fermion()]
        
        for iflow, flow in enumerate(matrix_element.get('color_flows')):
            for perm in all_perms:
                fermion_numbers = [perm[i] for i in external_fermions]
                iferm_list.append(helas_objects.HelasAmplitude.sign_flips_to_order(\
                    fermion_numbers))

        iferm_line = "DATA IFERM/" + \
                     ",".join(['%2r' % i for i in iferm_list]) + "/"
        
        return iferm_line
    
    def get_iferm_line(self, matrix_element, needed_perms):
        """Get the fermion factors for the needed permutations"""

        # The data line for iferm
        iferm_list = []
        all_perms = matrix_element.get('permutations')
        # misc.sprint(all_perms, needed_perms, type(matrix_element))
        external_fermions = [i for (i,w) in enumerate(matrix_element.\
                             get_external_wavefunctions()) if w.is_fermion()]
        for perm in needed_perms:
            fermion_numbers = [all_perms[perm][i] for i in external_fermions]
            iferm_list.append(helas_objects.HelasAmplitude.sign_flips_to_order(\
                fermion_numbers))

        iferm_line = "DATA IFERM/" + \
                     ",".join(['%2r' % i for i in iferm_list]) + "/"

        return iferm_line

    def get_flow_call_lines(self, needed_perms, perm_needed_flows, min_color_order,
                            me_flag = ''):
        """Write out the calls to all color flows. Need to multiply by
        fermion permutation factor for this color flow to get right
        sign."""

        # Generate the calls to all needed flows, by color order
        flow_call_lines = []
        
        # Write out the calls to the color flows, order by order
        for color_order in range(0, int(min_color_order) - 2, -2):
            # We only want to separate odd orders, since even
            # correspond to singlet gluon contributions only
            flow_call_lines.append("IF(ICO.EQ.%d) THEN" % \
                                       (1 - (int(color_order) // 2)))
            for iperm, perm in enumerate(needed_perms):
                orders = max([c/2 for (i,j,c) in perm_needed_flows[perm]])
                # Only include permutations with relevant flows
                if int(color_order)//2 > orders: continue

                # Set the perm needed in the flow calls
                flow_call_lines.append("DO I=1,NEXTERNAL")
                flow_call_lines.append("PERM(I)=PM(PERMS(I,%d))" % \
                                       (iperm + 1))
                flow_call_lines.append("ENDDO")
                # Now generate the flow calls
                for iflow, jmp, co in perm_needed_flows[perm]:
                    if co < color_order: 
                        continue
                    flow_call_lines.append(\
                        "JAMP(%d)=IFERM(%d)*FLOW%s%d(P,NHEL,PERM)" \
                        % (jmp, iperm+1, me_flag, iflow+1))

            if color_order % 2 == 0:
                flow_call_lines.append("ENDIF")

        return flow_call_lines





    def get_color_flow_lines(self, row_flow_factors, flow_jamp_dict, min_color_order, matrix_element):
        """Write summation of all color flows. Need to multiply by
        fermion permutation factor for this color flow to get right
        sign."""
        
### Alwaus use NLC setup
        min_color_order = -2

        # The color matrix summation lines for the basic color flows
        color_sum_lines = []

        nperms = len(matrix_element.get('permutations'))
        
        # Go through the rows and output the explicit color matrix
        # summation for this line

        rows = len(row_flow_factors)
        cols = len(row_flow_factors[0])

        for color_order in range(0, min_color_order - 2, -2):
            if color_order == 0:
                color_sum_lines.append("IF(ICO.EQ.1) THEN" )
            else:
                color_sum_lines.append("ELSE IF(ICO.EQ.%d) THEN" % \
                                       (1 - (color_order / 2)))

            # Loop over rows of colour matrix
            color_sum_lines.append("DO I=1, %d" % rows)
            color_sum_lines.append("DO J=1, %d" % cols)
      

            icol = 0

            #for irow in sorted(row_flow_factors.keys()):

            irow = 0

            orders = [n for (i,j,c,n) in row_flow_factors[irow] if \
                          int(n)//2 == int(color_order)//2]

            # Only include lines with relevant flows
            if not orders: continue

        # Get denominator and flows for this color_order

            den, factor_dict = self.organize_row(row_flow_factors[irow],
                                                     color_order)

        #        color_sum_lines.append(\
        #            'ZTEMP = ZTEMP+%(den)s*JAMP(JPERM(%(jamp)d,I))*DCONJG(%(flows)s)' % \
        #            {'den': self.fraction_to_string(den),
        #            # 'jamp': irow*nperms+1,
        #              'jamp': flow_jamp_dict[irow],
        #             'flows': "+".join(['%s*(%s)' % \
        #                            (self.fraction_to_string(fact),\
        #                             "+".join(["%d*JAMP(JPERM(%d,I))" % i for i in \
        #                                       factor_dict[fact]])) for fact \
        #                            in sorted(list(factor_dict.keys()), reverse=True)])})

            color_sum_lines.append(\
                    'ZTEMP = ZTEMP+JAMP(I)*DCONJG(SYM*CF(J,I)*JAMP(LOC(J,I)))')

            color_sum_lines[-1] = color_sum_lines[-1].replace('+-1*', '-')
            color_sum_lines[-1] = color_sum_lines[-1].replace('+1*', '+')
            color_sum_lines[-1] = color_sum_lines[-1].replace('(-1*', '(-')
            color_sum_lines[-1] = color_sum_lines[-1].replace('(1*', '(')
            color_sum_lines[-1] = color_sum_lines[-1].replace('+-', '-')
            color_sum_lines[-1] = color_sum_lines[-1].replace('+1D0*', '+')
            color_sum_lines[-1] = color_sum_lines[-1].replace('/1*', '*')

            color_sum_lines.append("ENDDO")
            color_sum_lines.append("ENDDO")


        if color_order <= min_color_order:
                color_sum_lines.append("ENDIF")

        return color_sum_lines

    
    def get_loc_data_lines(self, matrix_element, n=6):
        """Return the color matrix definition lines for this matrix element. Split
        rows in chunks of size n."""

        col_matrix = matrix_element.get('color_matrix')

        col_length = max(col_matrix.keys())[0]
        row_length = max(col_matrix.keys())[1]


        if not matrix_element.get('color_matrix'):
            return ["DATA Denom(1)/1/", "DATA (LOC(i,1),i=1,1) /1/"]
        else:
            ret_list = []
            my_cs = color.ColorString()
            index = 0

            for irow in range(row_length+1):
                
                col_list = [str(col_matrix[(i,irow)][1]+1) for i in range(col_length+1)]

                # First write the common denominator for this color matrix line
                #ret_list.append("DATA Denom(%i)/%i/" % (index + 1, denominator))
                # Then write the numerators for the matrix elements

                ret_list.append("DATA (LOC(I,%3r),i=%3r,%3r) /%s/" % \
                                 (irow+1, 1,len(col_list),','.join(col_list)))
                                 

            #    ret_list.append("C %s" % repr(my_cs))

            return ret_list


    def get_color_data_lines_nlc(self, matrix_element, n=6):
        """Return the color matrix definition lines for this matrix element. Split
        rows in chunks of size n."""

        col_matrix = matrix_element.get('color_matrix')

        col_length = max(col_matrix.keys())[0]
        row_length = max(col_matrix.keys())[1]

        if not col_matrix:
            return ["DATA Denom(1)/1/", "DATA (CF(i,1),i=1,1) /1/"]
        else:
            ret_list = []
            my_cs = color.ColorString()

            for irow in range(row_length+1):

                # First write the common denominator for this color matrix line
                #ret_list.append("DATA Denom(%i)/%i/" % (index + 1, denominator))
                # Then write the numerators for the matrix elements

                col_list = [str(float(col_matrix.col_matrix_fixed_Nc[(i, irow)][0])) for i in range(col_length+1)]

                ret_list.append("DATA (CF(I,%3r),I=%3r,%3r) /%s/" % \
                                    (irow+1, 1, len(col_list),','.join(col_list)))
                

            return ret_list


    @staticmethod
    def organize_row(flow_factors, color_order):
        """Organize the information for this row to get a nice output.
        The elements of flow_factors is Nc_power, jamp number, fraction.
        Return the common denominator and a dictionary from value to
        sorted list of jamp numbers. Only include jamps with the correct
        color order (color_order or color_order + 1)"""

        # First pick out only the relevant factors, based on color_order
        orders = [(i,j,c) for (i,j,c,n) in flow_factors if \
                  int(n)//2 == int(color_order)//2]
        assert(orders)
        co_flow_factors = orders
     
        # Then get common denominators for this row
        den = color_amp.ColorMatrix.lcmm(*[fact[2].denominator \
                                           for fact in co_flow_factors])
        if not den or den > 100000: den = 1
        return_dict = {}
        for facttuple in co_flow_factors:
            fact = facttuple[2]*den
            if fact == int(fact): fact = int(fact)
            return_dict.setdefault(abs(fact), []).append((fact/abs(fact),
                                                          facttuple[1]))
    #    misc.sprint(return_dict)

        return fractions.Fraction(1, int(den)), return_dict
        

    @staticmethod
    def Nc_power_to_fraction(Nc_power):
        """Given an Nc power, return the corresponding fraction"""
        if Nc_power < 0:
            return fractions.Fraction(1, 3**(-Nc_power))
        else:
            return fractions.Fraction(3**(Nc_power), 1)
        
    @staticmethod
    def fraction_to_string(fraction):
        """Return a Fortran string for a fraction"""
        if not fraction:
            return "0"
        return "%sD0" % str(fraction)
        
    def get_comp_data_line(self, matrix_element):
        """Write out the data line for the comp vector."""

        comp_list = matrix_element.get_comp_list()

        comp_data_line = 'DATA COMP/%s/' % ','.join([str(c) for c in comp_list])
        return comp_data_line
        
    def get_flow_functions_lines(self, matrix_element, me_flag = ''):
        """Write out function definition lines"""

        flow_functions = ",".join(["FLOW%s%d" % (me_flag, flow.get('number')) \
                                 for flow in matrix_element.get('color_flows')])

        return "COMPLEX*16 %(ff)s\nEXTERNAL %(ff)s" % {"ff": flow_functions}

    def get_call_flow_lines(self, matrix_element, me_flag = ''):
        """Write out calls to flow subroutine lines"""

        flow_calls=[]

        for flow in matrix_element.get('color_flows'):
            flow_calls.append('CALL FLOW%s%d(P,NHEL(:,IHEL), JAMP)'\
                                % (me_flag, flow.get('number')))

        return flow_calls    

    def combine_jamp_factors(self, allnums):
        """Combine all JAMPs with the same factors"""

        result = {}
        for fact, i in allnums:
            try:
                result[fact].append(i)
            except KeyError:
                result[fact] = [i]

        return result

#===============================================================================
# ProcessExporterFortranCOSA
#===============================================================================
class ProcessExporterFortranCOSA(export_v4.ProcessExporterFortranSA,
                                 ProcessExporterFortranCO):
    """Class to take care of exporting a set of color ordered matrix
    elements to MadGraph v4 StandAlone format."""

    super_matrix_file = "co_super_matrix_standalone_v4.inc"

    #===========================================================================
    # generate_subprocess_directory_v4
    #===========================================================================
    def generate_subprocess_directory(self, matrix_element,
                                         helas_call_writer, me=None):
        """Generate the Pxxxxx directory for a subprocess in MG4 standalone,
        including the necessary matrix.f and nexternal.inc files"""

        cwd = os.getcwd()

        # Create the directory PN_xx_xxxxx in the specified path
        dirpath = os.path.join(self.dir_path, 'SubProcesses', \
                       "P%s" % matrix_element.get('processes')[0].shell_string())

        try:
            os.mkdir(dirpath)
        except os.error as error:
            logger.warning(error.strerror + " " + dirpath)

        try:
            os.chdir(dirpath)
        except os.error:
            logger.error('Could not cd to directory %s' % dirpath)
            return 0

        logger.info('Creating files in directory %s' % dirpath)


        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()

        # Create the matrix.f file and the nexternal.inc file
        filename = 'matrix.f'
        self.write_matrix_element_v4(
            writers.FortranWriter(filename),
            matrix_element, None)

        # Create the flow.f files for each color flow
        calls = 0
        for flow in matrix_element.get('color_flows'):
        #    misc.sprint( flow.get('number'))
            filename = 'flow%d.f' % flow.get('number')
            calls += self.write_co_flow_v4(
                writers.FortranWriter(filename),
                flow,
                helas_call_writer,
                matrix_element)

        # Create the flow.ps files for each color flow
        calls = 0
        for flow in matrix_element.get('color_flows'):
            filename = 'flow%d.ps' % flow.get('number')
            plot = draw.MultiEpsDiagramDrawer(flow.get('base_amplitude').\
                                                 get('diagrams'),
                                              filename,
                                              model=matrix_element.get('processes')[0].\
                                                 get('model'),
                                              amplitude=True)
            logger.info("Generating Feynman diagrams for " + \
                         matrix_element.get('processes')[0].nice_string())
            plot.draw()

        filename = 'nexternal.inc'
        self.write_nexternal_file(writers.FortranWriter(filename),
                             nexternal, ninitial)

        filename = 'pmass.inc'
        self.write_pmass_file(writers.FortranWriter(filename),
                         matrix_element)

        linkfiles = ['check_sa.f', 'ipnext.f', 'coupl.inc', 'idenparts.f']

        for file in linkfiles:
            ln('../%s' % file)
        ln('../makefileP', name='makefile')
        # Return to original PWD
        os.chdir(cwd)

        if not calls:
            calls = -1
        return calls

    #===========================================================================
    # write_matrix_element_v4
    #===========================================================================
    def write_matrix_element_v4(self, writer, matrix_element, helascallwriter):
        """Export a matrix element to a matrix.f file in standalone
        color ordered amplitude format"""

        if not matrix_element.get('processes') \
               or not isinstance(matrix_element,
                              color_ordered_helas_objects.COHelasMatrixElement) \
               or not matrix_element.get('color_flows'):
            return 0

        if not isinstance(writer, writers.FortranWriter):
            raise writers.FortranWriter.FortranWriterError(\
                "writer not FortranWriter")

        # Set lowercase/uppercase Fortran code
        writers.FortranWriter.downcase = False

        replace_dict = {}

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(matrix_element)
        replace_dict['process_lines'] = process_lines

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
        replace_dict['nexternal'] = nexternal

        # Extract ncomb
        ncomb = matrix_element.get_helicity_combinations()
        replace_dict['ncomb'] = ncomb

        # Extract helicity lines
        helicity_lines = self.get_helicity_lines(matrix_element)
        replace_dict['helicity_lines'] = helicity_lines

        # Extract overall denominator
        # Averaging initial state color, spin, and identical FS particles
        den_factor_line = self.get_den_factor_line(matrix_element)
        replace_dict['den_factor_line'] = den_factor_line


        ##### Check if exporting color matrix still works!!!

        color_data_lines = self.get_color_data_lines_nlc(matrix_element)
        
        loc_lines =  self.get_loc_data_lines(matrix_element)

        replace_dict['loc_data_lines'] = "\n".join(loc_lines)

        reduced_rows = True

        if reduced_rows:
            nrows = len(list_NLC.list_rows('gluons',nexternal,ninitial))
        else:
            nrows = len(list_NLC.all_perms('gluons',nexternal,ninitial))

        replace_dict['nrows'] = nrows
        cols = list_NLC.permutations(list_NLC.list_rows('gluons',nexternal,ninitial)[0])
        replace_dict['ncols'] = len(cols)

        replace_dict['color_data_lines'] = "\n".join(color_data_lines)

        if reduced_rows:
            sym = nexternal -2
        else:
            sym = 1

        replace_dict['symmetry'] = self.factorial(sym)

        # Extract flow function definition lines
        flow_functions_lines = self.get_flow_functions_lines(matrix_element)
        replace_dict['flow_functions_lines'] = flow_functions_lines

        # Extract nperms
        replace_dict['nperms'] = len(matrix_element.get('permutations'))

        # Extract total number of Jamps
        replace_dict['njampsAL'] = len(matrix_element.get('permutations'))\
                                 *len(matrix_element.get('color_flows'))

        if reduced_rows:
            njamps = len(list_NLC.list_jamps('gluon',nexternal,ninitial))
        else:
            njamps = len(list_NLC.list_jamps_all_rows('gluon',nexternal,ninitial))

        replace_dict['njampsAL'] = njamps


    #    misc.sprint(replace_dict['njampsAL'])
    #    misc.sprint(len(matrix_element.get('color_flows'))   )
    #    misc.sprint(len(matrix_element.get('permutations')))


        # Extract call lines and color sum lines

        nflows, needed_perms, perm_needed_flows, row_flow_factors, \
                flow_jamp_dict = self.get_flow_info(matrix_element)

        nflowperms = len(needed_perms)

        flow_perms_data_lines = self.get_perm_lines(matrix_element,
                                                    needed_perms)
        flow_iferm_data_line  = self.get_iferm_line(matrix_element,
                                                   needed_perms)

    #    misc.sprint(sum([[n for (i,j,c,n) in row_flow_factors[key]]
    #                                for key in row_flow_factors.keys()], []))

        min_color_order = int(min(sum([[n for (i,j,c,n) in row_flow_factors[key]]
                                    for key in row_flow_factors.keys()], [])))
        flow_call_lines = self.get_flow_call_lines(needed_perms,perm_needed_flows, 
                                                   min_color_order)
    #    misc.sprint(flow_call_lines)

        color_sum_lines = self.get_color_flow_lines(row_flow_factors,
                                                    flow_jamp_dict, min_color_order,
                                                    matrix_element)
        # misc.sprint(color_sum_lines)

        call_flow_lines = self.get_call_flow_lines(matrix_element)

    
    #    jamp_perm_lines = self.get_jamp_data_lines(matrix_element)


        replace_dict['flow_perms_data_lines'] = '\n'.join(flow_perms_data_lines)
        replace_dict['flow_iferm_data_line'] = flow_iferm_data_line

        replace_dict['nflowperms'] = nflowperms
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)
        # AL: replaced color_sum_lines for new version
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)
        replace_dict['nflows'] = nflows
        replace_dict['color_order'] = matrix_element.get('color_order')
        replace_dict['call_flow_lines'] = '\n'.join(call_flow_lines)

    #    replace_dict['jamp_perm_lines'] = '\n'.join(jamp_perm_lines)
        replace_dict['jamp_perm_lines'] = '\n'


        # Extract the info about which particles should be permuted
        comp_data_line = self.get_comp_data_line(matrix_element)
        replace_dict['comp_data_line'] = comp_data_line

        # Create file contents
        file = open(os.path.join(_file_path, \
                                 'color_ordering/template_files/%s' % \
                                 self.super_matrix_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return 0

    def factorial(self,n):

        factorial = 1
        if int(n) >= 1:
          for i in range (1,int(n)+1):
            factorial = factorial * i

        return factorial


#===============================================================================
# ProcessExporterFortranCOME
#===============================================================================
class ProcessExporterFortranCOME(export_v4.ProcessExporterFortranME,
                                 ProcessExporterFortranCO):
    """Class to take care of exporting a set of matrix elements to
    MadEvent format."""

    matrix_file = "co_super_matrix_madevent_v4.inc"


    #===========================================================================
    # generate_subprocess_directory_v4 
    #===========================================================================
    def generate_subprocess_directory(self, matrix_element,
                                         co_helas_call_writer,
                                         me_number):
        """Generate the Pxxxxx directory for a subprocess in MG4 madevent,
        including the necessary matrix.f and various helper files"""

        cwd = os.getcwd()
        path = os.path.join(self.dir_path, 'SubProcesses')

        try:
             os.chdir(path)
        except OSError as error:
            error_msg = "The directory %s should exist in order to be able " % path + \
                        "to \"export\" in it. If you see this error message by " + \
                        "typing the command \"export\" please consider to use " + \
                        "instead the command \"output\". "
            raise MadGraph5Error(error_msg) 


        pathdir = os.getcwd()

        # Create the directory PN_xx_xxxxx in the specified path
        subprocdir = "P%s" % matrix_element.get('processes')[0].shell_string()
        try:
            os.mkdir(subprocdir)
        except os.error as error:
            logger.warning(error.strerror + " " + subprocdir)

        try:
            os.chdir(subprocdir)
        except os.error:
            logger.error('Could not cd to directory %s' % subprocdir)
            return 0

        logger.info('Creating files in directory %s' % subprocdir)


        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()

        # Create the matrix.f file, auto_dsig.f file and all inc files
        filename = 'matrix.f'
        calls,nc = self.write_matrix_element_v4(writers.FortranWriter(filename),
                                             matrix_element,
                                             co_helas_call_writer)



        # Create the flow.f files for each color flow
        calls = 0
        for flow in matrix_element.get('color_flows'):
            filename = 'flow%d.f' % flow.get('number')
            calls += self.write_co_flow_v4(
                writers.FortranWriter(filename),
                flow,
                co_helas_call_writer,
                matrix_element)

        filename = 'auto_dsig.f'
        self.write_auto_dsig_file(writers.FortranWriter(filename),
                             matrix_element)

        filename = 'configs.inc'
        mapconfigs, (s_and_t_channels, nqcd_list) = self.write_configs_file(\
            writers.FortranWriter(filename),
            matrix_element)

        filename = 'config_nqcd.inc'
        self.write_config_nqcd_file(writers.FortranWriter(filename),
                               nqcd_list)

        filename = 'config_subproc_map.inc'
        self.write_config_subproc_map_file(writers.FortranWriter(filename),
                                           s_and_t_channels)

        filename = 'coloramps.inc'
        self.write_coloramps_file(writers.FortranWriter(filename),
                             mapconfigs,
                             matrix_element)

        filename = 'get_color.f'
        self.write_colors_file(writers.FortranWriter(filename),
                               matrix_element)

        filename = 'decayBW.inc'
        self.write_decayBW_file(writers.FortranWriter(filename),
                           s_and_t_channels)

        filename = 'dname.mg'
        self.write_dname_file(writers.FortranWriter(filename),
                         matrix_element.get('processes')[0].shell_string())

        filename = 'iproc.dat'
        self.write_iproc_file(writers.FortranWriter(filename),
                         me_number)

        filename = 'leshouche.inc'
        self.write_leshouche_file(writers.FortranWriter(filename),
                             matrix_element)

        filename = 'maxamps.inc'
        self.write_maxamps_file(writers.FortranWriter(filename),
                           len(matrix_element.get('diagrams')),
                           len(matrix_element.get('color_flows')),
                           len(matrix_element.get('processes')),
                           1)

        filename = 'mg.sym'
        self.write_mg_sym_file(writers.FortranWriter(filename),
                          matrix_element)

        filename = 'ncombs.inc'
        self.write_ncombs_file(writers.FortranWriter(filename),
                          nexternal)

        filename = 'nexternal.inc'
        self.write_nexternal_file(writers.FortranWriter(filename),
                             nexternal, ninitial)

        filename = 'ngraphs.inc'
        self.write_ngraphs_file(writers.FortranWriter(filename),
                           len(mapconfigs))


        filename = 'pmass.inc'
        self.write_pmass_file(writers.FortranWriter(filename),
                         matrix_element)

        filename = 'props.inc'
        self.write_props_file(writers.FortranWriter(filename),
                         matrix_element,
                         s_and_t_channels)

        # Find config symmetries and permutations
        symmetry, perms, ident_perms = \
                  diagram_symmetry.find_symmetry(matrix_element)

        filename = 'symswap.inc'
        self.write_symswap_file(writers.FortranWriter(filename),
                                ident_perms)

        filename = 'symfact_orig.dat'
        self.write_symfact_file(writers.FortranWriter(filename),
                           symmetry)

        # Generate diagrams
        filename = "matrix.ps"
        plot = draw.MultiEpsDiagramDrawer(matrix_element.get('base_amplitude').\
                                             get('diagrams'),
                                          filename,
                                          model=matrix_element.get('processes')[0].\
                                             get('model'),
                                          amplitude='')
        logger.info("Generating Feynman diagrams for " + \
                     matrix_element.get('processes')[0].nice_string())
        plot.draw()

        # Create the flow.ps files for each color flow
        for flow in matrix_element.get('color_flows'):
            filename = 'flow%d.ps' % flow.get('number')
            plot = draw.MultiEpsDiagramDrawer(flow.get('base_amplitude').\
                                                 get('diagrams'),
                                              filename,
                                              model=matrix_element.get('processes')[0].\
                                                 get('model'),
                                              amplitude=True)
            logger.info("Generating Feynman diagrams for " + \
                         flow.get('processes')[0].nice_string())
            plot.draw()

        linkfiles = ['addmothers.f',
                     'cluster.f',
                     'cluster.inc',
                     'coupl.inc',
                     'cuts.f',
                     'cuts.inc',
                     'driver.f',
                     'genps.f',
                     'genps.inc',
                     'initcluster.f',
                     'idenparts.f',
                     'ipnext.f',
                     'makefile',
                     'message.inc',
                     'myamp.f',
                     'reweight.f',
                     'run.inc',
                     'maxconfigs.inc',
                     'run_config.inc',
                     'maxparticles.inc',
                     'setcuts.f',
                     'setscales.f',
                     'sudakov.inc',
                     'symmetry.f',
                     'unwgt.f']

        for file in linkfiles:
            ln('../' + file , '.')

        #import nexternal/leshouch in Source
        ln('nexternal.inc', '../../Source', log=False)
        ln('leshouche.inc', '../../Source', log=False)
        ln('maxamps.inc', '../../Source', log=False)

        # Return to SubProcesses dir
        os.chdir(pathdir)

        # Add subprocess to subproc.mg
        filename = 'subproc.mg'
        files.append_to_file(filename,
                             self.write_subproc,
                             subprocdir)

        # Return to original dir
        os.chdir(cwd)

        # Generate info page
        gen_infohtml.make_info_html(self.dir_path)


        if not calls:
            calls = 0
        return calls

    def finalize_v4_directory(self, matrix_elements, history, makejpg = False,
                              online = False, compiler = 'gfortran'):
        """Finalize ME v4 directory by creating jpeg diagrams, html
        pages,proc_card_mg5.dat and madevent.tar.gz."""

        #proc_charac
        self.create_proc_charac()
        self.create_run_card(matrix_elements, history)

        # Write maxconfigs.inc based on max of ME's/subprocess groups
        filename = os.path.join(self.dir_path,'Source','maxconfigs.inc')
        self.write_maxconfigs_file(writers.FortranWriter(filename),
                                   matrix_elements)
        
        # Write maxparticles.inc based on max of ME's/subprocess groups
        filename = os.path.join(self.dir_path,'Source','maxparticles.inc')
        self.write_maxparticles_file(writers.FortranWriter(filename),
                                     matrix_elements)
      
        # Add the combine_events.f modify param_card path/number of @X
        filename = pjoin(self.dir_path,'Source','combine_events.f')
        try:
            nb_proc =[p.get('id') for me in matrix_elements for m in me.get('matrix_elements') for p in m.get('processes')]
        except AttributeError:
            nb_proc =[p.get('id') for m in matrix_elements.get('matrix_elements') for p in m.get('processes')]
        nb_proc = len(set(nb_proc))
        self.write_combine_events(writers.FortranWriter(filename), nb_proc) # already formatted
                 
        # Touch "done" file
        os.system('touch %s/done' % os.path.join(self.dir_path,'SubProcesses'))

        # Check for compiler
        self.set_compiler(compiler)

        old_pos = os.getcwd()
        os.chdir(os.path.join(self.dir_path, 'SubProcesses'))
        P_dir_list = [proc for proc in os.listdir('.') if os.path.isdir(proc) and \
                                                                    proc[0] == 'P']

        devnull = os.open(os.devnull, os.O_RDWR)
        # Convert the poscript in jpg files (if authorize)
        if makejpg:
            logger.info("Generate jpeg diagrams")
            for Pdir in P_dir_list:
                os.chdir(Pdir)
                subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'internal', 'gen_jpeg-pl')],
                                stdout = devnull)
                os.chdir(os.path.pardir)

        logger.info("Generate web pages")
        # Create the WebPage using perl script

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'internal', 'gen_cardhtml-pl')], \
                                                                stdout = devnull)

        os.chdir(os.path.pardir)

        gen_infohtml.make_info_html(self.dir_path)

        [mv(name, './HTML/') for name in os.listdir('.') if \
                            (name.endswith('.html') or name.endswith('.jpg')) and \
                            name != 'index.html']               

        # Write command history as proc_card_mg5
        if os.path.isdir('Cards'):
            output_file = os.path.join('Cards', 'proc_card_mg5.dat')
            output_file = open(output_file, 'w')
            text = ('\n'.join(history) + '\n') % misc.get_time_info()
            output_file.write(text)
            output_file.close()

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'internal', 'gen_cardhtml-pl')],
                        stdout = devnull)

        # Generate madevent.tar.gz file
        if os.path.exists(os.path.join('SubProcesses', 'subproc.mg')):
            if os.path.exists('madevent.tar.gz'):
                os.remove('madevent.tar.gz')
            subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'internal', 'make_madevent_tar')],
                        stdout = devnull)


        if online:
            # Touch "Online" file
            os.system('touch %s/Online' % self.dir_path)

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin',
                                      'internal', 'gen_cardhtml-pl')],
                        stdout = devnull)

        #return to the initial dir
        os.chdir(old_pos)               


    #=============================================================================
    #  write_matrix_element_v4
    #=============================================================================
    def write_matrix_element_v4(self, writer, matrix_element, helas_call_writer,
                                proc_id = "", config_map = [], subproc_number=''):
        """Export a matrix element to a matrix.f file in MadEvent
        color ordered amplitude format"""


        if not matrix_element.get('processes') \
               or not isinstance(matrix_element,
                              color_ordered_helas_objects.COHelasMatrixElement) \
               or not matrix_element.get('color_flows'):
            return (0, 0)

        if not isinstance(writer, writers.FortranWriter):
            raise writers.FortranWriter.FortranWriterError(\
                "writer not FortranWriter")

        # Set lowercase/uppercase Fortran code
        writers.FortranWriter.downcase = False

        replace_dict = {}

        # Set proc_id
        replace_dict['proc_id'] = proc_id

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(matrix_element)
        replace_dict['process_lines'] = process_lines

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
        replace_dict['nexternal'] = nexternal

        # Extract ncomb
        ncomb = matrix_element.get_helicity_combinations()
        replace_dict['ncomb'] = ncomb

        # Extract IC data line
        ic_data_line = \
                     self.get_ic_data_line(matrix_element.get('color_flows')[0])
        replace_dict['ic_data_line'] = ic_data_line

        # Single permutation data line
        replace_dict['ip_data_line'] = "DATA IP/%s/" % \
                                       ','.join([str(i) for i in \
                                                 range(1, nexternal + 1)])

        # Extract helicity lines
        helicity_lines = self.get_helicity_lines(matrix_element)
        replace_dict['helicity_lines'] = helicity_lines

        # Extract overall denominator
        # Averaging initial state color, spin, and identical FS particles
        den_factor_line = self.get_den_factor_line(matrix_element)
        replace_dict['den_factor_line'] = den_factor_line

        # Extract ngraphs
        ngraphs = matrix_element.get_number_of_amplitudes()
        replace_dict['ngraphs'] = ngraphs

        # Extract ndiags
        ndiags = len(matrix_element.get('diagrams'))
        replace_dict['ndiags'] = ndiags

        # Set define_iconfigs_lines
        replace_dict['define_iconfigs_lines'] = \
             """INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
             COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG"""

        if proc_id:
            # Set lines for subprocess group version
            # Set define_iconfigs_lines
            replace_dict['define_iconfigs_lines'] += \
                 """\nINTEGER SUBDIAG(MAXSPROC),IB(2)
                 COMMON/TO_SUB_DIAG/SUBDIAG,IB"""    
            # Set set_amp2_line
            replace_dict['set_amp2_line'] = "ANS=ANS*AMP2(SUBDIAG(%s))/XTOT" % \
                                            proc_id
        else:
            # Standard running
            # Set set_amp2_line
            replace_dict['set_amp2_line'] = "ANS=ANS*AMP2(MAPCONFIG(ICONFIG))/XTOT"

        # Extract nwavefuncs
        nwavefuncs = matrix_element.get_number_of_wavefunctions()
        replace_dict['nwavefuncs'] = nwavefuncs

        # Extract helas calls
        helas_calls = helas_call_writer.get_matrix_element_calls(\
                    matrix_element)
        replace_dict['helas_calls'] = "\n".join(helas_calls)

        # Set the size of Wavefunction
        if not self.model or any([p.get('spin')==5 for p in self.model.get('particles') if p]):
            replace_dict['wavefunctionsize'] = 18
        else:
            replace_dict['wavefunctionsize'] = 6

        # Process specific flag for flow calls
        me_flag = ""
        if proc_id:
            me_flag = proc_id + "_"

        # Extract flow function definition lines
        flow_functions_lines = self.get_flow_functions_lines(matrix_element,
                                                             me_flag)
        replace_dict['flow_functions_lines'] = flow_functions_lines

        # Extract nperms
        replace_dict['nperms'] = len(matrix_element.get('permutations'))        


        # Extract flow call lines and color sum lines
        njamps, needed_perms, perm_needed_flows, row_flow_factors, \
                flow_jamp_dict = self.get_flow_info(matrix_element)
        nflowperms = len(needed_perms)


        flow_perms_data_lines = self.get_perm_lines(matrix_element,
                                                    needed_perms)
        flow_iferm_data_line  = self.get_iferm_line(matrix_element,
                                                   needed_perms)
        min_color_order = int(min(sum([[n for (i,j,c,n) in row_flow_factors[key]] \
                                    for key in row_flow_factors.keys()], [])))
        flow_call_lines = self.get_flow_call_lines(needed_perms,
                                                   perm_needed_flows, min_color_order,
                                                   me_flag)
        color_sum_lines = self.get_color_flow_lines(row_flow_factors,
                                                    flow_jamp_dict, min_color_order,
                                                    matrix_element)

    #    jamp_perm_lines = self.get_jamp_data_lines(matrix_element)

        call_flow_lines = self.get_call_flow_lines(matrix_element)
        replace_dict['flow_perms_data_lines'] = '\n'.join(flow_perms_data_lines)
        replace_dict['flow_iferm_data_line'] = flow_iferm_data_line
        replace_dict['nflowperms'] = nflowperms
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)
        # AL: replaced color_sum_lines for new version
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)

    #    replace_dict['jamp_perm_lines'] = '\n'.join(jamp_perm_lines)
        replace_dict['jamp_perm_lines'] = '\n'

        replace_dict['njamps'] = njamps
        replace_dict['call_flow_lines'] = '\n'.join(call_flow_lines)


        # Extract JAMP2 summation lines, using only leading color flows
        nflows, jamp2_lines = self.get_jamp2_lines(matrix_element,
                                                   perm_needed_flows)
        nflows = 1
    #    jamp2_lines = '\n'
        replace_dict['nflows'] = nflows
        replace_dict['jamp2_lines'] = '\n'.join(jamp2_lines)

        # Extract the info about which particles should be permuted
        comp_data_line = self.get_comp_data_line(matrix_element)
        replace_dict['comp_data_line'] = comp_data_line

        # Set color order
        replace_dict['color_order'] = matrix_element.get('color_order')

        #adding the support for the fake width (forbidding too small width)
        mass_width = matrix_element.get_all_mass_widths()
        width_list = set([e[1] for e in mass_width])
        
        replace_dict['fake_width_declaration'] = \
            ('  double precision fk_%s \n' * len(width_list)) % tuple(width_list)
        replace_dict['fake_width_declaration'] += \
            ('  save fk_%s \n' * len(width_list)) % tuple(width_list)
        fk_w_defs = []
        one_def = ' IF(%(w)s.ne.0d0) fk_%(w)s = SIGN(MAX(ABS(%(w)s), ABS(%(m)s*small_width_treatment)), %(w)s)'    
        for m, w in mass_width:
            if w == 'zero':
                if ' fk_zero = 0d0' not in fk_w_defs: 
                    fk_w_defs.append(' fk_zero = 0d0')
                continue    
            fk_w_defs.append(one_def %{'m':m, 'w':w})
        replace_dict['fake_width_definitions'] = '\n'.join(fk_w_defs)

        # Create file contents
        file = open(os.path.join(_file_path, \
                                 'color_ordering/template_files/%s' % \
                                 self.matrix_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return 0, max([len(flow.get_color_amplitudes()) for flow in \
                       matrix_element.get('color_flows')])

    def get_jamp2_lines(self, matrix_element, perm_needed_flows):
        """Extract the summation lines for JAMP2's, including only the
        leading color flows"""

        nflows = 0
        flows = matrix_element.get('color_flows')
        jamp2_lines = []


    #    max_Nc = max([flows[iflow].get('color_string').Nc_power \
    #                  for (iflow, j, co) in perm_needed_flows[0]])

        max_Nc= 6

        # The present color flows are in the first permutation
        for iflow, jamp, co in perm_needed_flows[0]:
            if flows[iflow].get('color_string').Nc_power == max_Nc:
                nflows += 1
                jamp2_lines.append(('JAMP2(%(nflow)d)=JAMP2(%(nflow)d)+' + \
                                    'JAMP(%(jamp)d)*DCONJG(JAMP(%(jamp)d))') % \
                                    {'nflow':nflows, 'jamp':jamp})

        return nflows, jamp2_lines

    def get_icolamp_lines(self, mapconfigs, matrix_element, num_matrix_element):
        """Return the ICOLAMP matrix, showing which JAMPs contribute to
        which configs (diagrams)."""

        ret_list = []

        booldict = {False: ".false.", True: ".true."}

        # No color, so only one color factor. Simply write a ".true." 
        # for each config (i.e., each diagram with only 3 particle
        # vertices
        configs = len(mapconfigs)
        ret_list.append("DATA(icolamp(1,i,%d),i=1,%d)/%s/" % \
                        (num_matrix_element, configs,
                         ','.join([".true." for i in range(configs)])))
        return ret_list

        # LATER ON MAKE SURE TO IMPLEMENT THIS PROPERLY!

#===============================================================================
# ProcessExporterFortranCOMEGroup
#===============================================================================
class ProcessExporterFortranCOMEGroup(export_v4.ProcessExporterFortranMEGroup,
                                      ProcessExporterFortranCOME):

    """Class to take care of exporting a set of matrix elements to
    MadEvent format."""

    matrix_file = "co_super_matrix_madevent_group_v4.inc"

    #===========================================================================
    # generate_subprocess_directory_v4
    #===========================================================================
    def generate_subprocess_directory(self, subproc_group,
                                         co_helas_call_writer,
                                         group_number):

        matrix_elements = subproc_group.get('matrix_elements')
        
        # First generate all files needed except for the flow files
        export_v4.ProcessExporterFortranMEGroup.\
                      generate_subprocess_directory(self,
                                                       subproc_group,
                                                       co_helas_call_writer,
                                                       group_number)

        # Create the flow.f files for each color flow
        calls = 0
        subprocdir = "P%d_%s" % (subproc_group.get('number'),
                                 subproc_group.get('name'))
        for ime,me in enumerate(subproc_group.get('matrix_elements')):
            # Process specific flag for flow calls
            me_flag = "%d_" % (ime + 1)
            for flow in me.get('color_flows'):
                filename = os.path.join(self.dir_path,
                                        'SubProcesses',
                                        subprocdir,
                                        'flow%s%d.f' % (me_flag,
                                                        flow.get('number')))
                calls += self.write_co_flow_v4(
                    writers.FortranWriter(filename),
                    flow,
                    co_helas_call_writer,
                    me_flag,
                    matrix_elements)

        # Rewrite maxamps.inc with correct values for flow
        filename = os.path.join(self.dir_path,
                                'SubProcesses',
                                subprocdir,
                                'maxamps.inc')

        self.write_maxamps_file(writers.FortranWriter(filename),
                           max([len(me.get('diagrams')) for me in \
                                        matrix_elements]),
                           max([len(me.get('color_flows')) for me in \
                                        matrix_elements]),
                           max([len(me.get('processes')) for me in \
                                matrix_elements]),
                           len(matrix_elements))



        return calls

    def write_matrix_element_v4(self, *args, **opts):
        """Just pass on to ProcessExporterFortranCOME"""

        return ProcessExporterFortranCOME.write_matrix_element_v4(\
                                                            self, *args, **opts)

    def get_icolamp_lines(self, *args, **opts):
        """Just pass on to ProcessExporterFortranCOME"""

        return ProcessExporterFortranCOME.get_icolamp_lines(self, *args, **opts)
        

#===============================================================================
# COFortranUFOHelasCallWriter
#===============================================================================
class COFortranUFOHelasCallWriter(helas_call_writers.FortranUFOHelasCallWriter):
    """The class for writing Helas calls in Fortran, starting from
    HelasWavefunctions and HelasAmplitudes. Include permutations for
    external wavefunctions, and BGHelasCurrent calls."""

    def generate_helas_call(self, argument):
        """Routine for automatic generation of Fortran Helas calls
        according to just the spin structure of the interaction.
        """

        if not isinstance(argument, helas_objects.HelasWavefunction) and \
           not isinstance(argument, helas_objects.HelasAmplitude):
            raise self.PhysicsObjectError("get_helas_call must be called with wavefunction or amplitude")
        
        call = "CALL "

        call_function = None

        if isinstance(argument, color_ordered_helas_objects.BGHelasCurrent):
            # Create call for wavefunction summation sumVN(fact1,W1,fact2,W2,...,Wres)
            call += "sum%s%d(" % \
                    (color_ordered_helas_objects.spin_dict[\
                                    argument.get('mothers')[0].get('spin')],
                                    len(argument.get('mothers')))
            call += "%s,W(1,%d)," * len(argument.get('mothers')) + \
                    "W(1,%d))"
            call_function = lambda wf: call % \
                (tuple(sum([[self.write_factor(mother.get('factor')),
                             mother.get('me_id')] for \
                            mother in wf.get('mothers')], []) + \
                [wf.get('me_id')]))
            self.add_wavefunction(argument.get_call_key(), call_function)
            return

        if isinstance(argument, helas_objects.HelasWavefunction) and \
               not argument.get('mothers'):
            # String is just IXXXXX, OXXXXX, VXXXXX or SXXXXX
            call = call + self.mother_dict[\
                argument.get_spin_state_number()]
            # Fill out with X up to 6 positions
            call = call + 'X' * (11 - len(call))
            call = call + "(P(0,IP(%d)),"
            if argument.get('spin') != 1:
                # For non-scalars, need mass and helicity
                call = call + "%s,NHEL(IP(%d)),"
            call = call + "%d*IC(IP(%d)),W(1,%d))"
            if argument.get('spin') == 1:
                call_function = lambda wf: call % \
                                (wf.get('number_external'),
                                 1,
                                 wf.get('number_external'),
                                 wf.get('me_id'))
            elif argument.is_boson():
                call_function = lambda wf: call % \
                                (wf.get('number_external'),
                                 wf.get('mass'),
                                 wf.get('number_external'),
                                 1,
                                 wf.get('number_external'),
                                 wf.get('me_id'))
            else:
                call_function = lambda wf: call % \
                                (wf.get('number_external'),
                                 wf.get('mass'),
                                 wf.get('number_external'),
                                 wf.get('fermionflow'),
                                 wf.get('number_external'),
                                 wf.get('me_id'))
        # Add the constructed function to wavefunction or amplitude dictionary
            self.add_wavefunction(argument.get_call_key(), call_function)
            return

        # By default, call mother function
        super(COFortranUFOHelasCallWriter, self).generate_helas_call(argument)
