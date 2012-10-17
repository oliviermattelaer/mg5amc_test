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

import aloha.create_aloha as create_aloha
import models.write_param_card as param_writer
from madgraph import MadGraph5Error, MG5DIR
from madgraph.iolibs.files import cp, ln, mv
_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0] + '/'
logger = logging.getLogger('madgraph.export_v4')

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
    def write_co_flow_v4(self, writer, flow, helas_call_writer, me_number = ''):
        """Export a matrix element to a flow.f file in MG4 standalone format"""

        if not flow.get('processes') or \
               not flow.get('diagrams'):
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

        # Extract nwavefuncs
        nwavefuncs = flow.get_number_of_wavefunctions()
        replace_dict['nwavefuncs'] = nwavefuncs

        # Extract IC data line
        ic_data_line = self.get_ic_data_line(flow)
        replace_dict['ic_data_line'] = ic_data_line

        # Extract helas calls
        helas_calls = helas_call_writer.get_matrix_element_calls(flow)
        replace_dict['helas_calls'] = "\n".join(helas_calls)

        # Extract JAMP lines
        jamp_lines = self.get_JAMP_lines(flow)
        replace_dict['jamp_lines'] = '\n'.join(jamp_lines)

        file = open(os.path.join(_file_path, \
                                 'color_ordering/template_files/%s' % \
                                 self.co_flow_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return len(filter(lambda call: call.find('#') != 0, helas_calls))

    def get_ic_data_line(self, flow):
        """Get the IC line, giving sign for external HELAS wavefunctions"""

        ret_line = "DATA IC/"
        for wf in flow.get_external_wavefunctions():
            if wf.is_fermion():
                # For fermions, need particle/antiparticle
                ret_line += "%d" % (- (-1) ** wf.get('is_part'))
            else:
                # For boson, need initial/final
                ret_line += "%d" % ((-1) ** (wf.get('state') == 1))
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
        needed_perms = sorted(list(set([icol / nflows for (icol, irow) in \
                                        color_matrix.keys()])))
        nperms = len(needed_perms)
        all_perms = matrix_element.get('permutations')

        # Get maximum color factor Nc
        max_Nc = max([max([c.Nc_power for c in color_matrix[(icol, irow)]]) \
                      for (icol, irow) in color_matrix.keys()])
        
        # Get maximum color factor for each flow in each perm
        # (for comments see below)
        perm_flow_factors = {}
        for (icol, irow) in sorted(color_matrix.keys()):
            iperm = icol / nflows
            iflow = icol % nflows
            flow_Nc = max([c.Nc_power for c in \
                                  color_matrix[(icol, irow)]]) - max_Nc
            perm_flow_factors[(iperm,iflow)] = \
              max(flow_Nc, perm_flow_factors.setdefault((iperm,iflow), flow_Nc))

        # Each row in the color matrix corresponds to one of the basic
        # flows, while each column corresponds to a flow. Keep track
        # of the needed flows for each permutation
        perm_needed_flows = {}
        row_flow_factors = {}
        # Mapping from basic flows (in 1st perm) to jamp number
        flow_jamp_dict = {}
        # nflows_needed keeps track of the present JAMP number
        jamp = 0
        for (icol, irow) in sorted(color_matrix.keys()):
            # irow is the number of the basic flow (from first permutation)
            # iperm is the permutation (among the full set, all_perms)
            iperm = icol / nflows
            # iflow is the flow number (for this permutation)
            iflow = icol % nflows

            # Calculate Nc for this flow in this row
            row_Nc = max([c.Nc_power for c in color_matrix[(icol, irow)]])
            flow_Nc = row_Nc - max_Nc

            # Add this flow to the needed flows for this permutation
            # (used for the flow call lines generated below)
            if not iflow in [i for (i,n,c) in \
                             perm_needed_flows.setdefault(iperm, [])]:
                jamp += 1
                perm_needed_flows[iperm].append((iflow, jamp,
                                         perm_flow_factors[(iperm, iflow)]))
                if iperm == 0: flow_jamp_dict[iflow] = jamp

            # Make sure that also the basic flow is included
            if not irow in [i for (i,n,c) in perm_needed_flows[0]]:
                jamp += 1
                perm_needed_flows[0].append((irow, jamp, 
                                             perm_flow_factors[(0, irow)]))
                flow_jamp_dict[irow] = jamp
            # Add the factor needed for this JAMP
            row_flow_factors.setdefault(irow, []).append(\
                            (row_Nc,
                             jamp if iperm > 0 else flow_jamp_dict[iflow],
                             color_matrix.col_matrix_fixed_Nc[(icol, irow)][0],
                             flow_Nc))


        return jamp, needed_perms, perm_needed_flows, row_flow_factors, \
               flow_jamp_dict

    def get_perm_lines(self, matrix_element, needed_perms):
        """Get the permutations needed to calculate JAMPs for this
        color order"""

        all_perms = matrix_element.get('permutations')
        nexternal = len(all_perms[0])

        # The data lines giving the needed permutations
        iperm_line_list = []
        for iperm, perm in enumerate(needed_perms):
            int_list = [iperm+1, nexternal]
            int_list.extend(all_perms[perm])
            iperm_line_list.append(\
                ("DATA (PERMS(I,%4r),I=1,%d) /" + \
                 ",".join(['%2r'] * nexternal) + "/") % tuple(int_list))

        return iperm_line_list
    
    def get_iferm_line(self, matrix_element, needed_perms):
        """Get the fermion factors for the needed permutations"""

        # The data line for iferm
        iferm_list = []
        all_perms = matrix_element.get('permutations')
        external_fermions = [i for (i,w) in enumerate(matrix_element.\
                             get_external_wavefunctions()) if w.is_fermion()]
        for perm in needed_perms:
            fermion_numbers = [all_perms[perm][i] for i in external_fermions]
            iferm_list.append(helas_objects.HelasAmplitude.sign_flips_to_order(\
                fermion_numbers))

        iferm_line = "DATA IFERM/" + \
                     ",".join(['%2r' % i for i in iferm_list]) + "/"

        return iferm_line

    def get_flow_call_lines(self, needed_perms, perm_needed_flows,
                            me_flag = ''):
        """Write out the calls to all color flows. Need to multiply by
        fermion permutation factor for this color flow to get right
        sign."""

        # Generate the calls to all needed flows, by color order
        flow_call_lines = []
        # First calculate minimum color order from color orders in 
        # perm_needed_flows
        min_color_order = min(sum([[c for (i,j,c) in perm_needed_flows[key]] \
                                    for key in perm_needed_flows.keys()], []))

        # Write out the calls to the color flows, order by order
        for color_order in range(0, min_color_order - 1, -2):
            # We only want to separate odd orders, since even
            # correspond to singlet gluon contributions only
            flow_call_lines.append("IF(ICO.EQ.%d) THEN" % \
                                       (1 - (color_order / 2)))
            for iperm, perm in enumerate(needed_perms):
                orders = max([c/2 for (i,j,c) in perm_needed_flows[perm]])
                # Only include permutations with relevant flows
                if color_order/2 > orders: continue

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

    def get_color_flow_lines(self, row_flow_factors, flow_jamp_dict):
        """Write summation of all color flows. Need to multiply by
        fermion permutation factor for this color flow to get right
        sign."""
        
        # The color matrix summation lines for the basic color flows
        color_sum_lines = []
        
        min_color_order = min(sum([[n for (i,j,c,n) in row_flow_factors[key]] \
                                    for key in row_flow_factors.keys()], []))

        # Go through the rows and output the explicit color matrix
        # summation for this line
        for color_order in range(0, min_color_order - 1, -2):
            color_sum_lines.append("IF(ICO.EQ.%d) THEN" % \
                                       (1 - (color_order / 2)))
            for irow in sorted(row_flow_factors.keys()):
                orders = [n for (i,j,c,n) in row_flow_factors[irow] if \
                          n/2 == color_order/2]
                # Only include lines with relevant flows
                if not orders: continue
                # Get denominator and flows for this color_order
                den, factor_dict = self.organize_row(row_flow_factors[irow],
                                                     color_order)
                color_sum_lines.append(\
                    'ZTEMP = ZTEMP+%(den)s*JAMP(%(jamp)d)*DCONJG(%(flows)s)' % \
                    {'den': self.fraction_to_string(den),
                     'jamp': flow_jamp_dict[irow],
                     'flows': "+".join(['%s*(%s)' % \
                                    (self.fraction_to_string(fact),\
                                     "+".join(["%d*JAMP(%d)" % i for i in \
                                               factor_dict[fact]])) for fact \
                                    in sorted(factor_dict.keys(), reverse=True)])})
                color_sum_lines[-1] = color_sum_lines[-1].replace('+-1*', '-')
                color_sum_lines[-1] = color_sum_lines[-1].replace('+1*', '+')
                color_sum_lines[-1] = color_sum_lines[-1].replace('(-1*', '(-')
                color_sum_lines[-1] = color_sum_lines[-1].replace('(1*', '(')
                color_sum_lines[-1] = color_sum_lines[-1].replace('+-', '-')
                color_sum_lines[-1] = color_sum_lines[-1].replace('+1D0*', '+')
                color_sum_lines[-1] = color_sum_lines[-1].replace('/1*', '*')
            color_sum_lines.append("ENDIF")
        return color_sum_lines

    @staticmethod
    def organize_row(flow_factors, color_order):
        """Organize the information for this row to get a nice output.
        The elements of flow_factors is Nc_power, jamp number, fraction.
        Return the common denominator and a dictionary from value to
        sorted list of jamp numbers. Only include jamps with the correct
        color order (color_order or color_order + 1)"""

        # First pick out only the relevant factors, based on color_order
        orders = [(i,j,c) for (i,j,c,n) in flow_factors if \
                  n/2 == color_order/2]
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
    def generate_subprocess_directory_v4(self, matrix_element,
                                         helas_call_writer):
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
            matrix_element)

        # Create the flow.f files for each color flow
        calls = 0
        for flow in matrix_element.get('color_flows'):
            filename = 'flow%d.f' % flow.get('number')
            calls += self.write_co_flow_v4(
                writers.FortranWriter(filename),
                flow,
                helas_call_writer)

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

        linkfiles = ['check_sa.f', 'ipnext.f', 'coupl.inc', 'makefile']


        for file in linkfiles:
            ln('../%s' % file)

        # Return to original PWD
        os.chdir(cwd)

        if not calls:
            calls = 0
        return calls

    #===========================================================================
    # write_matrix_element_v4
    #===========================================================================
    def write_matrix_element_v4(self, writer, matrix_element):
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

        # Extract flow function definition lines
        flow_functions_lines = self.get_flow_functions_lines(matrix_element)
        replace_dict['flow_functions_lines'] = flow_functions_lines

        # Extract nperms
        replace_dict['nperms'] = len(matrix_element.get('permutations'))       

        # Extract call lines and color sum lines
        nflows, needed_perms, perm_needed_flows, row_flow_factors, \
                flow_jamp_dict = self.get_flow_info(matrix_element)
        nflowperms = len(needed_perms)
        flow_perms_data_lines = self.get_perm_lines(matrix_element,
                                                    needed_perms)
        flow_iferm_data_line = self.get_iferm_line(matrix_element,
                                                   needed_perms)
        flow_call_lines = self.get_flow_call_lines(needed_perms,
                                                   perm_needed_flows)
        color_sum_lines = self.get_color_flow_lines(row_flow_factors,
                                                    flow_jamp_dict)

        replace_dict['flow_perms_data_lines'] = '\n'.join(flow_perms_data_lines)
        replace_dict['flow_iferm_data_line'] = flow_iferm_data_line
        replace_dict['nflowperms'] = nflowperms
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)
        replace_dict['nflows'] = nflows
        replace_dict['color_order'] = matrix_element.get('color_order')

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
    def generate_subprocess_directory_v4(self, matrix_element,
                                         co_helas_call_writer,
                                         me_number):
        """Generate the Pxxxxx directory for a subprocess in MG4 madevent,
        including the necessary matrix.f and various helper files"""

        cwd = os.getcwd()
        path = os.path.join(self.dir_path, 'SubProcesses')

        try:
             os.chdir(path)
        except OSError, error:
            error_msg = "The directory %s should exist in order to be able " % path + \
                        "to \"export\" in it. If you see this error message by " + \
                        "typing the command \"export\" please consider to use " + \
                        "instead the command \"output\". "
            raise MadGraph5Error, error_msg 


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
                co_helas_call_writer)

        filename = 'auto_dsig.f'
        self.write_auto_dsig_file(writers.FortranWriter(filename),
                             matrix_element)

        filename = 'configs.inc'
        mapconfigs, s_and_t_channels = self.write_configs_file(\
            writers.FortranWriter(filename),
            matrix_element)

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

        filename = 'symfact.dat'
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
                     'ipnext.f',
                     'makefile',
                     'message.inc',
                     'myamp.f',
                     'reweight.f',
                     'run.inc',
                     'maxconfigs.inc',
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

        # Write maxconfigs.inc based on max of ME's/subprocess groups
        filename = os.path.join(self.dir_path,'Source','maxconfigs.inc')
        self.write_maxconfigs_file(writers.FortranWriter(filename),
                                   matrix_elements)
        
        # Write maxparticles.inc based on max of ME's/subprocess groups
        filename = os.path.join(self.dir_path,'Source','maxparticles.inc')
        self.write_maxparticles_file(writers.FortranWriter(filename),
                                     matrix_elements)
        
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

    def write_matrix_element_v4(self, writer, matrix_element, helas_call_writer,
                                proc_id = "", config_map = []):
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
        flow_iferm_data_line = self.get_iferm_line(matrix_element,
                                                   needed_perms)
        flow_call_lines = self.get_flow_call_lines(needed_perms,
                                                   perm_needed_flows,
                                                   me_flag)
        color_sum_lines = self.get_color_flow_lines(row_flow_factors,
                                                    flow_jamp_dict)
        replace_dict['flow_perms_data_lines'] = '\n'.join(flow_perms_data_lines)
        replace_dict['flow_iferm_data_line'] = flow_iferm_data_line
        replace_dict['nflowperms'] = nflowperms
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)
        replace_dict['njamps'] = njamps

        # Extract JAMP2 summation lines, using only leading color flows
        nflows, jamp2_lines = self.get_jamp2_lines(matrix_element,
                                                   perm_needed_flows)
        replace_dict['nflows'] = nflows
        replace_dict['jamp2_lines'] = '\n'.join(jamp2_lines)

        # Extract the info about which particles should be permuted
        comp_data_line = self.get_comp_data_line(matrix_element)
        replace_dict['comp_data_line'] = comp_data_line

        # Set color order
        replace_dict['color_order'] = matrix_element.get('color_order')

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
        max_Nc = max([flows[iflow].get('color_string').Nc_power \
                      for (iflow, j, co) in perm_needed_flows[0]])
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
    def generate_subprocess_directory_v4(self, subproc_group,
                                         co_helas_call_writer,
                                         group_number):

        matrix_elements = subproc_group.get('matrix_elements')

        # First generate all files needed except for the flow files
        export_v4.ProcessExporterFortranMEGroup.\
                      generate_subprocess_directory_v4(self,
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
                    me_flag)

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
            raise self.PhysicsObjectError, \
                  "get_helas_call must be called with wavefunction or amplitude"
        
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

