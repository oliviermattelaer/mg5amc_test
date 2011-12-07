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
import madgraph.iolibs.misc as misc
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.gen_infohtml as gen_infohtml
import madgraph.iolibs.template_files as iolibs_template_files
import madgraph.color_ordering.color_ordered_amplitudes as \
       color_ordered_amplitudes
import madgraph.color_ordering.template_files as template_files
import madgraph.iolibs.ufo_expression_parsers as parsers
import madgraph.various.diagram_symmetry as diagram_symmetry

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

    co_flow_file = "co_flow_standalone_v4.inc"

    #===========================================================================
    # write_co_flow_v4
    #===========================================================================
    def write_co_flow_v4(self, writer, flow, fortran_model):
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
        replace_dict['number'] = flow.get('number')

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
        helas_calls = fortran_model.get_matrix_element_calls(flow)
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
                ret_line += "%d" % ((-1) ** (wf.get('state') == 'initial'))
            ret_line += ','

        ret_line = ret_line[:-1] + '/'
        
        return ret_line

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
    def generate_subprocess_directory_v4(self, matrix_element, fortran_model):
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
                fortran_model)

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

        linkfiles = ['check_sa.f', 'coupl.inc', 'makefile']


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
                              color_ordered_amplitudes.COHelasMatrixElement) \
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

        # Extract flow call lines and color sum lines
        nflows, flow_call_lines, color_sum_lines, nflowperms, \
                flow_perms_data_lines, flow_iferm_data_line = \
                self.get_flow_call_lines(matrix_element)
        replace_dict['flow_perms_data_lines'] = '\n'.join(flow_perms_data_lines)
        replace_dict['flow_iferm_data_line'] = flow_iferm_data_line
        replace_dict['nflowperms'] = nflowperms
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)
        replace_dict['nflows'] = nflows


        # Extract permutation data lines
        perms_data_lines, nperms = self.get_perms_data_lines(matrix_element)
        replace_dict['nperms'] = nperms
        replace_dict['perms_data_lines'] = "\n".join(perms_data_lines)

        # Create file contents
        file = open(os.path.join(_file_path, \
                                 'color_ordering/template_files/%s' % \
                                 self.super_matrix_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return 0


    def get_perms_data_lines(self, matrix_element):
        """Write out the data lines defining the permutations"""

        iperm_line_list = []
        i = 0
        for perm in matrix_element.get('permutations'):
            i = i + 1
            int_list = [i, len(perm)]
            int_list.extend(perm)
            iperm_line_list.append(\
                ("DATA (PERMS(I,%4r),I=1,%d) /" + \
                 ",".join(['%2r'] * len(perm)) + "/") % tuple(int_list))

        return iperm_line_list, len(iperm_line_list), 
        
    def get_flow_functions_lines(self, matrix_element):
        """Write out function definition lines"""

        flow_functions = ",".join(["FLOW%d" % flow.get('number') for flow in \
                                   matrix_element.get('color_flows')])

        return "COMPLEX*16 %(ff)s\nEXTERNAL %(ff)s" % {"ff": flow_functions}
        
    def get_flow_call_lines(self, matrix_element):
        """Write out the calls to all color flows. Need to multiply by
        fermion permutation factor for this color flow to get right
        sign. Also return the lines for defining all needed
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
        nexternal, ninitial = matrix_element.get_nexternal_ninitial()
        # The data lines giving the needed permutations
        iperm_line_list = []
        for iperm, perm in enumerate(needed_perms):
            int_list = [iperm+1, nexternal]
            int_list.extend(all_perms[perm])
            iperm_line_list.append(\
                ("DATA (PERMS(I,%4r),I=1,%d) /" + \
                 ",".join(['%2r'] * nexternal) + "/") % tuple(int_list))

        # The data line for iferm
        iferm_list = []
        external_fermions = [i for (i,w) in enumerate(matrix_element.\
                             get_external_wavefunctions()) if w.is_fermion()]
        for perm in needed_perms:
            fermion_numbers = [all_perms[perm][i] for i in external_fermions]
            iferm_list.append(helas_objects.HelasAmplitude.sign_flips_to_order(\
                fermion_numbers))

        iferm_line = "DATA IFERM/" + \
                     ",".join(['%2r' % i for i in iferm_list]) + "/"

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
            # Add this flow to the needed flows for this permutation
            # (used for the flow call lines generated below)
            if not iflow in [i for (i,n) in \
                             perm_needed_flows.setdefault(iperm, [])]:
                jamp += 1
                perm_needed_flows[iperm].append((iflow, jamp))
                if iperm == 0: flow_jamp_dict[iflow] = jamp

            # Make sure that also the basic flow is included
            if not irow in [i for (i,n) in perm_needed_flows[0]]:
                jamp += 1
                perm_needed_flows[0].append((irow, jamp))
                flow_jamp_dict[irow] = jamp
            # Add the factor needed for this JAMP
            row_flow_factors.setdefault(irow, []).append(\
                               (max([c.Nc_power for c in \
                                     color_matrix[(icol, irow)]]),
                                jamp if iperm > 0 else flow_jamp_dict[iflow],
                                color_matrix.col_matrix_fixed_Nc[(icol, irow)]))

        # Generate the calls to all needed flows
        flow_call_lines = []
        for iperm, perm in enumerate(needed_perms):
            # Set the perm needed in the flow calls
            flow_call_lines.append("DO I=1,NEXTERNAL")
            flow_call_lines.append("PERM(I)=PM(PERMS(I,%d))" % (iperm + 1))
            flow_call_lines.append("ENDDO")
            # Now generate the flow calls
            for iflow, jmp in perm_needed_flows[perm]:
                flow_call_lines.append(\
                    "JAMP(%d)=IFERM(%d)*FLOW%d(P,NHEL,PERM)" \
                    % (jmp, iperm+1, iflow+1))

        # Now output the color matrix summation lines for the basic color flows
        color_sum_lines = []
        
        # Go through the rows and output the explicit color matrix
        # summation for this line
        for irow in sorted(row_flow_factors.keys()):
            row_flow_factors[irow].sort(lambda c1,c2: c2[0]-c1[0] if \
                                        c1[0]-c2[0] != 0 else c1[1]-c2[1])
            color_sum_lines.append(\
                'ZTEMP = ZTEMP+JAMP(%(jamp)d)*DCONJG(%(flows)s)' % \
                {'jamp': flow_jamp_dict[irow],
                 'flows': "+".join(['%s*(%s)' % \
                                    (self.fraction_to_string(fact[0]),\
                                     "JAMP(%d)" % jamp) for (n, jamp, fact) \
                                     in row_flow_factors[irow]])})
            color_sum_lines[-1] = color_sum_lines[-1].replace('+-', '-')
            color_sum_lines[-1] = color_sum_lines[-1].replace('+1D0*', '+')
            color_sum_lines[-1] = color_sum_lines[-1].replace('/1*', '*')

        return jamp, flow_call_lines, color_sum_lines, nperms, \
               iperm_line_list, iferm_line

    def get_jamp_factors_for_row(self, matrix_element, color_matrix, irow,
                                 row_Nc_power, numerators):
        """Get the series of (numerator, JAMP) for this flow number"""

        color_basis = matrix_element.get('color_basis')
        flows = matrix_element.get('color_flows')
        nperms = len(matrix_element.get('permutations'))
        allnums = []
        for icol, col_bas_elem in enumerate(sorted(color_basis.keys())):
            for diag_tuple in color_basis[col_bas_elem]:
                iflow = diag_tuple[0]
                iperm = diag_tuple[1][0]
                nflow = 1 + iflow*nperms + iperm
                column_Nc_power = flows[diag_tuple[0]].get('color_string').\
                                  Nc_power
                if numerators[icol]:
                    # Only add if Nc_power large enough. Count only
                    # half of negative Nc powers, since the singlet contribution
                    # is put at -2 but should be -1
                    Nc_power = color_matrix[(irow,icol)][0].Nc_power + \
                           (row_Nc_power + column_Nc_power)//2
                    if Nc_power >= matrix_element.get('min_Nc_power'):
                        Nc_pow = self.Nc_power_to_fraction(column_Nc_power)
                        allnums.append((numerators[icol] * Nc_pow,
                                        nflow))
        
        return allnums

    def combine_jamp_factors(self, allnums):
        """Combine all JAMPs with the same factors"""

        result = {}
        for fact, i in allnums:
            try:
                result[fact].append(i)
            except KeyError:
                result[fact] = [i]

        return result
    
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
        return "%sd0" % str(fraction)
        

#===============================================================================
# ProcessExporterFortranCOME
#===============================================================================
class ProcessExporterFortranCOME(export_v4.ProcessExporterFortranME,
                                 ProcessExporterFortranCO):
    """Class to take care of exporting a set of matrix elements to
    MadEvent format."""

    matrix_file = "matrix_madevent_v4.inc"

    #===========================================================================
    # generate_subprocess_directory_v4 
    #===========================================================================
    def generate_subprocess_directory_v4(self, matrix_element,
                                         fortran_model,
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
        calls, ncolor = \
               self.write_matrix_element_v4(writers.FortranWriter(filename),
                                                matrix_element,
                                                fortran_model)

        filename = 'auto_dsig.f'
        self.write_auto_dsig_file(writers.FortranWriter(filename),
                             matrix_element)

        filename = 'configs.inc'
        mapconfigs, s_and_t_channels = self.write_configs_file(\
            writers.FortranWriter(filename),
            matrix_element)

        filename = 'coloramps.inc'
        self.write_coloramps_file(writers.FortranWriter(filename),
                             mapconfigs,
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
                           ncolor,
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
                     'makefile',
                     'message.inc',
                     'myamp.f',
                     'reweight.f',
                     'run.inc',
                     'maxconfigs.inc',
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
                              online = False):
        """Finalize ME v4 directory by creating jpeg diagrams, html
        pages,proc_card_mg5.dat and madevent.tar.gz."""

        # Write maxconfigs.inc based on max of ME's/subprocess groups
        filename = os.path.join(self.dir_path,'Source','maxconfigs.inc')
        self.write_maxconfigs_file(writers.FortranWriter(filename),
                                   matrix_elements)
        
        # Touch "done" file
        os.system('touch %s/done' % os.path.join(self.dir_path,'SubProcesses'))

        if not misc.which('g77'):
            logger.info('Change makefiles to use gfortran')
            subprocess.call(['python','./bin/Passto_gfortran.py'], cwd=self.dir_path, \
                            stdout = open(os.devnull, 'w')) 

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
                subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'gen_jpeg-pl')],
                                stdout = devnull)
                os.chdir(os.path.pardir)

        logger.info("Generate web pages")
        # Create the WebPage using perl script

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'gen_cardhtml-pl')], \
                                                                stdout = devnull)

        os.chdir(os.path.pardir)

        gen_infohtml.make_info_html(self.dir_path)
        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'gen_crossxhtml-pl')],
                        stdout = devnull)
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

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'gen_cardhtml-pl')],
                        stdout = devnull)

        # Run "make" to generate madevent.tar.gz file
        if os.path.exists(os.path.join('SubProcesses', 'subproc.mg')):
            if os.path.exists('madevent.tar.gz'):
                os.remove('madevent.tar.gz')
            subprocess.call(['make'], stdout = devnull)


        if online:
            # Touch "Online" file
            os.system('touch %s/Online' % self.dir_path)

        subprocess.call([os.path.join(old_pos, self.dir_path, 'bin', 'gen_cardhtml-pl')],
                        stdout = devnull)

        #return to the initial dir
        os.chdir(old_pos)               

#===============================================================================
# FortranUFOHelasCallWriter
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

        if isinstance(argument, color_ordered_amplitudes.BGHelasCurrent):
            # Create call for wavefunction summation sumVN(fact1,W1,fact2,W2,...,Wres)
            call += "sum%s%d(" % \
                    (color_ordered_amplitudes.spin_dict[\
                                    argument.get('mothers')[0].get('spin')],
                                    len(argument.get('mothers')))
            call += "%s,W(1,%d)," * len(argument.get('mothers')) + \
                    "W(1,%d))"
            call_function = lambda wf: call % \
                (tuple(sum([[self.write_factor(mother.get('factor')),
                             mother.get('number')] for \
                            mother in wf.get('mothers')], []) + \
                [wf.get('number')]))
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
                                 wf.get('number'))
            elif argument.is_boson():
                call_function = lambda wf: call % \
                                (wf.get('number_external'),
                                 wf.get('mass'),
                                 wf.get('number_external'),
                                 1,
                                 wf.get('number_external'),
                                 wf.get('number'))
            else:
                call_function = lambda wf: call % \
                                (wf.get('number_external'),
                                 wf.get('mass'),
                                 wf.get('number_external'),
                                 wf.get('fermionflow'),
                                 wf.get('number_external'),
                                 wf.get('number'))
        # Add the constructed function to wavefunction or amplitude dictionary
            self.add_wavefunction(argument.get_call_key(), call_function)
            return

        # By default, call mother function
        super(COFortranUFOHelasCallWriter, self).generate_helas_call(argument)
