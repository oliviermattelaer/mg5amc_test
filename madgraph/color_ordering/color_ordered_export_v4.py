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
        """Export a matrix element to a matrix.f file in MG4 standalone format"""

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
        """Export a matrix element to a matrix.f file in MG4 standalone format"""

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

        # Extract permutation data lines
        perms_data_lines, nperms = self.get_perms_data_lines(matrix_element)
        replace_dict['nperms'] = nperms
        replace_dict['perms_data_lines'] = "\n".join(perms_data_lines)
        
        # Extract nflows
        nflows = len(matrix_element.get('color_flows')) * nperms
        replace_dict['nflows'] = nflows

        # Extract flow function definition lines
        flow_functions_lines = self.get_flow_functions_lines(matrix_element)
        replace_dict['flow_functions_lines'] = flow_functions_lines

        # Extract flow call lines
        flow_call_lines = self.get_flow_call_lines(matrix_element)
        replace_dict['flow_call_lines'] = '\n'.join(flow_call_lines)

        # Extract color sum lines
        color_sum_lines = self.get_color_sum_lines(matrix_element)
        replace_dict['color_sum_lines'] = '\n'.join(color_sum_lines)

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
        """Write out the calls to all color flows"""

        return_lines = []

        nperms = len(matrix_element.get('permutations'))
        for iperm in range(nperms):
            for iflow, flow in enumerate(matrix_element.get('color_flows')):
                return_lines.append("JAMP(%d) = FLOW%d(P,NHEL,PERMS(1,%d))" \
                           % (1 + iperm + iflow * nperms, iflow + 1, iperm + 1))

        return return_lines
    def get_color_sum_lines(self, matrix_element):
        """Write out the data lines defining the permutations"""

        # First calculate the color matrix, since we want to use the full
        # color matrix for the standalone version
        col_mat = color_amp.ColorMatrix(matrix_element.get('color_basis'),
                                        Nc_power_min = \
                                        matrix_element.get('min_Nc_power'))
        color_basis = matrix_element.get('color_basis')
        denoms = col_mat.get_line_denominators()
        flows = matrix_element.get('color_flows')
        nperms = len(matrix_element.get('permutations'))
        res_lines = []

        # Now go row by row in the color matrix, and collect all
        # non-zero entries with Nc_power >= min_Nc_power
        for icol, col_basis_elem in enumerate(sorted(color_basis.keys())):
            numerators = col_mat.get_line_numerators(icol, denoms[icol])
            if not any([n for n in numerators]):
                continue
            for diag_tuple in color_basis[col_basis_elem]:
                iflow = diag_tuple[0]
                iperm = diag_tuple[1][0]
                nflow = 1 + iflow*nperms + iperm
                Nc_power = flows[iflow].get('color_string').Nc_power
                allnums = []
                # Find the factors for all JAMPs that should be included in
                # this multiplication
                for ic, cbe in enumerate(sorted(color_basis.keys())):
                    for dt in color_basis[cbe]:
                        tflow = 1 + dt[0]*nperms + dt[1][0]
                        if numerators[ic]:
                            # Only add if Nc_power large enough
                            Nc_pow = self.Nc_power_to_fraction(\
                                flows[dt[0]].get('color_string').Nc_power)
                            allnums.append((self.fraction_to_string(\
                                                 numerators[ic] * Nc_pow),
                                            tflow))
                if not allnums:
                    continue
                Nc_power = self.Nc_power_to_fraction(Nc_power)
                res_lines.append(\
                    'ZTEMP = ZTEMP+%(fact)s*%(flow)s/%(den)s*DCONJG(%(flows)s)' % \
                    {'fact': self.fraction_to_string(Nc_power),
                     'flow': 'JAMP(%d)' % nflow,
                     'den': str(denoms[icol]),
                     'flows': "+".join(['%s*JAMP(%d)' % allnums[i] for i in \
                                        range(len(allnums))])})
                res_lines[-1] = res_lines[-1].replace('+-', '-')
                res_lines[-1] = res_lines[-1].replace('+1D0*', '+')
                res_lines[-1] = res_lines[-1].replace('/1*', '*')

        return res_lines

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
            # Create call for wavefunction
            call += "sum%s%d(" % \
                    (self.spin_dict[argument.get('mothers')[0].get('spin')],
                                    len(argument.get('mothers')))
            call += "W(1,%d),%s," * len(argument.get('mothers')) + \
                    "W(1,%d))"
            call_function = lambda wf: call % \
                (tuple(sum([[mother.get('number'),
                             self.write_factor(mother.get('factor'))] for \
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
