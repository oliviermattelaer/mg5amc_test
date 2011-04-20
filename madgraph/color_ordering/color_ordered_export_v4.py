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
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.drawing_eps as draw
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.files as files
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.misc as misc
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.gen_infohtml as gen_infohtml
import madgraph.iolibs.template_files as iolibs_template_files
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
class ProcessExporterFortranCO(ProcessExporterFortranME):
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

    #===========================================================================
    # export the helas routine
    #===========================================================================
    def export_helas(self, helas_path):
        """Configure the files/link of the process according to the model"""

        # Import helas routine
        for filename in os.listdir(helas_path):
            filepos = os.path.join(helas_path, filename)
            if os.path.isfile(filepos):
                if filepos.endswith('Makefile.template'):
                    cp(filepos, self.dir_path + '/Source/DHELAS/Makefile')
                elif filepos.endswith('Makefile'):
                    pass
                else:
                    cp(filepos, self.dir_path + '/Source/DHELAS')

    #===========================================================================
    # write_matrix_element_v4
    #===========================================================================
    def write_matrix_element_v4(self, writer, matrix_element, fortran_model,
                                proc_id = "", config_map = []):
        """Export a matrix element to a matrix.f file in MG4 madevent format"""

        if not matrix_element.get('processes') or \
               not matrix_element.get('diagrams'):
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

        # Set proc_id
        replace_dict['proc_id'] = proc_id

        # Extract ncomb
        ncomb = matrix_element.get_helicity_combinations()
        replace_dict['ncomb'] = ncomb

        # Extract helicity lines
        helicity_lines = self.get_helicity_lines(matrix_element)
        replace_dict['helicity_lines'] = helicity_lines

        # Extract IC line
        ic_line = self.get_ic_line(matrix_element)
        replace_dict['ic_line'] = ic_line

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

        # Extract ncolor
        ncolor = max(1, len(matrix_element.get('color_basis')))
        replace_dict['ncolor'] = ncolor

        # Extract color data lines
        color_data_lines = self.get_color_data_lines(matrix_element)
        replace_dict['color_data_lines'] = "\n".join(color_data_lines)

        # Extract helas calls
        helas_calls = fortran_model.get_matrix_element_calls(\
                    matrix_element)
        replace_dict['helas_calls'] = "\n".join(helas_calls)

        # Extract amp2 lines
        amp2_lines = self.get_amp2_lines(matrix_element, config_map)
        replace_dict['amp2_lines'] = '\n'.join(amp2_lines)

        # Extract JAMP lines
        jamp_lines = self.get_JAMP_lines(matrix_element)
        replace_dict['jamp_lines'] = '\n'.join(jamp_lines)

        file = open(os.path.join(_file_path, \
                          'iolibs/template_files/%s' % self.matrix_file)).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return len(filter(lambda call: call.find('#') != 0, helas_calls)), ncolor

    #===========================================================================
    # write_auto_dsig_file
    #===========================================================================
    def write_auto_dsig_file(self, writer, matrix_element, proc_id = ""):
        """Write the auto_dsig.f file for the differential cross section
        calculation, includes pdf call information"""

        if not matrix_element.get('processes') or \
               not matrix_element.get('diagrams'):
            return 0

        nexternal, ninitial = matrix_element.get_nexternal_ninitial()

        if ninitial < 1 or ninitial > 2:
            raise writers.FortranWriter.FortranWriterError, \
                  """Need ninitial = 1 or 2 to write auto_dsig file"""

        replace_dict = {}

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(matrix_element)
        replace_dict['process_lines'] = process_lines

        # Set proc_id
        replace_dict['proc_id'] = proc_id
        replace_dict['numproc'] = 1

        # Set dsig_line
        if ninitial == 1:
            # No conversion, since result of decay should be given in GeV
            dsig_line = "pd(IPROC)*dsiguu"
        else:
            # Convert result (in GeV) to pb
            dsig_line = "pd(IPROC)*conv*dsiguu"

        replace_dict['dsig_line'] = dsig_line

        # Extract pdf lines
        pdf_lines = self.get_pdf_lines(matrix_element, ninitial, proc_id != "")
        replace_dict['pdf_lines'] = pdf_lines

        # Lines that differ between subprocess group and regular
        if proc_id:
            replace_dict['numproc'] = int(proc_id)
            replace_dict['passcuts_begin'] = ""
            replace_dict['passcuts_end'] = ""
            # Set lines for subprocess group version
            # Set define_iconfigs_lines
            replace_dict['define_subdiag_lines'] = \
                 """\nINTEGER SUBDIAG(MAXSPROC),IB(2)
                 COMMON/TO_SUB_DIAG/SUBDIAG,IB"""    
        else:
            replace_dict['passcuts_begin'] = "IF (PASSCUTS(PP)) THEN"
            replace_dict['passcuts_end'] = "ENDIF"
            replace_dict['define_subdiag_lines'] = ""

        file = open(os.path.join(_file_path, \
                          'iolibs/template_files/auto_dsig_v4.inc')).read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

    #===========================================================================
    # write_coloramps_file
    #===========================================================================
    def write_coloramps_file(self, writer, mapconfigs, matrix_element):
        """Write the coloramps.inc file for MadEvent"""

        lines = self.get_icolamp_lines(mapconfigs, matrix_element, 1)
        lines.insert(0, "logical icolamp(%d,%d,1)" % \
                        (max(len(matrix_element.get('color_basis').keys()), 1),
                         len(mapconfigs)))


        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_maxconfigs_file
    #===========================================================================
    def write_maxconfigs_file(self, writer, matrix_elements):
        """Write the maxconfigs.inc file for MadEvent"""

        if isinstance(matrix_elements, helas_objects.HelasMultiProcess):
            maxconfigs = max([me.get_num_configs() for me in \
                              matrix_elements.get('matrix_elements')])
        else:
            maxconfigs = max([me.get_num_configs() for me in matrix_elements])

        lines = "integer lmaxconfigs\n"
        lines += "parameter(lmaxconfigs=%d)" % maxconfigs

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_configs_file
    #===========================================================================
    def write_configs_file(self, writer, matrix_element):
        """Write the configs.inc file for MadEvent"""

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()

        configs = [(i+1, d) for i,d in enumerate(matrix_element.get('diagrams'))]
        mapconfigs = [c[0] for c in configs]
        return mapconfigs, self.write_configs_file_from_diagrams(writer,
                                                            [c[1] for c in configs],
                                                            mapconfigs,
                                                            nexternal, ninitial)

    #===========================================================================
    # write_configs_file_from_diagrams
    #===========================================================================
    def write_configs_file_from_diagrams(self, writer, configs, mapconfigs,
                                         nexternal, ninitial):
        """Write the actual configs.inc file.
        configs is the diagrams corresponding to configs,
        mapconfigs gives the diagram number for each config."""

        lines = []

        s_and_t_channels = []

        minvert = min([max(diag.get_vertex_leg_numbers()) for diag in configs])

        nconfigs = 0

        for iconfig, helas_diag in enumerate(configs):
            if any([vert > minvert for vert in
                    helas_diag.get_vertex_leg_numbers()]):
                # Only 3-vertices allowed in configs.inc
                continue
            nconfigs += 1

            # Need to reorganize the topology so that we start with all
            # final state external particles and work our way inwards

            schannels, tchannels = helas_diag.get('amplitudes')[0].\
                                              get_s_and_t_channels(ninitial)

            s_and_t_channels.append([schannels, tchannels])

            allchannels = schannels
            if len(tchannels) > 1:
                # Write out tchannels only if there are any non-trivial ones
                allchannels = schannels + tchannels

            # Write out propagators for s-channel and t-channel vertices

            lines.append("# Diagram %d" % (mapconfigs[iconfig]))
            # Correspondance between the config and the diagram = amp2
            lines.append("data mapconfig(%d)/%d/" % (nconfigs,
                                                     mapconfigs[iconfig]))

            for vert in allchannels:
                daughters = [leg.get('number') for leg in vert.get('legs')[:-1]]
                last_leg = vert.get('legs')[-1]
                lines.append("data (iforest(i,%d,%d),i=1,%d)/%s/" % \
                             (last_leg.get('number'), nconfigs, len(daughters),
                              ",".join([str(d) for d in daughters])))
                if vert in schannels:
                    lines.append("data sprop(%d,%d)/%d/" % \
                                 (last_leg.get('number'), nconfigs,
                                  last_leg.get('id')))
                    lines.append("data tprid(%d,%d)/0/" % \
                                 (last_leg.get('number'), nconfigs))
                elif vert in tchannels[:-1]:
                    lines.append("data tprid(%d,%d)/%d/" % \
                                 (last_leg.get('number'), nconfigs,
                                  abs(last_leg.get('id'))))
                    lines.append("data sprop(%d,%d)/0/" % \
                                 (last_leg.get('number'), nconfigs))

        # Write out number of configs
        lines.append("# Number of configs")
        lines.append("data mapconfig(0)/%d/" % nconfigs)

        # Write the file
        writer.writelines(lines)

        return s_and_t_channels

    #===========================================================================
    # write_decayBW_file
    #===========================================================================
    def write_decayBW_file(self, writer, s_and_t_channels):
        """Write the decayBW.inc file for MadEvent"""

        lines = []

        booldict = {False: ".false.", True: ".true."}

        for iconf, config in enumerate(s_and_t_channels):
            schannels = config[0]
            for vertex in schannels:
                # For the resulting leg, pick out whether it comes from
                # decay or not, as given by the from_group flag
                leg = vertex.get('legs')[-1]
                lines.append("data gForceBW(%d,%d)/%s/" % \
                             (leg.get('number'), iconf + 1,
                              booldict[leg.get('from_group')]))

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_dname_file
    #===========================================================================
    def write_dname_file(self, writer, dir_name):
        """Write the dname.mg file for MG4"""

        line = "DIRNAME=%s" % dir_name

        # Write the file
        writer.write(line + "\n")

        return True

    #===========================================================================
    # write_iproc_file
    #===========================================================================
    def write_iproc_file(self, writer, me_number):
        """Write the iproc.dat file for MG4"""
        line = "%d" % (me_number + 1)

        # Write the file
        for line_to_write in writer.write_line(line):
            writer.write(line_to_write)
        return True

    #===========================================================================
    # write_leshouche_file
    #===========================================================================
    def write_leshouche_file(self, writer, matrix_element):
        """Write the leshouche.inc file for MG4"""

        # Write the file
        writer.writelines(self.get_leshouche_lines(matrix_element, 0))

        return True

    #===========================================================================
    # get_leshouche_lines
    #===========================================================================
    def get_leshouche_lines(self, matrix_element, numproc):
        """Write the leshouche.inc file for MG4"""

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()

        lines = []
        for iproc, proc in enumerate(matrix_element.get('processes')):
            legs = proc.get_legs_with_decays()
            lines.append("DATA (IDUP(i,%d,%d),i=1,%d)/%s/" % \
                         (iproc + 1, numproc+1, nexternal,
                          ",".join([str(l.get('id')) for l in legs])))
            if iproc == 0 and numproc == 0:
                for i in [1, 2]:
                    lines.append("DATA (MOTHUP(%d,i),i=1,%2r)/%s/" % \
                             (i, nexternal,
                              ",".join([ "%3r" % 0 ] * ninitial + \
                                       [ "%3r" % i ] * (nexternal - ninitial))))

            # Here goes the color connections corresponding to the JAMPs
            # Only one output, for the first subproc!
            if iproc == 0:
                # If no color basis, just output trivial color flow
                if not matrix_element.get('color_basis'):
                    for i in [1, 2]:
                        lines.append("DATA (ICOLUP(%d,i,1,%d),i=1,%2r)/%s/" % \
                                 (i, numproc+1,nexternal,
                                  ",".join([ "%3r" % 0 ] * nexternal)))

                else:
                    # First build a color representation dictionnary
                    repr_dict = {}
                    for l in legs:
                        repr_dict[l.get('number')] = \
                            proc.get('model').get_particle(l.get('id')).get_color()\
                            * (-1)**(1+l.get('state'))
                    # Get the list of color flows
                    color_flow_list = \
                        matrix_element.get('color_basis').color_flow_decomposition(repr_dict,
                                                                                   ninitial)
                    # And output them properly
                    for cf_i, color_flow_dict in enumerate(color_flow_list):
                        for i in [0, 1]:
                            lines.append("DATA (ICOLUP(%d,i,%d,%d),i=1,%2r)/%s/" % \
                                 (i + 1, cf_i + 1, numproc+1, nexternal,
                                  ",".join(["%3r" % color_flow_dict[l.get('number')][i] \
                                            for l in legs])))

        return lines

    #===========================================================================
    # write_maxamps_file
    #===========================================================================
    def write_maxamps_file(self, writer, maxamps, maxflows,
                           maxproc,maxsproc):
        """Write the maxamps.inc file for MG4."""

        file = "       integer    maxamps, maxflow, maxproc, maxsproc\n"
        file = file + "parameter (maxamps=%d, maxflow=%d)\n" % \
               (maxamps, maxflows)
        file = file + "parameter (maxproc=%d, maxsproc=%d)" % \
               (maxproc, maxsproc)

        # Write the file
        writer.writelines(file)

        return True

    #===========================================================================
    # write_mg_sym_file
    #===========================================================================
    def write_mg_sym_file(self, writer, matrix_element):
        """Write the mg.sym file for MadEvent."""

        lines = []

        # Extract process with all decays included
        final_legs = filter(lambda leg: leg.get('state') == True,
                       matrix_element.get('processes')[0].get_legs_with_decays())

        ninitial = len(filter(lambda leg: leg.get('state') == False,
                              matrix_element.get('processes')[0].get('legs')))

        identical_indices = {}

        # Extract identical particle info
        for i, leg in enumerate(final_legs):
            if leg.get('id') in identical_indices:
                identical_indices[leg.get('id')].append(\
                                    i + ninitial + 1)
            else:
                identical_indices[leg.get('id')] = [i + ninitial + 1]

        # Remove keys which have only one particle
        for key in identical_indices.keys():
            if len(identical_indices[key]) < 2:
                del identical_indices[key]

        # Write mg.sym file
        lines.append(str(len(identical_indices.keys())))
        for key in identical_indices.keys():
            lines.append(str(len(identical_indices[key])))
            for number in identical_indices[key]:
                lines.append(str(number))

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_mg_sym_file
    #===========================================================================
    def write_default_mg_sym_file(self, writer):
        """Write the mg.sym file for MadEvent."""

        lines = "0"

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_ncombs_file
    #===========================================================================
    def write_ncombs_file(self, writer, nexternal):
        """Write the ncombs.inc file for MadEvent."""

        # ncomb (used for clustering) is 2^nexternal
        file = "       integer    n_max_cl\n"
        file = file + "parameter (n_max_cl=%d)" % (2 ** nexternal)

        # Write the file
        writer.writelines(file)

        return True

    #===========================================================================
    # write_props_file
    #===========================================================================
    def write_props_file(self, writer, matrix_element, s_and_t_channels):
        """Write the props.inc file for MadEvent. Needs input from
        write_configs_file."""

        lines = []

        particle_dict = matrix_element.get('processes')[0].get('model').\
                        get('particle_dict')

        for iconf, configs in enumerate(s_and_t_channels):
            for vertex in configs[0] + configs[1][:-1]:
                leg = vertex.get('legs')[-1]
                if leg.get('id') == 21 and 21 not in particle_dict:
                    # Fake propagator used in multiparticle vertices
                    mass = 'zero'
                    width = 'zero'
                    pow_part = 0
                else:
                    particle = particle_dict[leg.get('id')]
                    # Get mass
                    if particle.get('mass').lower() == 'zero':
                        mass = particle.get('mass')
                    else:
                        mass = "abs(%s)" % particle.get('mass')
                    # Get width
                    if particle.get('width').lower() == 'zero':
                        width = particle.get('width')
                    else:
                        width = "abs(%s)" % particle.get('width')

                    pow_part = 1 + int(particle.is_boson())

                lines.append("pmass(%d,%d)  = %s" % \
                             (leg.get('number'), iconf + 1, mass))
                lines.append("pwidth(%d,%d) = %s" % \
                             (leg.get('number'), iconf + 1, width))
                lines.append("pow(%d,%d) = %d" % \
                             (leg.get('number'), iconf + 1, pow_part))

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_processes_file
    #===========================================================================
    def write_processes_file(self, writer, subproc_group):
        """Write the processes.dat file with info about the subprocesses
        in this group."""

        lines = []

        for ime, me in \
            enumerate(subproc_group.get('matrix_elements')):
            lines.append("%s %s" % (str(ime+1) + " " * (7-len(str(ime+1))),
                                    ",".join(p.base_string() for p in \
                                             me.get('processes'))))
            if me.get('has_mirror_process'):
                mirror_procs = [copy.copy(p) for p in me.get('processes')]
                for proc in mirror_procs:
                    legs = copy.copy(proc.get('legs'))
                    legs.insert(0, legs.pop(1))
                    proc.set("legs", legs)
                lines.append("mirror  %s" % ",".join(p.base_string() for p in \
                                                     mirror_procs))
            else:
                lines.append("mirror  none")

        # Write the file
        writer.write("\n".join(lines))

        return True

    #===========================================================================
    # write_symswap_file
    #===========================================================================
    def write_symswap_file(self, writer, ident_perms):
        """Write the file symswap.inc for MG4 by comparing diagrams using
        the internal matrix element value functionality."""

        lines = []

        # Write out lines for symswap.inc file (used to permute the
        # external leg momenta
        for iperm, perm in enumerate(ident_perms):
            lines.append("data (isym(i,%d),i=1,nexternal)/%s/" % \
                         (iperm+1, ",".join([str(i+1) for i in perm])))
        lines.append("data nsym/%d/" % len(ident_perms))

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_symfact_file
    #===========================================================================
    def write_symfact_file(self, writer, symmetry):
        """Write the files symfact.dat for MG4 by comparing diagrams using
        the internal matrix element value functionality."""


        # Write out lines for symswap.inc file (used to permute the
        # external leg momenta
        lines = [ "%3r %3r" %(i+1, s) for i,s in enumerate(symmetry) if s != 0] 

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_symperms_file
    #===========================================================================
    def write_symperms_file(self, writer, perms):
        """Write the symperms.inc file for subprocess group, used for
        symmetric configurations"""

        lines = []
        for iperm, perm in enumerate(perms):
            lines.append("data (perms(i,%d),i=1,nexternal)/%s/" % \
                         (iperm+1, ",".join([str(i+1) for i in perm])))

        # Write the file
        writer.writelines(lines)

        return True

    #===========================================================================
    # write_subproc
    #===========================================================================
    def write_subproc(self, writer, subprocdir):
        """Append this subprocess to the subproc.mg file for MG4"""

        # Write line to file
        writer.write(subprocdir + "\n")

        return True

