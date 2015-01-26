#!/usr/local/bin/python2.7
# encoding: utf-8
"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""
from oncotator.utils.RunSpecification import RunSpecification
from oncotator.utils.RunSpecificationFactory import RunSpecificationFactory

'''
Oncotator -- An annotation engine for Cancer

Oncotator is a description

It defines classes_and_methods

@author:     aramos, mgupta, and lichtens
        
@copyright:  2012 Broad Institute. All rights reserved.
        
@license:    TODO: license

@contact:    oncotator@broadinstitute.org
@deffield    updated: Updated
'''
import sys
from oncotator import NGSLIB_INSTALLED
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.utils.MutUtils import MutUtils

if not (sys.version_info[0] == 2  and sys.version_info[1] in [ 7]):
    raise "Oncotator requires Python 2.7.x : " + str(sys.version_info)

import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import logging
from oncotator.utils.version import VERSION 
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils
from oncotator.utils.OptionConstants import OptionConstants

__version__ = VERSION
__all__ = []

__date__ = '2012-12-29'
__updated__ = '2012-12-29'

DEBUG = 1
TESTRUN = 0
PROFILE = 1

#TODO: These need to be dynamic from a config file.
DEFAULT_DB_DIR = '/xchip/cga/reference/annotation/db/oncotator_v1_ds_gencode_current/'
DEFAULT_DEFAULT_ANNOTATIONS = '/xchip/cga/reference/annotation/db/tcgaMAFManualOverrides2.4.config'
DEFAULT_TX_MODE = TranscriptProvider.TX_MODE_CANONICAL


class CLIError(Exception):
    """Generic exception to raise and log different fatal errors."""

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


def parseOptions(program_version_message):
    # Setup argument parser
    description = program_version_message + '''

    Oncotator is a tool for annotating human genomic point mutations and indels with data relevant to cancer researchers.

    '''
    epilog = '''
    Example usage
    -------------
    oncotator -v --input_format=MAFLITE --output_format=TCGAMAF myInputFile.maflite myOutputFile.maf.annotated hg19
    
    IMPORTANT NOTE:  hg19 is only supported genome build for now.

    Default values specified by -d or --default_annotation_values are used when an annotation does not exist or is populated with an empty string ("")

    Both default and override config files and command line specifications stack.

    Example of an override_config or default_config file:

    # Create center, source, sequencer, and score annotations, with the values broad.mit.edu, WXS, Illumina GAIIx, and <blank> for all mutations.
    #  This will overwrite all mutations.
    [manual_annotations]
    override:center=broad.mit.edu,source=WXS,sequencer=Illumina GAIIx,score=

    Example of cache urls:

    # Use a file (/home/user/myfile.cache) ... note the three forward slashes after "file:" for absolute path.
    -u file:///home/user/myfile.cache
    -u file://relative_file.cache

    # memcache
    -u memcache://localhost:11211

    Please note that only VCF input will populate the alt_allele_seen annotation.  All other inputs assume that the alternate is present if it appears at all.
        This feature is to allow users to include or exclude GT of 0/0 or ./. variants when converting VCFs to MAF.

        If --skip-no-alt is specified, VCF input processing will remove mutations with alt_allele_seen of False entirely (the mutations will not even seen when output format is SIMPLE_TSV).

    -----
    Copyright 2012 Broad Institute. All rights reserved.  Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.
    Oncotator is free for non-profit use.  See LICENSE for complete licensing information.
    '''
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: 5]", default=5)
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-i', '--input_format', type=str, default="MAFLITE", choices=OncotatorCLIUtils.getSupportedInputFormats(), help='Input format.  Note that MAFLITE will work for any tsv file with appropriate headers, so long as all of the required headers (or an alias -- see maflite.config) are present.  [default: %s]' % "MAFLITE")
    parser.add_argument('--db-dir', dest='dbDir', default=DEFAULT_DB_DIR,
                        help='Main annotation database directory. [default: %s]' % DEFAULT_DB_DIR)
    parser.add_argument('-o', '--output_format', type=str, default="TCGAMAF",choices=OncotatorCLIUtils.getSupportedOutputFormats(), help='Output format. [default: %s]' % "TCGAMAF")
    parser.add_argument('--override_config', type=str, 
                        help="File path to manual annotations in a config file format (section is 'manual_annotations' and annotation:value pairs).")
    parser.add_argument('--default_config', type=str,
                        help="File path to default annotation values in a config file format (section is 'manual_annotations' and annotation:value pairs).")
    parser.add_argument('--no-multicore', dest="noMulticore", action='store_true', default=False, help="Disables all multicore functionality.")
    parser.add_argument('input_file', type=str, help='Input file to be annotated.  Type is specified through options.')
    parser.add_argument('output_file', type=str, help='Output file name of annotated file.')
    parser.add_argument('genome_build', metavar='genome_build', type=str, help="Genome build.  For example: hg19", choices=["hg19"])
    parser.add_argument('-a', '--annotate-manual', dest="override_cli",type=str, action='append', default=[], help="Specify annotations to override.  Can be specified multiple times.  E.g. -a 'name1:value1' -a 'name2:value2' ")
    parser.add_argument('-d', '--annotate-default', dest="default_cli",type=str, action='append', default=[], help="Specify default values for annotations.  Can be specified multiple times.  E.g. -d 'name1:value1' -d 'name2:value2' ")
    parser.add_argument('-u', '--cache-url', dest="cache_url", type=str, default=None, help=" URL to use for cache.  See help for examples.")
    parser.add_argument('-r', '--read_only_cache', action='store_true', dest="read_only_cache", default=False, help="Makes the cache read-only")
    parser.add_argument('--tx-mode', dest="tx_mode", default=DEFAULT_TX_MODE, choices=TranscriptProvider.TX_MODE_CHOICES, help="Specify transcript mode for transcript providing datasources that support multiple modes.  [default: %s]" % DEFAULT_TX_MODE)
    parser.add_argument('--infer_genotypes', dest='infer_genotypes', default="false", choices=["yes", "true", "t", "1", "y", "no", "false", "f", "0", "n"],
                        help="Forces the VCF output renderer to populate the output genotypes as heterozygous.  This option should only be used when converting a MAFLITE to a VCF; otherwise, the option has no effect.  [default: %s]" % "false")
    parser.add_argument('--skip-no-alt', dest="skip_no_alt", action='store_true', help="If specified, any mutation with annotation alt_allele_seen of 'False' will not be annotated or rendered.  Do not use if output format is a VCF.  If alt_allele_seen annotation is missing, render the mutation.")
    parser.add_argument('--log_name', dest='log_name', default="oncotator.log", help="Specify log output location.  Default: oncotator.log")
    parser.add_argument('--prepend', dest="prepend", action='store_true', help="If specified for TCGAMAF output, will put a 'i_' in front of fields that are not directly rendered in Oncotator TCGA MAFs")
    parser.add_argument('--infer-onps', dest="infer_onps", action='store_true', help="Will merge adjacent SNPs,DNPs,TNPs,etc if they are in the same sample.  This assumes that the input file is position sorted.  This may cause problems with VCF -> VCF conversion, and does not guarantee input order is maintained.")
    parser.add_argument('-c', '--canonical-tx-file', dest="canonical_tx_file", type=str, help="Simple text file with list of transcript IDs (one per line) to always select where possible for variants.  Transcript IDs must match the ones used by the transcript provider in your datasource (e.g. gencode ENST00000123456).  If more than one transcript can be selected for a variant, uses the method defined by --tx-mode to break ties.  Using this list means that a transcript will be selected from this list first, possibly superseding a best-effect.  Note that transcript version number is not considered, whether included in the list or not.")

    # Process arguments
    args = parser.parse_args()
    
    return args


def main(argv=None):  # IGNORE:C0111
    """Command line options."""
    from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils
    from oncotator.Annotator import Annotator

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_version = "%s" % __version__
    program_version_message = '%%(prog)s %s' % program_version

    try:
        args = parseOptions(program_version_message)
        verbose = args.verbose
        if verbose > 0:
            print("Verbose mode on")
        
        logFilename = args.log_name  # 'oncotator.log'

        # Create a basic logger to a file
        loggingFormat = '%(asctime)s %(levelname)s [%(name)s:%(lineno)d] %(message)s'
        logging.basicConfig(filename=logFilename, level=logging.INFO, format=loggingFormat)
        
        
        # Add a console logger to the root logger, which means that all loggers generated will have the console dump.  
        #    Output on the console will be the same as what is in the log file. 
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARN)
        formatter = logging.Formatter(loggingFormat)
        ch.setFormatter(formatter)
        
        if verbose:
            ch.setLevel(logging.INFO)
            print("Path:")
            print(sys.path)
            print(" ")
        
        logging.getLogger('').addHandler(ch)
        
        logger = logging.getLogger(__name__)
        logger.info("Oncotator " + program_version)
        logger.info("Args: " + str(args))
        logger.info('Log file: ' + os.path.abspath(logFilename))
        
        if DEBUG:
            logger.setLevel(logging.DEBUG)

        if not NGSLIB_INSTALLED:
            logger.warn("ngslib module not installed.  Will be unable to annotate with BigWig datasources.")
        
        # Initiate an Oncotator session.
        inputFilename = os.path.expanduser(args.input_file)
        outputFilename = os.path.expanduser(args.output_file)
        inputFormat = args.input_format.upper()
        outputFormat = args.output_format.upper()

        datasourceDir = os.path.expanduser(args.dbDir)
        cache_url = args.cache_url
        read_only_cache = args.read_only_cache
        tx_mode = args.tx_mode
        is_skip_no_alts = args.skip_no_alt
        genome_build = args.genome_build
        is_no_prepend = not args.prepend

        # Parse annotation overrides
        commandLineManualOverrides = args.override_cli
        overrideConfigFile = args.override_config
        if overrideConfigFile is not None and not os.path.exists(overrideConfigFile):
            logger.warn("Could not find " + overrideConfigFile + "   ... proceeding anyway.")
            overrideConfigFile = None
        manualOverrides = OncotatorCLIUtils.determineAllAnnotationValues(commandLineManualOverrides, overrideConfigFile)

        # Parse default overrides
        commandLineDefaultValues = args.default_cli
        defaultConfigFile = args.default_config
        if defaultConfigFile is not None and not os.path.exists(defaultConfigFile):
            if defaultConfigFile != DEFAULT_DEFAULT_ANNOTATIONS:
                logger.warn("Could not find " + defaultConfigFile + "   ... proceeding anyway.")
            else:
                logger.info("Could not find Broad-specific " + defaultConfigFile + "   ... proceeding without any default annotations.  __UNKNOWN__ may appear in TCGA MAF outputs.")
            defaultConfigFile = None
        defaultValues = OncotatorCLIUtils.determineAllAnnotationValues(commandLineDefaultValues, defaultConfigFile)

        # Create a run configuration to pass to the Annotator class.
        annotating_type = None
        if inputFormat == "SEG_FILE":
            annotating_type = RunSpecification.ANNOTATE_SEGMENTS
        runConfig = RunSpecificationFactory.create_run_spec(inputFormat, outputFormat, inputFilename, outputFilename,
                                                      globalAnnotations=manualOverrides, datasourceDir=datasourceDir,
                                                      isMulticore=(not args.noMulticore),
                                                      defaultAnnotations=defaultValues, cacheUrl=cache_url,
                                                      read_only_cache=read_only_cache, tx_mode=tx_mode,
                                                      is_skip_no_alts=is_skip_no_alts, genomeBuild=genome_build,
                                                      other_opts=determineOtherOptions(args), annotating_type=annotating_type)

        annotator = Annotator()
        annotator.initialize(runConfig)
        annotator.annotate()
        
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0


def determineOtherOptions(args):
    opts = dict()
    opts[OptionConstants.NO_PREPEND] = not args.prepend
    opts[OptionConstants.VCF_OUT_INFER_GENOTYPES] = MutUtils.str2bool(args.infer_genotypes)
    opts[OptionConstants.INFER_ONPS] = args.infer_onps
    opts[OptionConstants.CUSTOM_CANONICAL_TX_LIST_FILE] = args.canonical_tx_file
    return opts


def main_profile():
    import cProfile
    import pstats
    
    print("Profiling enabled...")
    profile_filename = 'Oncotator_profile.bin'
    cProfile.run('main()', profile_filename)
    statsfile = open("profile_stats.txt", "wb")
    p = pstats.Stats(profile_filename, stream=statsfile)
    stats = p.strip_dirs().sort_stats('cumulative')
    stats.print_stats()
    statsfile.close()
    sys.exit(0)    

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
        
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        main_profile()
    #sys.exit(main())
    main()
    
