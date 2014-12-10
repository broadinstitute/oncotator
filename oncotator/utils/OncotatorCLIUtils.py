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
import os
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.input.InputMutationCreator import InputMutationCreatorOptions
from oncotator.input.OnpCombiner import OnpCombiner
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.output.GeneListOutputRenderer import GeneListOutputRenderer
from oncotator.output.VcfOutputRenderer import VcfOutputRenderer
import logging

from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from ConfigUtils import ConfigUtils

from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.output.SimpleBedOutputRenderer import SimpleBedOutputRenderer
from oncotator.output.TcgaVcfOutputRenderer import TcgaVcfOutputRenderer
from oncotator.utils.OptionConstants import OptionConstants
from oncotator.utils.RunSpecification import RunSpecification


class OncotatorCLIUtils(object):
    """
    Utility methods for implementing a command line interface (or any presentation layer class).
    
    Provides utilities for starting an Oncotator session given a set of parameters.
    
    IMPORTANT: If more input/output formats are needed, simply add to the dicts created in  
        createInputFormatNameToClassDict()
        createOutputFormatNameToClassDict()

        First parameter is the class name and the second is the default config file.  If no config file is needed, then
            just specify empty string.  Constructors of the classes must be able to take a config file, even if it
            is ignored.
        
        That will enable support for a new type.  The CLI drives the available 
            formats from the above methods. key is a unique string (e.g. TCGAMAF)
            and the value is the class name (NOT in quotes).  See the method code.
        
    """
    
    def __init__(self, params):
        """
        Constructor -- Never use this.  All methods should be called from a static context.  Throws an exception.
        """
        raise NotImplementedError('This class should not be instantiated.  All methods are static.')

    @staticmethod
    def createInputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more input formats."""
        return {'MAFLITE': (MafliteInputMutationCreator, 'maflite_input.config'),
                "VCF": (VcfInputMutationCreator, 'vcf.in.config'),
                "SEG_FILE": (MafliteInputMutationCreator, 'seg_file_input.config')}

    @staticmethod
    def createOutputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more output formats."""
        return {'TCGAMAF': (TcgaMafOutputRenderer, "tcgaMAF2.4_output.config"),
                "SIMPLE_TSV": (SimpleOutputRenderer, ""),
                'SIMPLE_BED': (SimpleBedOutputRenderer, ""),
                'TCGAVCF': (TcgaVcfOutputRenderer, "tcgaVCF1.1_output.config"),
                'VCF': (VcfOutputRenderer, "vcf.out.config"),
                'GENE_LIST': (GeneListOutputRenderer, "gene_list_output.config")}

    @staticmethod
    def getSupportedOutputFormats():
        """ Lists the supported output formats """
        tmp = OncotatorCLIUtils.createOutputFormatNameToClassDict()
        return tmp.keys()     
    
    @staticmethod
    def getSupportedInputFormats():
        """ Lists the supported input formats """
        tmp = OncotatorCLIUtils.createInputFormatNameToClassDict()
        return tmp.keys()

    @staticmethod
    def create_input_creator(inputFilename, inputFormat, genome_build="hg19", input_creator_options=None):
        inputCreatorDict = OncotatorCLIUtils.createInputFormatNameToClassDict()
        if inputFormat not in inputCreatorDict.keys():
            raise NotImplementedError("The inputFormat specified: " + inputFormat + " is not supported.")
        else:
            inputConfig = inputCreatorDict[inputFormat][1]
            inputCreator = inputCreatorDict[inputFormat][0](inputFilename, inputConfig, genome_build, input_creator_options)
            if input_creator_options.get(OptionConstants.INFER_ONPS):  #If we're combinging ONPs, wrap the input creater
                inputCreator = OnpCombiner(inputCreator)
        return inputCreator

    @staticmethod
    def create_output_renderer(outputFilename, outputFormat, otherOptions):
        outputRendererDict = OncotatorCLIUtils.createOutputFormatNameToClassDict()
        if outputFormat not in outputRendererDict.keys():
            raise NotImplementedError("The outputFormat specified: " + outputFormat + " is not supported.")
        else:
            outputConfig = outputRendererDict[outputFormat][1]
            outputRenderer = outputRendererDict[outputFormat][0](outputFilename, outputConfig, otherOptions)
        return outputRenderer
    
    @staticmethod
    def createManualAnnotationsGivenConfigFile(configFile):
        """
        Assumes a config file as:
        [manual_annotations]
        # annotation3 has blank ('') value
        override:annotation1=value1,annotation2=value2,annotation3=,annotation4=value4
        
        Returns a dictionary:
        {annotation1:value1,annotation2:value2,annotation3:'',annotation4=value4}
        """
        if (configFile is None) or (configFile.strip() == ''):
            return dict()

        if not os.path.exists(configFile):
            logging.getLogger(__name__).warn("Could not find annotation config file: " + configFile + "  ... Proceeding without it.")
            return dict()

        config = ConfigUtils.createConfigParser(configFile)
        opts = config.get('manual_annotations','override').split(',')
        result = dict()
        for optTmp in opts:
            opt = optTmp.split('=')
            if (len(opt) == 1) or (opt[1] is None):
                opt[1] = '' 
            result[opt[0]] = opt[1]
        return result

    @staticmethod
    def determineAllAnnotationValues(commandLineManualOverrides, overrideConfigFile):
        manualOverrides = OncotatorCLIUtils.createManualAnnotationsGivenConfigFile(overrideConfigFile)
        for clmo in commandLineManualOverrides:
            if clmo.find(":") == -1:
                logging.getLogger(__name__).warn("Could not parse annotation: " + str(clmo) + "   ... skipping")
                continue
            keyval = clmo.split(':', 1)
            manualOverrides[keyval[0]] = keyval[1]

        return manualOverrides