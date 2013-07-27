"""
# By downloading the PROGRAM you agree to the following terms of use:
#
# BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
# FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
#
# This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
# WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
# WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
# NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
#
# 1. DEFINITIONS
# 1.1	"PROGRAM" shall mean copyright in the object code and source code known as Oncotator and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.
#
# 2. LICENSE
# 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
# LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
# The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
# 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
# 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
#
# 3. OWNERSHIP OF INTELLECTUAL PROPERTY
# LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
#
# Copyright 2012 Broad Institute, Inc.
# Notice of attribution:  The Oncotator program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc.
#
# LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
#
# 4. INDEMNIFICATION
# LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
#
# 5. NO REPRESENTATIONS OR WARRANTIES
# THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
# IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
#
# 6. ASSIGNMENT
# This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
#
# 7. MISCELLANEOUS
# 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
# 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
# 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
# 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
# 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
# 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
# 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
#"""
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.output.VcfOutputRenderer import VcfOutputRenderer


"""
Class file for OncotatorCLIUtils


Created on Dec 19, 2012

@author: lichtens


"""

# Step 1
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from ConfigUtils import ConfigUtils

from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.output.SimpleBedOutputRenderer import SimpleBedOutputRenderer
from oncotator.output.TcgaVcfOutputRenderer import TcgaVcfOutputRenderer


class RunSpecification(object):
    """ This class contains a specification for an oncotator run.  
    
    This class is quite simple and simply holds parameters.  Annotator classes will be expected to be able to initialize and annotate
        from an instance of this class.
        
        inputCreator -- instance
        outputCreator -- instance
        manualAnnotations -- A dict of name-value pairs that will be annotated onto all mutations.
            name is the name of the annotation.
            value is the value of the annotation.
            The source will always be listed as MANUAL
            TODO: The Annotator takes care of some of this.  Move the relevant documentation.
        datasources -- A list of datasources (instance of Datasource).
        isMulticore -- use multicore processing, where available (True/False)
        numCores -- number of cores to use if isMulticore is True.  Otherwise, this is ignored.
        """
    def __init__(self):
        self.__inputCreator = None
        self.__outputRenderer = None
        self.__inputFilename = None
        self.__outputFilename = None
        self.__manualAnnotations = None
        self.__datasources = None
        self.__isMulticore = False
        self.__numCores = None
        pass

    def get_is_multicore(self):
        return self.__isMulticore


    def get_num_cores(self):
        return self.__numCores


    def set_is_multicore(self, value):
        self.__isMulticore = value


    def set_num_cores(self, value):
        self.__numCores = value


    def del_is_multicore(self):
        del self.__isMulticore


    def del_num_cores(self):
        del self.__numCores


    def get_datasources(self):
        return self.__datasources


    def set_datasources(self, value):
        self.__datasources = value


    def del_datasources(self):
        del self.__datasources


    def get_input_creator(self):
        return self.__inputCreator


    def get_output_renderer(self):
        return self.__outputRenderer


    def get_manual_annotations(self):
        return self.__manualAnnotations


    def set_input_creator(self, value):
        self.__inputCreator = value


    def set_output_renderer(self, value):
        self.__outputRenderer = value

    def set_manual_annotations(self, value):
        self.__manualAnnotations = value

    def del_input_creator(self):
        del self.__inputCreator

    def del_output_renderer(self):
        del self.__outputRenderer

    def del_manual_annotations(self):
        del self.__manualAnnotations
    
    

    def initialize(self, inputCreator, outputRenderer, manualAnnotations=dict(), datasources=[], isMulticore=False, numCores=4, defaultAnnotations=dict()):
        self.inputCreator = inputCreator
        self.outputRenderer = outputRenderer
        self.manualAnnotations = manualAnnotations
        self.datasources = datasources
        self.isMulticore = isMulticore
        self.numCores = numCores
        self.defaultAnnotations = defaultAnnotations

    
    inputCreator = property(get_input_creator, set_input_creator, del_input_creator, "inputCreator's docstring")
    outputRenderer = property(get_output_renderer, set_output_renderer, del_output_renderer, "outputRenderer's docstring")
    manualAnnotations = property(get_manual_annotations, set_manual_annotations, del_manual_annotations, "manualAnnotations's docstring")
    datasources = property(get_datasources, set_datasources, del_datasources, "datasources's docstring")
    isMulticore = property(get_is_multicore, set_is_multicore, del_is_multicore, "isMulticore's docstring")
    numCores = property(get_num_cores, set_num_cores, del_num_cores, "numCores's docstring")

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
    def createRunConfig(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations=dict(), datasourceDir=None, genomeBuild="hg19", isMulticore=False, numCores = 4, defaultAnnotations=dict()):
        ds = DatasourceCreator.createDatasources(datasourceDir, genomeBuild, isMulticore=isMulticore, numCores=numCores)
        return OncotatorCLIUtils.createRunConfigGivenDatasources(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations, ds, genomeBuild, isMulticore, numCores)
    
    @staticmethod
    def createInputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more input formats."""
        return {'MAFLITE':(MafliteInputMutationCreator, 'maflite_input.config'), "VCF":(VcfInputMutationCreator, 'vcf.in.config')}

    @staticmethod
    def createOutputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more output formats."""
        return {'TCGAMAF':(TcgaMafOutputRenderer, 'tcgaMAF2.4_output.config'),"SIMPLE_TSV":(SimpleOutputRenderer, ''),'SIMPLE_BED':(SimpleBedOutputRenderer, ""),'TCGAVCF':(TcgaVcfOutputRenderer, 'tcgaVCF1.1_output.config'), 'VCF':(VcfOutputRenderer, 'vcf.out.config')}

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
    def createRunConfigGivenDatasources(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations=dict(), datasourceList=[], genomeBuild="hg19", isMulticore=False, numCores=4, defaultAnnotations=dict()):
        """ This is a very simple interface to start an Oncotator session.  As a warning, this interface may notbe supported in future versions.
        
        If datasourceDir is None, then the default location is used.  TODO: Define default location.
        
        IMPORTANT: Current implementation attempts to annotate using a default set of datasources.
        
        TODO: Make sure that this note above is no longer the case.  Current implementation attempts to annotate using a default set of datasources
        TODO: This method may get refactored into a separate class that handles RunConfigutaion objects. 
        """  
        # TODO: Use dependency injection for list of name value pairs?  Otherwise, set it up as an attribute on this class.
        # TODO: Use dependency injection to return instance of the input/output classes
        # TODO: Support more than the default configs.
        # TODO: On error, list the supported formats (both input and output) 
        # TODO: Make sure that we can pass in both a class and a config file, not just a class.
        inputCreator = None
        outputRenderer = None
        
        # Step 2
        inputCreatorDict = OncotatorCLIUtils.createInputFormatNameToClassDict()
        if inputFormat not in inputCreatorDict.keys():
            raise NotImplementedError("The inputFormat specified: " + inputFormat + " is not supported.")
        else:
            inputConfig = inputCreatorDict[inputFormat][1]
            inputCreator = inputCreatorDict[inputFormat][0](inputFilename, inputConfig)

        outputRendererDict = OncotatorCLIUtils.createOutputFormatNameToClassDict()   
        if outputFormat not in outputRendererDict.keys():
            raise NotImplementedError("The outputFormat specified: " + outputFormat + " is not supported.")
        else:
            outputConfig = outputRendererDict[outputFormat][1]
            outputRenderer = outputRendererDict[outputFormat][0](outputFilename, outputConfig)
            
        result = RunSpecification()
        result.initialize(inputCreator, outputRenderer, manualAnnotations=globalAnnotations, datasources=datasourceList, isMulticore=isMulticore, numCores=numCores, defaultAnnotations=defaultAnnotations)
        return result
    
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
         
        config = ConfigUtils.createConfigParser(configFile)
        opts = config.get('manual_annotations','override').split(',')
        result = dict()
        for optTmp in opts:
            opt = optTmp.split('=')
            if (len(opt) == 1) or (opt[1] is None):
                opt[1] = '' 
            result[opt[0]] = opt[1]
        return result
        