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


"""
Created on Nov 9, 2012

@author: lichtens
"""
from ConfigParser import SafeConfigParser
import os.path
import os
import re
import logging
from os.path import expanduser

class ConfigUtils(object):
    """
    Collection of utility methods.
    """


    def __init__(self,params):
        """
        Constructor
        """
        pass
    
    @staticmethod
    def buildAlternateKeyDictionaryFromConfig(configParser, sectionKey="alternatives"):
        """ Creates a dictionary of columnName to list of alternative names.  This reads from the internal configparser.
        Parses alternatives as a comma-separated list"""
        result = dict()
        keys = configParser.options(sectionKey)
        for k in keys:
            result[k] = configParser.get(sectionKey, k).split(',')
        return result
    
    @staticmethod
    def buildReverseAlternativeDictionary(alternativeDict):
        """ Creates a dictionary where the values become the keys and the keys become the values.  Note that this method takes a snapshot at time of creation.
        Important:  This method also flattens lists.  For example:
        {'a':['1','2']} --> {'1':'a', '2','a'}
        
        Important: This method assumes that the values are all lists.
        """
        result = dict()
        for key in alternativeDict.keys():
            alts = alternativeDict[key]
            for a in alts:
                result[a] = key
        return result
    
    @staticmethod
    def buildReverseAlternativeDictionaryFromConfig(configParser, sectionKey="alternatives"):
        """ Creates a dictionary where the values become the keys and the keys become the values.  Note that this method takes a snapshot at time of creation.
        Important:  This method also flattens lists.  For example:
        {'a':['1','2']} --> {'1':'a', '2','a'}
        """
        result = dict()
        keys = configParser.options(sectionKey)
        for k in keys:
            alts = configParser.get(sectionKey, k).split(',')
            for a in alts:
                result[a] = k
        return result
    
    @staticmethod
    def hasSectionKey(configParser, sectionKey=None):
        """ Checks whether the config file has a section with sectionKey name or not. """
        if sectionKey is not None:
            return configParser.has_section(sectionKey)
        return False
    
    @staticmethod
    def createConfigParser(sourceConfigFile, ignoreCase=True, additional_config_dir=""):
        """ Creates a config parser instance using the Oncotator conventions for config files.  
        


        Note:  ALL config files must end with .config
        
        In order of precedence:
            basicConfig file (.config)
            localBasicConfig file (.local.config)
            
            The filename convention for a local basic config is as follows:
                myConfig.config (the basicConfig)
                myConfig.local.config 
                
                All sections, keys, and values in the local config supersede the basic config.  No local configs should be 
                committed to git.
        
        IMPORTANT: If the config file is not found.  Attempts to find it again by prepending "configs/" and then checking
         the additional_config_dir.
        If still not found, assumes that the file is part of the install package. 
           
           
        Precedence: 
                1) File name relative to current dir (CWD)
                2) Filename configs/CWD
                3) Additional directory (additional_config_dir)
                4) Package resources 
                
                .local.config can be found in any of these locations.
                
                In other words, test.config can be package resources and test.local.config can be in ~/.oncotator, 
                    and it will be as if these two files are local to each other.   

        :param sourceConfigFile:
        :param ignoreCase:
        :param additional_config_dir: Lowest precedence directory to search.  Does NOT look for a config subdirectory here


            """
       
        # TODO: Must be more graceful way of handling this.
        # TODO: Add ~/.oncotator to areas that are searched for config files.
        
        sourceConfigFP = ConfigUtils._locateConfigFile(sourceConfigFile, additional_config_dir)
        config = SafeConfigParser()
        if not ignoreCase:
            config.optionxform = str
        config.readfp(sourceConfigFP)
        
        # Create the local filename string by grabbing the ending of .config and 
        configBaseFilename = os.path.basename(sourceConfigFP.name)
        
        pattern = "\.config$"
        configLocalFilename = re.sub(pattern, ".local.config", configBaseFilename)
        
        dirName = os.path.dirname(sourceConfigFile)
        if dirName <> "":
            dirName = dirName + '/'
        configLocalSourceFilename = dirName + configLocalFilename
        configLocalSourceFP = ConfigUtils._locateConfigFile(configLocalSourceFilename, isRelaxedLogging=True, additional_config_dir=additional_config_dir)
        if configLocalSourceFP <> None:
            config.readfp(configLocalSourceFP)
        return config
    
    @staticmethod
    def _locateConfigFile(sourceConfigFile, isRelaxedLogging=False, additional_config_dir=None):
        """ This method implements the searching as described in createConfigParser.
        
        IMPORTANT: This returns a file pointer, not a filename. 
        
        isRelaxedLogging of True disables warning messages of missing files.  This is useful for when looking for config files that are not likely to exist. 
        """
        if not os.path.exists(sourceConfigFile):
            if not isRelaxedLogging:
                logging.getLogger(__name__).debug("Could not find config file (" + sourceConfigFile + ").  Trying configs/ prepend.")
            if os.path.exists("configs/"+sourceConfigFile):
                logging.getLogger(__name__).info("Found config file (" + sourceConfigFile + ") using configs/ prepend.")
                sourceConfigFile = "configs/" + sourceConfigFile
                sourceConfigFP = file(sourceConfigFile, 'r')
            elif additional_config_dir is not None and additional_config_dir.strip() != "" and os.path.exists(additional_config_dir + "/" + sourceConfigFile):
                logging.getLogger(__name__).info("Found config file (" + sourceConfigFile + ") using " + additional_config_dir + ".")
                sourceConfigFile = additional_config_dir + "/" + sourceConfigFile
                sourceConfigFP = file(sourceConfigFile, 'r')
            else:
                logging.getLogger(__name__).debug("Attempting to get " + sourceConfigFile + " from package resources")
                sourceConfigFP = None
                try:
                    from pkg_resources import resource_stream
                    sourceConfigFP = resource_stream("oncotator.configs", sourceConfigFile)
                except IOError as ioe:
                    if not isRelaxedLogging:
                        logging.getLogger(__name__).warn("Could not find " + sourceConfigFile + " in installed resources ")   
        else:
            sourceConfigFP = file(sourceConfigFile,'r')       
        return sourceConfigFP

    @staticmethod
    def createTemplateFP(templateName):
        sourceConfigFP = None
        if not os.path.exists(templateName):
            try:
                from pkg_resources import resource_stream

                sourceConfigFP = resource_stream("oncotator.configs", templateName)
            except IOError as ioe:
                print("Could not find template file in installed resources ")
        else:
            sourceConfigFP = file('ds_config.template', 'r')
        return sourceConfigFP

    @staticmethod
    def hasSectionKeys(configParser, sectionKeys=[]):
        """ Checks whether the config file has a section with sectionKey name or not. """
        for sectionKey in sectionKeys:
            if sectionKey is not None and not configParser.has_section(sectionKey):
                return False
        return True