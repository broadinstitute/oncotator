# LICENSE_GOES_HERE


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
                logging.getLogger(__name__).warn("Could not find config file (" + sourceConfigFile + ").  Trying configs/ prepend.")
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