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
from functools import wraps
from nose.tools.nontrivial import nottest
from oncotator.utils.MultiprocessingUtils import MyManager
from ConfigParser import SafeConfigParser
import os
from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.datasources import Gaf
from oncotator.datasources import dbSNP
from oncotator.utils import ConfigUtils
import logging


def data_provider_decorator(fn_data_provider):
    """Data provider decorator, allows another callable to provide the data for the test.
    Modified from https://pypi.python.org/pypi/unittest-data-provider/1.0.0
    to work with nose and to accumulate assertion errors."""

    def test_decorator(fn):
        @wraps(fn)
        def repl(self, *args):
            assertion_errors = []
            ctr = 0
            for i in fn_data_provider():
                try:
                    ctr += 1
                    fn(self, *i)
                except AssertionError as ae:
                    assertion_errors.append("Assertion error data set %s -- %s" % (str(i), ae.message))
            if len(assertion_errors) > 0:
                raise AssertionError("\n"+"\n".join(assertion_errors) + "\n" + str(len(assertion_errors)) + " of " + str(ctr) + " tests failed.")
        return repl
    return test_decorator


class TestUtils(object):
    """
    Stateless class implementing methods that are useful in the unit tests.
    """


    def __init__(self,params):
        """
        Nothing to do when initializing this class.  
        """
        pass
    
    @staticmethod
    def createUnitTestConfig():
        config = SafeConfigParser()
        config.readfp(file('configs/default-test.config','r'))
        config.read(['configs/personal-test.config'])
        return config
    
    @staticmethod
    def createGafDatasource(config, tx_mode="CANONICAL", protocol="file"):
        """ Creates a Gaf 3.0 datasource from a config file.
            Assumes a gaf3.0 section with keys: gaf_fname and gaf_transcript_seqs_fname
            """
        gaf_fname = config.get("gaf3.0", "gaf_fname")
        gaf_transcripts_fname = config.get("gaf3.0", "gaf_transcript_seqs_fname")
        gafDatasource = Gaf(gaf_fname, gaf_transcripts_fname,tx_mode=tx_mode, protocol=protocol)
        return gafDatasource

    @staticmethod
    def createGafDatasourceProxy(config, tx_mode="CANONICAL", protocol="file"):
        """ Creates a Gaf 3.0 datasource from a config file and make a proxy.
            Assumes a gaf3.0 section with keys: gaf_fname and gaf_transcript_seqs_fname
            """
        MyManager.register('Gaf', Gaf)
        manager = MyManager()
        manager.start()
        gaf_fname = config.get("gaf3.0", "gaf_fname")
        gaf_transcripts_fname = config.get("gaf3.0", "gaf_transcript_seqs_fname")
        gafDatasource = manager.Gaf(gaf_fname, gaf_transcripts_fname, tx_mode=tx_mode, protocol="file")
        return gafDatasource

    @staticmethod
    def createDbSnpDatasource(config):
        dbsnpFilename = config.get("dbSNP","dbSNPFilename")
        return dbSNP(dbsnpFilename)

    @staticmethod
    def createCosmicDatasource(config):
        """ Creates a Cosmic datasource from a config file.
            """
        cosmic_dirname = config.get("COSMIC", "CosmicDir")
        cosmicDatasource = DatasourceCreator.createDatasource(cosmic_dirname + "/cosmic.config", cosmic_dirname)
        return cosmicDatasource

    @staticmethod
    def setupLogging(filename, package_name):
        """Static method to set up logging for a unit test.  Argument is usually __file__ from within the unit
        test file.  Should be called globally in a unit test.
        Usually called as follows::
                ...
                        TestUtils.setupLogging(__file__, __name__)
                        class MyClassTest(unittest.TestCase):
                ...

        """
        #TODO: Low priority: Use decorator instead of TestUtils.setupLogging(...)
        curdir = os.path.dirname(filename) + '/'
        logging.basicConfig(filemode='w', filename=(curdir + 'out/oncotator_unitTest_' + package_name + '.log'),
                            level=logging.DEBUG, format='%(asctime)s %(levelname)s [%(name)s:%(lineno)d]  %(message)s')