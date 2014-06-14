# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 9, 2013

@author: lichtens
'''
import unittest

from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotator import Annotator
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator 
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.MutationData import MutationData
import logging

TestUtils.setupLogging(__file__, __name__)
class GenericGeneDataSourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()


    def tearDown(self):
        pass


        
    def testBasicAnnotation(self):
        ''' Annotate from a basic tsv gene file.  Use the Gaf to annotate before trying the tsv -- required since the gene annotation must be populated.
        Using trimmed CancerGeneCensus as basis for this test.
        ''' 
        
        # cut -f 1 oncotator/test/testdata/small_tsv_ds/CancerGeneCensus_Table_1_full_2012-03-15_trim.txt | egrep -v Symbol | sed -r "s/^/'/g" | sed ':a;N;$!ba;s/\n/,/g' | sed -r "s/,'/','/g"
        genesAvailable = ['ABL1','ABL2','ACSL3','AF15Q14','AF1Q','AF3p21','AF5q31','AKAP9','AKT1','AKT2','ALDH2','ALK','ALO17','APC','ARHGEF12','ARHH','ARID1A','ARID2','ARNT','ASPSCR1','ASXL1','ATF1','ATIC','ATM','ATRX','BAP1','BCL10','BCL11A','BCL11B']
        
        # We need a gaf data source to annotate gene

        gafDatasource = TestUtils.createTranscriptProviderDatasource(config=self.config)
        geneDS = DatasourceFactory.createDatasource("testdata/small_tsv_ds/small_tsv_ds.config", "testdata/small_tsv_ds/")
        outputFilename = 'out/genericGeneTest.out.tsv'
        
        annotator = Annotator()
        annotator.setInputCreator(MafliteInputMutationCreator('testdata/maflite/Patient0.snp.maf.txt'))
        annotator.setOutputRenderer(SimpleOutputRenderer(outputFilename))
        annotator.addDatasource(gafDatasource)
        annotator.addDatasource(geneDS)
        annotator.annotate()
        
        # Check that there were actual annotations performed.
        tsvReader = GenericTsvReader(outputFilename)
        
        fields = tsvReader.getFieldNames()
        self.assertTrue('CGC_Abridged_Other Syndrome/Disease' in fields, "'CGC_Other Syndrome/Disease' was not present in the header")
        self.assertTrue('CGC_Abridged_Mutation Type' in fields, "'CGC_Abridged_Mutation Type' was not present in the header")
        
        ctr = 1
        linesThatShouldBeAnnotated = 0
        for lineDict in tsvReader:
            self.assertTrue('gene' in lineDict.keys())
            if lineDict['gene'] in genesAvailable:
                self.assertTrue(lineDict['CGC_Abridged_GeneID'] <> '', "'CGC_Abridged_GeneID' was missing on a row that should have been populated.  Line: " + str(ctr))
                linesThatShouldBeAnnotated = linesThatShouldBeAnnotated + 1
            ctr = ctr + 1
        self.assertTrue((linesThatShouldBeAnnotated) > 0, "Bad data -- cannot test missed detects.")
    
    def testAnnotationSourceIsPopulated(self):
        ''' Tests that the annotation source is not blank for the example tsv datasource. '''
        geneDS = DatasourceFactory.createDatasource("testdata/small_tsv_ds/small_tsv_ds.config", "testdata/small_tsv_ds/")
        self.assertTrue(geneDS <> None, "gene indexed datasource was None.")
        
        m = MutationData()
        m.createAnnotation('gene',"ABL1")
        m = geneDS.annotate_mutation(m)
        self.assertTrue(m['CGC_Abridged_Name'] == "v-abl Abelson murine leukemia viral oncogene homolog 1","Test gene TSV datasource did not annotate properly.")
        self.assertTrue(m.getAnnotation('CGC_Abridged_Name').getDatasource() <> "Unknown", "Annotation source was unknown")
        self.assertTrue(m.getAnnotation('CGC_Abridged_Name').getDatasource().strip() <> "", "Annotation source was blank")
        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()