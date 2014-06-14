# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Oct 23, 2012

@author: gavlee
'''
import unittest
from oncotator.MutationData import MutationData
import logging
from oncotator.DuplicateAnnotationException import DuplicateAnnotationException

TestUtils.setupLogging(__file__, __name__)
class MutationDataTest(unittest.TestCase):
    

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        pass


    def tearDown(self):
        pass


    def testSetValues(self):
        m = MutationData()
        m.createAnnotation("fake1", "1")
        m.createAnnotation("fake2", "blah blah")
        self.assertTrue(m["fake1"] == "1", "Could not properly retrieve annotation using the dictionary interface.  " + str(m["fake1"]))
        self.assertTrue(m["fake2"] == "blah blah", "Could not properly retrieve annotation using the dictionary interface.  " + str(m["fake2"]))
        
        m["fake2"] = "Whoa"
        self.assertTrue(m["fake2"] == "Whoa", "Could not properly retrieve annotation using the dictionary interface, after a value change.")
        print(str(m))
        
    def testIter(self):
        m = MutationData()
        m.createAnnotation("fake1", "1")
        m.createAnnotation("fake2", "blah blah")
        for k in m:
            self.assertTrue((k in ["fake1", "fake2"]) or (k in MutationData.attributes), "Key not present: " + k)

    def testKeys(self):
        m = MutationData()
        m.createAnnotation("fake1", "1")
        m.createAnnotation("fake2", "blah blah")
        self.assertTrue("fake1" in m.keys() and "fake2" in m.keys(), "Could not find annotations added.")

    def testDuplicateException(self):
        ''' Check that a Duplicate Exception is raised by default when annotation value is changed through createAnnotation'''
        m = MutationData()
        m.createAnnotation("fake1", "1")
        with self.assertRaises(DuplicateAnnotationException):
            m.createAnnotation("fake1", "blah blah")
        
    def testDuplicateAnnotationWithSameValue(self):
        ''' Allow a duplicate annotation to be created if the value is the same.  No exception should be thrown.'''
        m = MutationData()
        m.createAnnotation("fake1", "1")
        m.createAnnotation("fake1", "1")
        
    def testAddTag(self):
        ''' Test adding a tag to an annotation '''
        m = MutationData()
        m.createAnnotation("fake1", "1")
        m.addTagToAnnotation("fake1", "fakeTag")
        self.assertTrue("fakeTag" in m.getAnnotation("fake1").getTags(), "Tag was not added properly.")

    def testPickleable(self):
        """Test that a near-empty MutationData can be pickled"""
        m = MutationData()
        m.chr = "2"
        m.createAnnotation("fake1", "1")
        m.addTagToAnnotation("fake1", "fakeTag")
        import cPickle
        cPickle.dump(m, open("out/testMDPickle.pkl", 'w'))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSetValues']
    unittest.main()