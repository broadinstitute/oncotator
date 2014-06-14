# LICENSE_GOES_HERE

import logging
import re
from shove import Shove
from oncotator.MissingAnnotationException import MissingAnnotationException
from oncotator.datasources.PositionTransformingDatasource import PositionTransformingDatasource


class TranscriptToUniProtProteinPositionTransformingDatasource(PositionTransformingDatasource):
    """ Given a transcript protein sequence, map it to the proper place in the uniprot gene.

    This datasource requires the following annotations to be populated:
        transcript_id
        protein_change

    transcript_id is required as the key into the Shove (sqlite) database.
    protein_change is configurable (i.e. a different annotation can be used), but this is the default.

    src_file is a URL.  Such as sqlite:///db.sqlite

    Annotates gene with [self.title]_aapos

    self.title + "_" will be prepended to the given outputPositionAnnotationName

    Requires a shove database url (as src_file).  The database itself can be created using
        the scripts/uniprot_utils/createUniprotProteinSeqsAlignments.py

    If a transcript is not found in the backing database, the newposition will be ""
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, inputPositionAnnotationName="protein_change", outputPositionAnnotationName="aapos"):
        super(PositionTransformingDatasource, self).__init__(src_file, title=title, version=version)
        self.db = Shove(src_file, "memory://")
        self.inputAnnotationName = inputPositionAnnotationName
        self.outputAnnotationName = self.title + "_" + outputPositionAnnotationName
        self.proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)")

    def _parsePosition(self, position):
        """ Utility class to strip decoration from the position itself.
        """
        tmp = position
        if position.find("p.") != -1:
            tmp = position.replace("p.", "")
        proteinPosition = self.proteinRegexp.match(tmp)
        if proteinPosition is not None:
            return proteinPosition.group(1)
        else:
            return ""

    def _get_uni_pos(self,fh, AA):
        new_pos = 0
        query_AA= ''
        uni_AA = ''
        pat1 = re.compile(r'Query:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
        pat2= re.compile(r'Sbjct:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
        q = 0
        qline = None
        for line in fh:
            if q == 1 and line.startswith('Sbjct: '):
                sline = pat2.search(line)
                new_pos, query_AA, uni_AA = self._map_uni_pos(qline, sline, AA)
                break
            if line.startswith('Query: '):
                qline = pat1.search(line)
                #print '{0}\t{1}'.format(qline.group(1), qline.group(3))
                if int(qline.group(1)) <= AA and int(qline.group(3)) >= AA:
                    q = 1
        return new_pos, query_AA, uni_AA

    def _map_uni_pos(self,qobj, sobj, aa_pos):
        qoff = 0
        soff = 0
        for i in range(len(qobj.group(2))):
            if qobj.group(2)[i] == '-':
                qoff -= 1
            qpos = i + int(qobj.group(1)) + qoff

            if sobj.group(2)[i] == '-':
                soff -= 1
            spos = i + int(sobj.group(1)) + soff

            if qpos == aa_pos:
                return spos, qobj.group(2)[i] , sobj.group(2)[i]

    def annotate_mutation(self, mutation):
        requiredAnnotations = [self.inputAnnotationName, 'transcript_id']

        # Check that all required annotations are present.
        if not all(field in mutation for field in requiredAnnotations):
            missingList = []
            for r in requiredAnnotations:
                if r not in mutation:
                    missingList.append(r)
            raise MissingAnnotationException("Missing required annotation: " + str(missingList))

        val = mutation[self.inputAnnotationName]
        transcript_id = mutation['transcript_id']
        positionOnly = self._parsePosition(val)
        newPos = ""
        if (positionOnly is not None) and (positionOnly != ""):
            try:
                adata = self.db[transcript_id]
                newPos = str(self._get_uni_pos(adata, int(positionOnly))[0])
            except KeyError:
                # Do nothing, we will end up just annotting with ""
                pass

        mutation.createAnnotation(self.outputAnnotationName, newPos, annotationSource=self.title)
        return mutation