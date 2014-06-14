# LICENSE_GOES_HERE


import csv
from string import Template
from datetime import datetime
from oncotator.MissingAnnotationException import MissingAnnotationException
from oncotator.output.OutputRenderer import OutputRenderer
from oncotator.utils.ConfigUtils import ConfigUtils
import logging
from oncotator.utils.MutUtils import MutUtils


class TcgaVcfOutputRenderer(OutputRenderer):
    """
    Create a TCGA VCF (version 1.1).

    This implementation is very basic.  Basically, this populates a barebones TCGA VCF, using a template.

    Notes:
    ONLY SUPPORTS HG19.  Makes hg19 assumption in several places.

    The required annotations may change in future versions of this code.
    The following annotations are required:
    center
    vcfProcessLog
    individual_barcode
    normal_barcode
    platform
    source
    normal_accession
    softwareName
    softwareVersion
    softwareParams
    normal_file
    normal_uuid
    tumor_barcode
    tumor_accession
    tumor_file
    tumor_uuid
    tumor_subtype

    """

    requiredHeaderAnnotations = ["center",
                                 "vcfProcessLog",
                                 "individual_barcode",
                                 "normal_barcode",
                                 "platform",
                                 "source",
                                 "normal_accession",
                                 "softwareName",
                                 "softwareVersion",
                                 "softwareParams",
                                 "normal_file",
                                 "normal_uuid",
                                 "tumor_barcode",
                                 "tumor_accession",
                                 "tumor_file",
                                 "tumor_uuid",
                                 "tumor_subtype",
                                 "geneAnno"]

    requiredMutAnnotations = ['t_lod_fstar',
                              'init_n_lod',
                              't_ref_count',
                              't_alt_count',
                              'n_ref_count',
                              'n_alt_count',
                              't_alt_sum',
                              'n_ref_sum',
                              'ref_context',
                              'dbSNP_RS',
                              'judgement']

    def __init__(self, filename, configFile="tcgaVCF1.1_output.config", otherOptions=None):
        self._filename = filename
        self.logger = logging.getLogger(__name__)
        self.config = ConfigUtils.createConfigParser(configFile)
        self.alternativeDictionary = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config)
        self.seenDbSNPs = dict()
        self.fieldMap = {}

    def createVcfHeader(self, m, commentString=""):
        """Create the VCF Header using a simple template. """
        sourceConfigFP = ConfigUtils.createTemplateFP("tcgaVCF1.1Header.template")
        sHeaderTemplate = Template(sourceConfigFP.read())

        missingAnnotations = []
        headerSubKeys = dict()
        for reqHdr in TcgaVcfOutputRenderer.requiredHeaderAnnotations:
            if reqHdr not in m.keys():
                missingAnnotations.append(reqHdr)
                headerSubKeys[reqHdr] = "."
            else:
                headerSubKeys[reqHdr] = m[reqHdr]

        headerSubKeys['date'] = str(datetime.now().date()).replace('-', '')
        headerSubKeys['comments'] = commentString
        headerSubKeys['tumor_subtype_upper'] = headerSubKeys['tumor_subtype'].upper()
        sHeader = sHeaderTemplate.safe_substitute(headerSubKeys)
        return sHeader

    def renderID(self, m):
        ID = '.'
        dbRS_rs = self._get_annotation_value(m, 'dbSNP_RS')
        if dbRS_rs != '':
            if len(dbRS_rs.split('|'))>1:
                dbRS_rs = sorted(dbRS_rs.split('|'))[0]
            if len(dbRS_rs.split(';'))>1:
                dbRS_rs = sorted(dbRS_rs.split(';'))[0]

            if dbRS_rs not in self.seenDbSNPs.keys():
                self.seenDbSNPs[dbRS_rs] = 0
            self.seenDbSNPs[dbRS_rs] += 1
            appendStr = ''
            if self.seenDbSNPs[dbRS_rs] > 1:
                appendStr = 't' + str(self.seenDbSNPs[dbRS_rs] - 1)
            ID = dbRS_rs + appendStr
        return ID

    def genotype(self, n_lod, t_lod):
        norm_lod_threshold = 0
        tumor_lod_threshold = -2.3
        # To determine the genotype for normal, bear in mind that the n_lod is related to reference.
        # TODO: Need to check normal_best_gt?
        gtN = '0/1'
        if n_lod > norm_lod_threshold:
            gtN = '0/0'
        gtT = '0/1'
        if t_lod < tumor_lod_threshold:
            gtT = '0/0'
        return gtN, gtT

    def determineSomaticStatus(self, gtN, gtT,qual):
        isSomatic = False
        ss = 'Germline'
        ssCode = '1'
        if (gtT == '0/0') and (gtN == '0/1'):
            isSomatic= True
            ss='LOH'
            ssCode= '3'
            qual = -qual
        elif (gtT == '0/1') and (gtN == '0/1'):
            isSomatic = False
            ss = 'Germline'
            ssCode= '1'
        elif (gtT == '0/1') and (gtN == '0/0'):
            isSomatic = True
            ss='Somatic'
            ssCode= '2'
        else:
            # Used to raise exception here.  Now we skip this case.  We should not be considering it anyway.
            #raise Exception('NO CALL: Normal=%s  (%s), Tumor=%s (%s) Line (%d)' % (gtN, n_lod, gtT, t_lod, ctr ))
            return None, None
        return ss, ssCode

    def _retrieve_alias_key(self, key):
        return self.fieldMap.get(key, key)

    def _get_annotation_value(self, m, key, default=None, is_blank_default=False):
        val = m.get(self._retrieve_alias_key(key), default)
        if is_blank_default and val.strip() == "":
            val = default
        return val

    def _generateFormatFieldWithValues(self,n_gt, n_alt, n_ref, mq0, bq='30', ss='2', ssc='.'):
        """ GT:AD:DP:FA:MQ0:BQ:SS:SSC
        0/1:17,1:18:0.056:0:30:2:255

        """
        ad = n_ref + ',' + n_alt
        dp = str(int(n_alt) + int(n_ref) )

        if int(dp) == 0:
            fa = '.'
        else:
            fa = '%.3f' % ( float(n_alt) / float(dp) )

        return '%s:%s:%s:%s:%s:%s:%s:%s' % (n_gt, ad, dp, fa, mq0, bq, ss, ssc)

    def _generateFilterField(self,m):
        judgement = "KEEP"
        if 'judgement' in m.keys():
            judgement = self._get_annotation_value(m, 'judgement')

        if judgement == 'KEEP':
            filter='PASS'
        else:
            filter='mf1'
        return filter

    def _generateRejectInfoField(self, isDbSnp, dp, mq0, isSomatic, vt, ss, vlsc ='.'):
        """
        Generate the actual info string for a given mutation call.  This method should only be called on mutation calls that were filtered out.

        Output:  String that follows the format for a passed mutation call.

        The info fields from a sample header...
    ##INFO=<ID=VT,Number=1,Type=String,Description="Somatic variant type">
    ##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of sample">
    ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic variant">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">
    ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
        """
        db = ''
        if isDbSnp:
            db='DB;'

        somatic = ''
        if isSomatic:
            somatic='SOMATIC;'

        #TODO: Need vlsc to be passed in properly.

        return '%sDP=%s;MQ0=%s;%sSS=%s;VT=%s;VLSC=%s' % (db, dp, mq0, somatic, ss, vt, vlsc)

    def _generatePassInfoField(self, isDbSnp, dp, mq0, isSomatic, vt, gene, vc, ss, transcriptID, vlsc='.'):
        """
        Generate the actual info string for a given mutation call.  This method should only be called on mutation calls that were NOT filtered out.

        Output:  String that follows the format for a passed mutation call.

         From the header...
    ##INFO=<ID=Gene,Number=1,Type=String,Description="Hugo Gene Symbol">
    ##INFO=<ID=VT,Number=1,Type=String,Description="Somatic variant type">
    ##INFO=<ID=VC,Number=1,Type=String,Description="Somatic variant classification">
    ##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of sample">
    ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic variant">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">
    ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
    ##INFO=<ID=TID,Number=1,Type=String,Description="Transcript ID">

    Example output (old -- no transcriptID):  DP=185;Gene=DDX20;MQ0=0;SOMATIC;SS=Somatic;VC=Missense;VT=SNP
        """
        db = ''
        if isDbSnp:
            db='DB;'

        if transcriptID == "":
            transcriptID = "."

        somatic = ''
        if isSomatic:
            somatic='SOMATIC;'

            #TODO: This needs to be corrected:  This is a dummy value that is not being calculated from any actual data.
            vlsc='255'

        return '%sDP=%s;Gene=%s;MQ0=%s;%sSS=%s;VC=%s;VT=%s;TID=%s;VLSC=%s' % (db, dp, gene, mq0, somatic, ss, vc, vt, transcriptID, vlsc)

    def _generateInfoField(self,m, filter, mq0, ss):
        dp = str( int(m['t_ref_count']) + int(m['t_alt_count']) + int(m.get('n_ref_count', '0')) + int(m.get('n_alt_count', '0')) )
        isSomatic = (ss == "LOH") or (ss == "Somatic")
        isDbSnp = (m['dbSNP_RS'] != '')
        if filter == 'PASS':
            gene = self._get_annotation_value(m, 'gene')
            if gene == "Unknown":
                gene = "."
            info = self._generatePassInfoField(isDbSnp, dp, mq0, isSomatic, self._get_annotation_value(m, 'variant_type'), gene, self._get_annotation_value(m, 'variant_classification'), ss, self._get_annotation_value(m, 'transcript_id'))
        else:
            info = self._generateRejectInfoField(isDbSnp, dp, mq0, isSomatic, self._get_annotation_value(m, 'variant_type'), ss)
        return info

    def _generateBQ(self,c, s):
        """ c = count, s = sum of base quality scores.  Both are strings"""
        if (int(c) == 0) or (int(s)==0):
            return "."
        return '%d' % (int(s)/int(c))

    def _generateFormatField(self):
        """
        :return: String that can be used for fields in the FORMAT field
        """
        return 'GT:AD:DP:FA:MQ0:BQ:SS:SSC'

    def _extract_lod(self, m, lod_annotation_name):
        """Determine lod, defaulting to 50 if unable to get necessary values.
        """
        rawLod = self._get_annotation_value(m, lod_annotation_name, '50')
        if rawLod == "." or rawLod == "":
            lod = 50
        else:
            lod = int(float(rawLod))
        return lod

    def _createMutRow(self,m):

        if m is None:
            return None

        # Calculate values
        # TODO: This is a sloppy way to do GT
        n_lod = self._extract_lod(m, 'init_n_lod')
        t_lod = self._extract_lod(m, 't_lod_fstar')


        # Nothing to call if both LODs are too low.
        if (t_lod < 0) and (n_lod < 0):
            return

        qual = max(t_lod, 0)

        gtN, gtT = self.genotype(n_lod, t_lod)

        ref,alt,new_start = MutUtils.retrievePrecedingBaseFromReference(m)

        ss,ssCode = self.determineSomaticStatus(gtN, gtT, qual)
        if ss is None or ssCode is None:
            return

        n_alt_count = self._get_annotation_value(m, 'n_alt_count', '0', is_blank_default=True)
        n_ref_count = self._get_annotation_value(m, 'n_ref_count', '0', is_blank_default=True)
        n_alt_sum = self._get_annotation_value(m, 'n_alt_sum', '0', is_blank_default=True)
        t_alt_count = self._get_annotation_value(m, 't_alt_count', '0', is_blank_default=True)
        t_ref_count = self._get_annotation_value(m, 't_ref_count', '0', is_blank_default=True)
        t_alt_sum = self._get_annotation_value(m, 't_alt_sum', '0', is_blank_default=True)
        mq0=0
        normalFormat = self._generateFormatFieldWithValues(gtN, n_alt_count, n_ref_count, mq0,
                                                           self._generateBQ(n_alt_count, n_alt_sum), ssCode)
        primaryFormat = self._generateFormatFieldWithValues(gtT, t_alt_count, t_ref_count, mq0,
                                                            self._generateBQ(t_alt_count, t_alt_sum), ssCode)

        filterVal = self._generateFilterField(m)
        info = self._generateInfoField(m,filterVal, mq0, ss)

        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	PRIMARY
        mutRow = dict()

        mutRow['CHROM'] = self._renderChrom(m.chr)
        mutRow['POS'] = new_start
        mutRow['ID'] = self.renderID(m)
        mutRow['REF'] = ref
        mutRow['ALT'] = alt
        mutRow['QUAL'] = qual
        mutRow['FILTER'] = filterVal
        mutRow['INFO'] = info
        mutRow['FORMAT'] = self._generateFormatField()
        mutRow['NORMAL'] = normalFormat
        mutRow['PRIMARY'] = primaryFormat

        return mutRow

    def _renderChrom(self, chr):
        """Render chromosome as used in TCGA VCF 1.1.

        NOTE: Only support hg19.  This method will need to be changed for other genome builds
        """
        chromList = map(str, range(0,23))
        chromList.extend(['X', 'Y', 'MT'])
        result = chr
        if chr not in chromList:
            result = '<' + chr + '>'
        return result

    def _handleMissingAnnotations(self, m):
        missingHeaderAnnotations = MutUtils.retrieveMissingAnnotations(m,
                                                                       TcgaVcfOutputRenderer.requiredHeaderAnnotations)
        missingMutAnnotations = MutUtils.retrieveMissingAnnotations(m, TcgaVcfOutputRenderer.requiredMutAnnotations)
        if len(missingHeaderAnnotations) > 0:
            sError = "The following annotations are required for rendering a TCGA VCF 1.1, but were not found: " + str(
                missingHeaderAnnotations)
            self.logger.error(sError)
            raise MissingAnnotationException(sError)
        if len(missingMutAnnotations) > 0:
            sError = "The following annotations important for rendering a TCGA VCF 1.1.  Proceeding... : " + str(
                missingMutAnnotations)
            self.logger.warn(sError)

    def _renderNoMutations(self, metadata, comments):
        """Render the case where we have no input mutations. """
        pass

    def _createVcfHeaderFilePtr(self, comments, m):
        self._handleMissingAnnotations(m)
        # fieldMapping = MutUtils.createFieldsMapping(headers=TcgaVcfOutputRenderer.requiredHeaderAnnotations, annotations=annotations,
        #                              alternativeDictionary=self.alternativeDictionary, isRenderInternalFields=False)
        # Initialize the output file and write a header.
        fp = file(self._filename, 'w')
        sHeader = self.createVcfHeader(m, "   ".join(comments))
        fp.write(sHeader)
        return fp

    def renderMutations(self, mutations, metadata=None, comments=None):
        if comments is None:
            comments = []

        outputHeaders = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'PRIMARY']

        # Create a list of annotation names and make sure to catch the case where there are no variants specified.
        try:
            m = mutations.next()
        except StopIteration as si:
            m = None

        if m is not None:
            fp = self._createVcfHeaderFilePtr(comments, m)
        else:
            fp = self._createVcfHeaderFilePtr(comments, metadata.asDict())

        if m is not None:
            fieldsUsed = self.alternativeDictionary.keys()

            annotations = MutUtils.getAllAttributeNames(m)
            self.fieldMap = MutUtils.createFieldsMapping(fieldsUsed, annotations, self.alternativeDictionary, True)

        # Write each row:
        ctr = 0
        unrenderableRows = 0
        tsvWriter = csv.DictWriter(fp, outputHeaders, delimiter="\t", lineterminator="\n")
        mutRow = self._createMutRow(m)

        if mutRow is not None:
            tsvWriter.writerow(mutRow)
            ctr += 1

        for m in mutations:
            if (ctr % 1000) == 0:
                self.logger.info("Processed " + str(ctr) + " mutations")
            mutRow = self._createMutRow(m)

            # We may not render all rows.
            if mutRow is not None:
                tsvWriter.writerow(mutRow)
            else:
                unrenderableRows += 1
            ctr += 1
        self.logger.info("Processed all " + str(ctr) + " mutations.  Could not render: " + str(unrenderableRows))

