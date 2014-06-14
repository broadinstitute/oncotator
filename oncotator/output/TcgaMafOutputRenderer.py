# LICENSE_GOES_HERE
from oncotator.utils.OptionConstants import OptionConstants


"""
Created on Nov 7, 2012

@author: lichtens
"""

from OutputRenderer import OutputRenderer
import logging
import csv
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.version import VERSION
from oncotator.utils.ConfigUtils import ConfigUtils
from collections import OrderedDict



class TcgaMafOutputRenderer(OutputRenderer):
    """
    
    Render a generator or list of mutations into a TCGA MAF file.  
    
    TCGA MAF specification can be found at: https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
    
    Version specified in the config file in the "general" section.
    """
    def getTcgaMafVersion(self):
        return self.config.get("general", "version")

    def __init__(self, filename, configFile="tcgaMAF2.4_output.config", other_options=None):
        """
        TODO: Need functionality for not prepending the i_ on internal fields.
        """
        options = dict() if other_options is None else other_options

        self._filename = filename
        self.logger = logging.getLogger(__name__)
        self.config = ConfigUtils.createConfigParser(configFile)

        self.logger.info("Building alternative keys dictionary...")
        self.alternativeDictionary = ConfigUtils.buildAlternateKeyDictionaryFromConfig(self.config)
        
        #TODO: Read missing options from the config file or specify that error should be thrown.
        self.options = options

        self._prepend = self.config.get("general", "prepend")
        if self.options.get(OptionConstants.NO_PREPEND, False):
            self._prepend = ""

        self.exposedColumns = set(self.config.get("general", "exposedColumns").split(','))
    
    def lookupNCBI_Build(self, build):
        """ If a build number exists in the config file, use that.  Otherwise, use the name specified. """
        if not self.config.has_option("genomeBuild", build):
            return build
        self.config.get("genomeBuild", build, vars={"genomeBuild":build})
    
    def _createMutationRow(self, m, headers, fieldMapping):
        """ Create a single mutation dictionary (i.e. render a line).  A dictionary as per the csv library.
        Headers will usually be the fieldMapping keys, but extra parameter is provided here in case subset is desired.
        Also, allows caching of the keys ahead of time. """
        row = dict()
        for h in headers:
            annotation = fieldMapping[h]
            value = m.get(annotation, "__UNKNOWN__")
            row[h] = value
        return row

    def _determine_new_allele_if_blank(self, d, allele_key, new_value):
        """

        :param d: dictionary of column names
        :param allele_key: key to replace if "" or does not exist.
        :param new_value: value to use if "" or does not exist
        :return:
        """
        result = d.get(allele_key, new_value)
        if result.strip() == "":
            result = new_value
        return result

    def _update_validation_values(self, row):
        """ If Validation_Status  == "Valid" then
          Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2 cannot  be null
         If Mutation_Status == "Somatic" and Validation_Status == "Valid", then
          Match_Norm_Validation_Allele1 == Match_Norm_Validation_Allele2 == Reference_Allele and (Tumor_Validation_Allele1 or Tumor_Validation_Allele2) != Reference_Allele

         If Validation_Status == "Invalid" then
          Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2 cannot be null AND Tumor_Validation_Allelle1 == Match_Norm_Validation_Allele1 AND Tumor_Validation_Allelle2 == Match_Norm_Validation_Allele2  (Added as a replacement for 8a as a result of breakdown)


        IMPORTANT: The input parameter is altered.

        :param row: dict with name value pairs that include the TCGA MAF columns.  This is usually not the mutation.
        """
        if row['Validation_Status'] == "Valid":
            if row['Mutation_Status'] == "Somatic":
                row['Tumor_Validation_Allele1'] = self._determine_new_allele_if_blank(row, 'Tumor_Validation_Allele1',
                                                                                      row['Reference_Allele'])
                row['Tumor_Validation_Allele2'] = self._determine_new_allele_if_blank(row, 'Tumor_Validation_Allele2',
                                                                                      row['Tumor_Seq_Allele2'])
                row['Match_Norm_Validation_Allele1'] = self._determine_new_allele_if_blank(row,
                                                                                           'Match_Norm_Validation_Allele1',
                                                                                           row['Reference_Allele'])
                row['Match_Norm_Validation_Allele2'] = self._determine_new_allele_if_blank(row,
                                                                                           'Match_Norm_Validation_Allele2',
                                                                                           row['Reference_Allele'])

        if row['Validation_Status'] == "Invalid":

            # Only valid mutation status value is None for an invalid mutation
            if row['Mutation_Status'] != "None":
                row['Mutation_Status'] = "None"

            # If the alleles are blank, populate properly for invalid mutation.  Basically, everything becomes reference
            row['Match_Norm_Validation_Allele1'] = self._determine_new_allele_if_blank(row,
                                                                                       'Match_Norm_Validation_Allele1',
                                                                                       row['Reference_Allele'])
            row['Match_Norm_Validation_Allele2'] = self._determine_new_allele_if_blank(row,
                                                                                       'Match_Norm_Validation_Allele2',
                                                                                       row['Reference_Allele'])
            row['Tumor_Validation_Allele1'] = self._determine_new_allele_if_blank(row, 'Tumor_Validation_Allele1',
                                                                                  row['Match_Norm_Validation_Allele1'])
            row['Tumor_Validation_Allele2'] = self._determine_new_allele_if_blank(row, 'Tumor_Validation_Allele2',
                                                                                  row['Match_Norm_Validation_Allele2'])


    def _writeMutationRow(self, dw, fieldMap, fieldMapKeys, m):
        """ If this row should be rendered, then write it to the given DictWriter

        Additionally, apply corrections needed to make this a valid TCGA MAF.

        This method must be called as a last step before writing the output, as it relies on the output row, as opposed
            to the annotated mutation.

        :param dw: DictWriter
        :param fieldMap:
        :param fieldMapKeys:
        :param m:
        :return:
        """
        row = self._createMutationRow(m, fieldMapKeys, fieldMap)
        if row['Entrez_Gene_Id'] == "":
            row['Entrez_Gene_Id'] = "0"

        if row['Entrez_Gene_Id'] == "0" and row['Hugo_Symbol'] != "Unknown":
            logging.getLogger(__name__).warn("Entrez Gene ID was zero, but Hugo Symbol was not Unknown.  Is the HGNC and/or Transcript datasource complete?")

        self._update_validation_values(row)

        dw.writerow(row)

    def renderMutations(self, mutations, metadata=None, comments=None):
        """ Returns a file name pointing to the maf file that is generated. """
        if metadata is None:
            metadata = OrderedDict()

        if comments is None:
            comments = []

        self.logger.info("TCGA MAF output file: " + self._filename)
        self.logger.info("Render starting...")

        requiredColumns = self.config.get("general", "requiredColumns").split(',')
        optionalColumns = self.config.get("general", "optionalColumns").split(',')

        # Create the header list, making sure to preserve order.
        headers = requiredColumns
        headers.extend(optionalColumns)

        # Create a list of annotation names
        try:
            m = mutations.next()
            annotations = MutUtils.getAllAttributeNames(m)
        except StopIteration as si:

            # There are no mutations, so use the config file and metadata to determine what columns to output
            metadataAnnotations = metadata.keys()
            annotations = set(headers).union(metadataAnnotations)
            m = None

        # Create a mapping between column name and annotation name
        fieldMap = MutUtils.createFieldsMapping(headers, annotations, self.alternativeDictionary,
                                                self.config.getboolean("general", "displayAnnotations"),
                                                exposedFields=self.exposedColumns, prepend=self._prepend)
        fieldMapKeys = fieldMap.keys()
        internalFields = sorted(list(set(fieldMapKeys).difference(headers)))
        headers.extend(internalFields)
        
        # Initialize the output file and write a header.
        fp = file(self._filename, 'w')
        fp.write("#version " + self.getTcgaMafVersion() + "\n")
        
        for c in comments:
            fp.write("## " + c + "\n")
        
        # Initialize a csv DictWriter
        # Remove headers that start with "_"
        dw = csv.DictWriter(fp, headers, delimiter="\t", lineterminator="\n")
        dw.writeheader()
        ctr = 0

        try:
            # Add the NCBI build
            if m is not None:
                m.createAnnotation('ncbi_build', self.lookupNCBI_Build(m.build), annotationSource="OUTPUT")
                self._writeMutationRow(dw, fieldMap, fieldMapKeys, m)
                ctr += 1

            for m in mutations:

                # Add the NCBI build
                m.createAnnotation('ncbi_build', self.lookupNCBI_Build(m.build), annotationSource="OUTPUT")
                self._writeMutationRow(dw, fieldMap, fieldMapKeys, m)
                
                # Update mutation count and log every 1000 mutations
                ctr += 1
                if (ctr % 1000) == 0:
                    self.logger.info("Rendered " + str(ctr) + " mutations.")
        except Exception as e:
            import traceback
            self.logger.error(traceback.format_exc())
            self.logger.error("Error at mutation " + str(ctr) + " " + str([m.chr,m.start,m.end]) + ": ")
        
        fp.close()
        self.logger.info("Rendered all " + str(ctr) + " mutations.")
        return self._filename