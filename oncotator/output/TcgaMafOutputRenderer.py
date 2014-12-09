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

        self._is_entrez_id_message_logged = False
    
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

        # Use HGNC Entrez Gene ID, if available and nothing else has populated it.,
        if row['Entrez_Gene_Id'] == "" and m.get('HGNC_Entrez Gene ID(supplied by NCBI)', "") != "":
            row['Entrez_Gene_Id'] = m.get('HGNC_Entrez Gene ID(supplied by NCBI)')

        if row['Entrez_Gene_Id'] == "":
            row['Entrez_Gene_Id'] = "0"

        if not self._is_entrez_id_message_logged and row['Entrez_Gene_Id'] == "0" and row['Hugo_Symbol'] != "Unknown":
            logging.getLogger(__name__).warn("Entrez Gene ID was zero, but Hugo Symbol was not Unknown.  Is the HGNC and/or Transcript datasource complete?")
            self._is_entrez_id_message_logged = True
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
            self.logger.error("Error at mutation " + str(ctr) + " " + str([m.chr,m.start,m.end,m.ref_allele,m.alt_allele]) + ": ")
            self.logger.error("Incomplete: rendered %d mutations." % (ctr))
            fp.close()
            raise e
        
        fp.close()
        if self._is_entrez_id_message_logged:
            logging.getLogger(__name__).warn("Some Entrez_Gene_IDs may be missing for valid Hugo Symbols in this TCGA MAF.")
        self.logger.info("Rendered all " + str(ctr) + " mutations.")
        return self._filename