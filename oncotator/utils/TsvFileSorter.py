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


import collections
from collections import OrderedDict
import heapq
import os
import tempfile
import itertools

__author__ = 'lichtens'

class TsvFileSorter(object):
    """ Static class for sorting a tsv file in place.

    This code is an adaptation of the cookbook code found at:
    http://code.activestate.com/recipes/576755-sorting-big-files-the-python-26-way/
    """



    def __init__(self, fieldNames = ["chrom","pos","sampleName"], delimiter = '\t'):
        self.fieldNames = fieldNames
        self.fieldNames.append("line")
        self._Pair = collections.namedtuple(typename="Pair", field_names=["key","value"])
        self.delimiter = delimiter
        self.headerList = None
        self.columnPos = dict()

    def __createLineDict(self, hdrList, line):
        toks = line.split(self.delimiter)
        result = OrderedDict()
        ctr = 0
        for hdr in hdrList:
            result[hdr] = toks[ctr]
            ctr += 1
        result['line'] = line
        return result

    def __merge(self, *partitions):
        """

        :param partitions:
        """
        # iterables = [(self._Pair(key=, value=line) for line in partition)
        #              for partition in partitions]
        iterables = []
        for partition in partitions:
            partitionList = []
            for line in partition:
                d = self.__createLineDict(hdrList=self.headerList, line=line)
                # TODO: key_list is not quite correct.  Modify to specify a callable (dict to tuple)
                key_list = [d[self.fieldNames[0]].lower(), int(d[self.fieldNames[1]]), int(d[self.fieldNames[2]])]
                partitionList.append(self._Pair(key=tuple(key_list), value=line))
            iterables.append(partitionList)

        # Merge multiple sorted inputs
        for pair in heapq.merge(*iterables):
            print pair.value
            yield pair.value


    def sortFile(self,readfilename, outputFilename, iswriteHeader=True, length=1280000):
        tempdirs = list([tempfile.gettempdir()])
        partitions = list()
        isHashHeader = False
        try:
            header = None
            partitionCtr = 0
            with open(name=readfilename, mode='rb', buffering=64 * 1024) as fp:
                iterable = iter(fp)
                for tempdir in itertools.cycle(tempdirs):
                    lines = list(itertools.islice(iterable, length)) # returns an iterator of size length
                    if header is None:
                        header = lines[0].rstrip() # first line is the header
                        if header.find("#") == 0:
                            header = header.replace("#", "")
                            isHashHeader = True
                        self.headerList = header.split(self.delimiter)
                        print(self.fieldNames[0])
                        for f in self.fieldNames:
                            if f == "line":
                                continue
                            self.columnPos[f] = self.headerList.index(f)
                        lines = lines[1:len(lines)]
                    for index in range(len(lines)):
                        lines[index] = self.__createLineDict(hdrList=self.headerList, line=lines[index].rstrip('\n'))

                    if not lines:
                        break
                    # lines = sorted(lines, key=operator.attrgetter('chrom', 'pos', 'sampleName')) # sort the list of lines
                    # TODO: Determine type of the input columns and then convert to lowercase str or use int/float
                    lines = sorted(lines, key=lambda t: (t[self.fieldNames[0]].lower(), int(t[self.fieldNames[1]]), int(t[self.fieldNames[2]]))) # sort the list of lines
                    partition = open(name=os.path.join(tempdir, '%06i' % len(partitions)), mode='w+b',
                                     buffering=64 * 1024)
                    partitions.append(partition)
                    partition.write('\n'.join([record['line'] for record in lines]) + '\n')
                    partition.flush()
                    partition.seek(0)
                    partitionCtr += 1
                    print(str(partitionCtr))
            with open(name=outputFilename, mode='wb', buffering=64 * 1024) as fp:
                if iswriteHeader:
                    if isHashHeader:
                        fp.write("#")
                    fp.write(header + "\n")
                print("Merging...")
                fp.writelines(self.__merge(*partitions)) # generators are allowed as inputs to writelines funtion
        finally:
            for partition in partitions:
                try:
                    partition.close()
                    os.remove(path=partition.name)
                except Exception:
                    pass
