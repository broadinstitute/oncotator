# LICENSE_GOES_HERE


import re

class AlignmentBlock(object):
    """Contains data for one block of a GAF alignment, representing one
    ungapped segment of one feature-composite pairwise alignment"""
    def __init__(self, featureStart, featureEnd, compositeStart, compositeEnd):
        self.featureStart = featureStart
        self.featureEnd = featureEnd
        self.compositeStart = compositeStart
        self.compositeEnd = compositeEnd

    def __str__(self):
        return(str.join("\t", (str(self.featureStart), str(self.featureEnd),
                               str(self.compositeStart),
                               str(self.compositeEnd))))
        
    def __repr__(self):
        return(self.__str__())


class Gaf(object):
    """Object wrapper for GAF data.  Each line of GAF data describes a pairwise
    alignment or mapping between a feature and a composite.  The composite is
    often the genome, but can also be a larger entity that the feature might be
    part of or might align to part of.

    the GAF fields are as follows:
        entryNumber: an index (1..N)
        featureId: name of the feature object
        featureType: type of the feature object
        featureDbSource: source of data on the feature
        featureDbVersion: version of the featureDbSource
        featureDbDate: release date of featureDbSource
        featureSeqFileName: name of a file containining feature sequence data
        compositeId: name of the composite object (often the genome)
        compositeDbSource: source of data on the composite
        compositeDbVersion: version of the compositeDbSource
        compositeDbDate: release date of compositeDbSource
        alignmentType: the type of alignment (often 'pairwise')
        featureCoordinates: a string of the feature coordinates in a
           comma-delimited list of blocks, in which each block is of the
           format <start>-<end>, and the coordinates are 1-based and increasing
           (i.e. if the composite is the genome and the feature is on the minue
           strand, start and end are listed w.r.t the positive strand with
           start < end).  For example: 1-5,7-12,18-200

        compositeCoordinates: a string of the same format as the feature
           coordinates, except that if the composite is the genome, then
           the string is of the format chrom:coords:strand, e.g. chr1:1-5:+

           Note that if the same feature and composite have more than one
           alignment, they have one Gaf object in which all of the
           featureCoordinates and compositeCoordinates are semicolon-delimited.
           
        gene: the gene or genes associated with this feature.  If there is
           more than one gene, they are delimited by semicolons.
        geneLocus: the coordinate locus string (chrom:chromStart-chromEnd:strand)
           for the gene associated with this feature.  If the feature is
           associated with multiple genes, their locus strings are concatenated
           and delimited with semicolons.
        featureAlias: contains other names for the feature, and is usually blank.
        featureInfo: contains auxiliary annotation data terms, depending on the
           feature type.
    """
    def __init__(self, line=None, entryNumber=0):
        self.entryNumber = entryNumber
        self.featureId = ""
        self.featureType = ""
        self.featureDbSource = ""
        self.featureDbVersion = ""
        self.featureDbDate = ""
        self.featureSeqFileName = ""
        self.compositeId = ""
        self.compositeType = ""
        self.compositeDbSource = ""
        self.compositeDbVersion = ""
        self.compositeDbDate = ""
        self.alignmentType = ""
        self.featureCoordinates = ""
        self.compositeCoordinates = ""
        self.gene = ""
        self.geneLocus = ""
        self.featureAliases = ""
        self.featureInfo = ""
        if line != None:
            tokens = line.rstrip().split("\t")
            if len(tokens) <= 1:
                print line
            assert(len(tokens) > 1)
            self._setFields(tokens)

    def _setFields(self, tokens):
        """Set the fields of this GAF object according to a tab-delimited
        line of GAF data"""
        self.entryNumber = tokens[0]
        self.featureId = tokens[1]
        self.featureType = tokens[2]
        self.featureDbSource = tokens[3]
        self.featureDbVersion = tokens[4]
        self.featureDbDate = tokens[5]
        self.featureSeqFileName = tokens[6]
        self.compositeId = tokens[7]
        self.compositeType = tokens[8]
        self.compositeDbSource = tokens[9]
        self.compositeDbVersion = tokens[10]
        self.compositeDbDate = tokens[11]
        self.alignmentType = tokens[12]
        self.featureCoordinates = tokens[13]
        self.compositeCoordinates = tokens[14]
        if len(tokens) > 15:
            self.gene = tokens[15]
        if len(tokens) > 16:
            self.geneLocus = tokens[16]
        if len(tokens) > 17:
            self.featureAliases = tokens[17]
        if len(tokens) > 18:
            self.featureInfo = tokens[18]

    def _getRow(self):
        """Generate a tuple containing all fields in string format"""
        row = [str(self.entryNumber),
               self.featureId, self.featureType, self.featureDbSource,
               self.featureDbVersion, self.featureDbDate,
               self.featureSeqFileName,
               self.compositeId, self.compositeType, self.compositeDbSource,
               self.compositeDbVersion, self.compositeDbDate,
               self.alignmentType, self.featureCoordinates,
               self.compositeCoordinates, self.gene, self.geneLocus,
               self.featureAliases, self.featureInfo]
        return(row)
               
    def __str__(self):
        """Generate a tab-delimited string representation of the data"""
        return(str.join("\t", self._getRow()))
        
    def __repr__(self):
        return(self.__str__())

    def write(self, fh):
        """Write this object to a file"""
        fh.write(str(self))
        fh.write('\n')

        
    def copy(self, gaf):
        """Copy the data from another GAF object"""
        self.entryNumber = gaf.entryNumber
        self.featureId = gaf.featureId
        self.featureType = gaf.featureType
        self.featureDbSource = gaf.featureDbSource
        self.featureDbVersion = gaf.featureDbVersion
        self.featureDbDate = gaf.featureDbDate
        self.featureSeqFileName = gaf.featureSeqFileName
        self.compositeId = gaf.compositeId
        self.compositeType = gaf.compositeType
        self.compositeDbSource = gaf.compositeDbSource
        self.compositeDbVersion = gaf.compositeDbVersion
        self.compositeDbDate = gaf.compositeDbDate
        self.alignmentType = gaf.alignmentType
        self.featureCoordinates = gaf.featureCoordinates
        self.compositeCoordinates = gaf.compositeCoordinates
        self.gene = gaf.gene
        self.geneLocus = gaf.geneLocus
        self.featureAliases = gaf.featureAliases
        self.featureInfo = gaf.featureInfo


    def gafToBlocks(self):
        """Generate a list of AlignmentBlock objects, in which each
        block represents an ungapped alignment segment between the
        feature and composite"""
        blockList = []
        (chrom, compositeCoords, strand) = self.compositeCoordinates.split(":")
        compositeSegments = compositeCoords.split(",")
        featureSegments = self.featureCoordinates.split(",")
        assert len(compositeSegments) == len(featureSegments)
        for ii in range(len(compositeSegments)):
            cBlocks = compositeSegments[ii].split("-")
            cStart = cBlocks[0]
            if len(cBlocks) > 1:
                cEnd = cBlocks[1]
            else:
                cEnd = cStart
            thisBlock = AlignmentBlock(0, 0, int(cStart), int(cEnd))
            blockList.append(thisBlock)
        for ii in range(len(blockList)):
            fBlocks = featureSegments[ii].split("-")
            fStart = fBlocks[0]
            if len(fBlocks) > 1:
                fEnd = fBlocks[1]
            else:
                fEnd = fStart
            if strand == '+':
                blockIdx = ii
            else:
                blockIdx = -1 - ii
            blockList[blockIdx].featureStart = int(fStart)
            blockList[blockIdx].featureEnd = int(fEnd)
            featureLength = blockList[blockIdx].featureEnd \
                - blockList[blockIdx].featureStart
            compositeLength = blockList[blockIdx].compositeEnd \
                - blockList[blockIdx].compositeStart
            assert featureLength == compositeLength
        return(blockList)

    def subsetCompositeCoordinates(self, start, end):
        """Given start and end coordinates that specify a segment of the
        feature coordinates (i.e. start >= feature start, end <= feature end),
        generate and return the composite string that corresponds to this
        segment
        """
        (chrom, blocks, strand) = self.compositeCoordinates.split(":")
        blockList = self.gafToBlocks()
        subsetBlockList = []
        for block in blockList:
            if block.featureStart >= start and block.featureEnd <= end:
                subsetBlockList.append(block)
            elif block.featureStart <= end and block.featureEnd >= start:
                if strand == '+':
                    if block.featureStart < start:
                        block.compositeStart = block.compositeStart \
                            + start - block.featureStart
                        block.featureStart = start
                    if block.featureEnd > end:
                        block.compositeEnd = block.compositeEnd \
                            - block.featureEnd + end
                        block.featureEnd = end
                else:
                    if block.featureStart < start:
                        block.compositeEnd = block.compositeEnd \
                            - start + block.featureStart
                        block.featureStart = start
                    if block.featureEnd > end:
                        block.compositeStart = block.compositeStart \
                            + block.featureEnd - end
                        block.featureEnd = end
                compositeLength = block.compositeEnd - block.compositeStart
                featureLength = block.featureEnd - block.featureStart
                assert compositeLength == featureLength
                subsetBlockList.append(block)
        compositeString = ""
        delimiter = ""
        for block in subsetBlockList:
            compositeString = "%s%s%d-%d" % (compositeString, delimiter, 
                                             block.compositeStart,
                                             block.compositeEnd)
            delimiter = ","
        compositeString = "%s:%s:%s" % (chrom, compositeString, strand)
        return(compositeString)

    def compositeCoordsToLocus(self):
        """Given a composite coordinate string, return a locus
        coordinates string by first parsing out the chromosome and
        strand, then splitting the coordinates themselves on dashes,
        and taking the first token (the overall start) and the last
        token (the overall end)
        """
        (chrom, coordinates, strand) = self.compositeCoordinates.split(":")
        coordinatePieces = coordinates.split("-")
        locus = "%s:%s-%s:%s" % (chrom, coordinatePieces[0], coordinatePieces[-1], strand)
        return(locus)

    def singleBaseCoordFixup(self, inputCoordinates):
        """Certain classes of GAF objects (SNPs, SNP probes) frequently
        represent a single base, such that their coordinate string is
        [<chrom>:]<start>-<start>[:<strand>].  Generate a representation of
        this string as [<chrom:>]<start>[:<strand>]
        """
        coordinateUnits = inputCoordinates.split(":")
        if len(coordinateUnits) == 1:
            coordinateBlocks = coordinateUnits[0].split(",")
        else:
            assert(len(coordinateUnits) == 3)
            coordinateBlocks = coordinateUnits[1].split(",")
        delimiter=""
        newCoordinates = ""
        for block in coordinateBlocks:
            (start, end) = block.split("-")
            if start == end:
                newCoordinates = "%s%s%s" % (newCoordinates,
                                             delimiter, start)
            else:
                newCoordinates = "%s%s%s" % (newCoordinates,
                                             delimiter, block)
            delimiter = ","
        if len(coordinateUnits) == 1:
            outputCoordinates = newCoordinates
        else:
            outputCoordinates = "%s:%s:%s" % (coordinateUnits[0],
                                              newCoordinates,
                                              coordinateUnits[2])
        return(outputCoordinates)
    


    def featureToComposite(self, feat, comp, inferStrand=True,
                           mergeAdjacent=False, debug=False):
        """Given two gafs with the same composite, return a gaf where the first
        is the feature and the second is the composite
        """
        if inferStrand:
            strand = feat.compositeCoordinates.split(":")[2]
        else:
            strand = '+'
        if debug:
            print "assigning", feat, "from", comp, " strand", strand
        # Do the feature and composite overlap at all?  If they don't (as
        # revealed by coordinate strings of length 0), return an empty class.
        self.alignFeatureToComposite(feat, comp, strand,
                                     mergeAdjacent=mergeAdjacent, debug=debug)
        if len(self.featureCoordinates) == 0:
            assert len(self.compositeCoordinates) == 0
            self.__init__()
        else:
            self.entryNumber = feat.entryNumber
            self.featureId = feat.featureId
            self.featureType = feat.featureType
            self.featureDbSource = feat.featureDbSource
            self.featureDbVersion = feat.featureDbVersion
            self.featureDbDate = feat.featureDbDate
            self.featureSeqFileName = feat.featureSeqFileName
            self.compositeId = comp.featureId
            self.compositeType = comp.featureType
            self.compositeDbSource = comp.featureDbSource
            self.compositeDbVersion = comp.featureDbVersion
            self.compositeDbDate = comp.featureDbDate
            self.alignmentType = feat.alignmentType
            self.gene = feat.gene
            self.geneLocus = feat.geneLocus
            self.featureAliases = feat.featureAliases
            self.featureInfo = feat.featureInfo

    def alignFeatureToComposite(self, feat, comp, strand,
                                mergeAdjacent=False, debug=False):
        """Given two GAFs with the same composite, generate the coordinate
        portion of a third GAF in which the first is the feature and the
        second is the composite
        """
        #
        # Step 1: go through all the feature blocks and all the composite 
        # blocks and make a list of any overlapping sub-blocks. 
        featureBlocks = feat.gafToBlocks()
        compositeBlocks = comp.gafToBlocks()
        featureToCompositeBlocks = list()
        for fb in featureBlocks:
            if debug:
                print("testing feature block from ", fb.featureStart,
                      "to", fb.featureEnd, "genomic", fb.compositeStart,
                      fb.compositeEnd)
            for cb in compositeBlocks:
                if fb.compositeStart <= cb.compositeEnd \
                       and cb.compositeStart <= fb.compositeEnd:
                    overlapStart = max(fb.compositeStart, cb.compositeStart)
                    overlapEnd = min(fb.compositeEnd, cb.compositeEnd)
                    if debug:
                        print("overlapping region from", overlapStart,
                              "to", overlapEnd)
                    if strand == '+':
                        thisFeatStart = fb.featureStart \
                                        + (overlapStart - fb.compositeStart)
                        thisFeatEnd = fb.featureEnd \
                                      + (overlapEnd - fb.compositeEnd)
                        if debug:
                            print("feature from ", fb.featureStart, "to",
                                  fb.featureEnd, "genomic", fb.compositeStart,
                                  "to", fb.compositeEnd, "coords", thisFeatStart,
                                  thisFeatEnd)
                        thisCompStart = cb.featureStart \
                                        + (overlapStart - cb.compositeStart)
                        thisCompEnd = cb.featureEnd \
                                      + (overlapEnd - cb.compositeEnd)
                        if debug:
                            print("plus strand: composite from ", cb.featureStart,
                                  "to", cb.featureEnd, "genomic",
                                  cb.compositeStart, "to", cb.compositeEnd,
                                  "coords", thisCompStart, thisCompEnd,
                                  "feature end", cb.featureEnd, "overlap end",
                                  overlapEnd, "composite end", cb.compositeEnd)
                    else:
                        thisFeatStart = fb.featureStart \
                                        + (fb.compositeEnd - overlapEnd)
                        thisFeatEnd = fb.featureEnd \
                                      + (overlapStart - fb.compositeStart)
                        thisCompStart = cb.featureStart \
                                        - (overlapEnd - cb.compositeEnd)
                        thisCompEnd = cb.featureEnd \
                                      - (overlapStart - cb.compositeStart)
                        if debug:
                            print("minus strand: composite from ",
                                  cb.featureStart, "to", cb.featureEnd, "genomic",
                                  cb.compositeStart, "to", cb.compositeEnd,
                                  "overlap", overlapStart, overlapEnd,"coords",
                                  thisCompStart, thisCompEnd, "feature end",
                                  cb.featureEnd, "composite end", cb.compositeEnd)
                            print("feature from ", fb.featureStart, "to",
                                  fb.featureEnd, "genomic", fb.compositeStart,
                                  "to", fb.compositeEnd, "overlap", overlapStart,
                                  overlapEnd,"coords (", thisFeatStart,
                                  thisFeatEnd,"), (", thisCompStart, thisCompEnd,
                                  ")", "feature end", fb.featureEnd,
                                  "composite end", fb.compositeEnd)
                    block = AlignmentBlock(thisFeatStart, thisFeatEnd,
                                           thisCompStart, thisCompEnd)
                    featureToCompositeBlocks.append(block)
        #
        # Step 2: sort the new block list by order of feature start coordinates
        featureToCompositeBlocks.sort(key=lambda AlignmentBlock:
                                      AlignmentBlock.featureStart)

        #
        # Step 3: Step through the list.  Anytime two blocks are directly adjacent
        # in both feature and composite space, merge them if mergeAdjacent is True
        if debug:
            print "working on list of ", len(featureToCompositeBlocks), "blocks"
        ii = 0
        while ii < len(featureToCompositeBlocks) - 1:
            thisBlock = featureToCompositeBlocks[ii]
            nextBlock = featureToCompositeBlocks[ii+1]
            removedNextBlock = False
            if debug:
                print "comparing blocks", ii, "and", ii+1,\
                      "endpoints", thisBlock.featureEnd, nextBlock.featureStart
            if nextBlock.featureStart == thisBlock.featureEnd + 1:
                if debug:
                    print("blocks adjacent in feature space.  Comparing",
                          "composite coords","(", thisBlock.compositeStart,
                          thisBlock.compositeEnd,") (", nextBlock.compositeStart,
                          nextBlock.compositeEnd,")" )
                if nextBlock.compositeStart == thisBlock.compositeEnd + 1:
                    if mergeAdjacent:
                        thisBlock.featureEnd = nextBlock.featureEnd
                        thisBlock.compositeEnd = nextBlock.compositeEnd
                        del featureToCompositeBlocks[ii+1]
                        if debug:
                            print "removing block", ii+1
                        removedNextBlock = True
                elif thisBlock.compositeStart == nextBlock.compositeEnd + 1:
                    if mergeAdjacent:
                        thisBlock.featureEnd = nextBlock.featureEnd
                        thisBlock.compositeStart = nextBlock.compositeStart
                        del featureToCompositeBlocks[ii+1]
                        if debug:
                            print "removing block", ii+1
                        removedNextBlock = True
            if not removedNextBlock:
                ii = ii + 1

        #
        # Step 4: format the list into GAF feature and composite
        # coordinate strings
        alignFeatCoords = ""
        alignCompCoords = ""
        delimiter = ""
        for block in featureToCompositeBlocks:
            #
            # If the feature is on the minus strand and the composite is the
            # genome, then build the strings in reverse orientation with
            # respect to feature and composite.  This reflects that the
            # first block of the feature will be the last block of the composite.
            #
            if strand == '-' and comp.featureType == 'genome':
                if block.featureStart != block.featureEnd:
                    alignFeatCoords = "%d-%d%s%s" % (block.featureStart,
                                                     block.featureEnd,
                                                     delimiter, alignFeatCoords)
                else:
                    alignFeatCoords = "%d%s%s" % (block.featureStart, delimiter,
                                                  alignFeatCoords)
                if block.compositeStart != block.compositeEnd:
                    alignCompCoords = "%d-%d%s%s" % (block.compositeStart,
                                                     block.compositeEnd,
                                                     delimiter, alignCompCoords)
                else:
                    alignCompCoords = "%d%s%s" % (block.compositeStart, delimiter,
                                                  alignCompCoords)
            else:
                if block.featureStart != block.featureEnd:
                    alignFeatCoords = "%s%s%d-%d" % (alignFeatCoords, delimiter,
                                                     block.featureStart,
                                                     block.featureEnd)
                else:
                    alignFeatCoords = "%s%s%d" % (alignFeatCoords, delimiter,
                                                  block.featureStart)
                if block.compositeStart != block.compositeEnd:
                    alignCompCoords = "%s%s%d-%d" % (alignCompCoords, delimiter,
                                                     block.compositeStart,
                                                     block.compositeEnd)
                else:
                    alignCompCoords = "%s%s%d" % (alignCompCoords, delimiter,
                                                  block.compositeStart)
            delimiter = ","
        self.featureCoordinates = alignFeatCoords
        self.compositeCoordinates = alignCompCoords

    def combine(self, newGaf):
        """Combine two GAF entries into one.  This is commonly done when
        one GAF data element turns out to have two or more genomic alignments.
        In such cases, separate GAF entries are built for each alignment,
        and then they are combined after building"""
        #
        # Append the feature and composite coordinate strings.
        self.featureCoordinates = "%s;%s" % (self.featureCoordinates,
                                             newGaf.featureCoordinates)
        self.compositeCoordinates = "%s;%s" % (self.compositeCoordinates,
                                             newGaf.compositeCoordinates)
        #
        # Append the gene, geneLocus, and featureInfo strings if the
        # contents of the new one are not already in the contents of self.
        # Note that this requires escaping some regular expression wildcard
        # characters that we might expect to find: ? in gene and + in geneLocus.
        if not re.search(re.sub("\?", "\?", newGaf.gene), self.gene):
            self.gene = "%s;%s" % (self.gene, newGaf.gene)
        if not re.search(re.sub("\+", "\+", newGaf.geneLocus), self.geneLocus):
            self.geneLocus = "%s;%s" % (self.geneLocus, newGaf.geneLocus)
        if not re.search(newGaf.featureInfo, self.featureInfo):
            self.featureInfo = "%s;%s" % (self.featureInfo, newGaf.featureInfo)

