import collections
import logging
import more_itertools
from oncotator.MutationData import MutationData
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.utils.SampleNameSelector import SampleNameSelector
from oncotator.utils.iterutils import flatmap

__author__ = 'louisb'


class OnpQueue(object):
    """
    Bookkeeping class to maintain the mutations waiting to be combined
    """

    def __init__(self, mutations):
        """
        Initialize an new queue with a MutationData iterator
        :param mutations: any MutationData producing Iterator
        """
        self.mutations = more_itertools.peekable(mutations)
        self.sns = SampleNameSelector(self.mutations.peek())
        self.queue = collections.defaultdict(list)
        self.indel_queue = []
        self.last = 0
        self.logger = logging.getLogger(__name__)
        self.warned_about_order = False

    @staticmethod
    def _create_start_position_dict(mutations):
        """
        Create a start_position -> mutation dict
        :param mutations: a collection of MutationData
        :return: a dictionary containing all the input MutationData grouped by start postion
        """
        assert (mutations is not None)
        starts = collections.defaultdict(list)
        for mut in mutations:
            starts[int(mut.start)] += [mut]
        return starts

    @staticmethod
    def _paths(finished_paths, path_so_far, start, muts):
        """Return all paths from the start position through the mutation graph
        :param finished_paths: completed paths
        :param path_so_far: the accumulated mutation->mutation path so far
        :param start: the start position to travers the muts from
        :param muts: a dictionary in the form {start_position: [Mutation]}
        :return: All paths through adjacent mutations starting with mutations at chromosome position start
        """
        if muts == [] or start not in muts:
            finished_paths.append(path_so_far)
        else:
            # return reduce(operator.concat, lambda mut: OnpCombiner._paths(path + [mut], mut.end+1, muts), [])
            # path =  map(lambda mut: OnpQueue._paths(path + [mut], int(mut.end)+1, muts), muts[start])
            for mut in muts[start]:
                OnpQueue._paths(finished_paths, path_so_far + [mut], int(mut.end) + 1, muts)
            return finished_paths
            #return reduce(operator.concat, path)

    def _add(self, mutation):
        variant_type = TranscriptProviderUtils.infer_variant_type(mutation.ref_allele, mutation.alt_allele)
        # only combine ONPs, not indels
        if not TranscriptProviderUtils.is_xnp(variant_type):
            self.indel_queue.append(mutation)
        else:
            self.queue[self.sns.getSampleName(mutation)].append(mutation)

    def _walk_mutation_paths(self, muts):
        """
        Find all paths through adjacent mutations and return those as combined mutations

        Find the first mutations by chromosome position and compute all paths through adjacent mutations reachable from them.
        If there are any nodes that were not reached, choose the first position with unreached nodes and repeat
        :param muts: a list of mutations to walk throught
        :return: a list of new mutations combined from
        """
        unreached = muts
        paths = []
        starts = self._create_start_position_dict(muts)
        while unreached:
            paths += self._paths([], [], min([int(mut.start) for mut in unreached]), starts)
            reached = [mut for path in paths for mut in path]
            unreached = [mut for mut in muts if mut not in reached]

        paths = [OnpQueue._combine_mutations(path) for path in paths]
        return paths

    def _dump_all(self):
        results = []
        for (sample, muts) in self.queue.iteritems():
            results += self._walk_mutation_paths(muts)

        #add all stored up indels
        results += self.indel_queue or []
        self.indel_queue = []
        results.sort(key=lambda x: (int(x.start), int(x.end)))

        self.queue.clear()
        return results


    def _get_all_values(self):
        return [j for i in self.queue.values() for j in i]

    def _is_adjacent_to_any_xnp(self, new_mutation):
        return self._is_adjacent(new_mutation, self._get_all_values())

    def _is_adjacent(self, new_mutation, mutations):
        if mutations:
            ends = [int(x.end) for x in mutations]
            return int(new_mutation.start) <= 1 + max(ends)
        else:
            return False

    @staticmethod
    def _combine_mutations(mutations):
        """
        Merge multiple adjacent mutations into a single new mutation.

        :param mutations: an ordered list of MutationData
        :returns a new MutationData

        :warning: _combine_mutations does not make any attempt to sanity check input mutations
        it will happily combine overlapping and non-adjacent mutations on disparate chromosomes
        """
        if len(mutations) == 0:
            return None
        if len(mutations) == 1:
            return mutations[0]

        # special logic for the attributes
        start = min([mut.start for mut in mutations])
        end = max([mut.end for mut in mutations])
        chr = mutations[0].chr
        ref = "".join([mut.ref_allele for mut in mutations])
        alt = "".join([mut.alt_allele for mut in mutations])
        build = "|".join(set([x.build for x in mutations]))

        #create the new mutation
        newmut = MutationData(chr=chr, start=start, end=end, ref_allele=ref, alt_allele=alt, build=build)

        #add annotations to the mutation
        allAnnotations = set(flatmap(lambda x: x.keys(), mutations))
        annotationNames = allAnnotations - set(mutations[0].getAttributeNames())
        for annotName in annotationNames:
            annotations = []
            for mut in mutations:
                try:
                    annotations.append(mut.getAnnotation(annotName))
                except KeyError:
                    pass

            values = sorted((set([x.getValue() for x in annotations if x.getValue()])))
            value = "|".join(values)
            tags = sorted(set(flatmap(lambda x: x.getTags(), annotations)))
            source = annotations[0].getDatasource()
            datatype = annotations[0].getDataType()
            number = annotations[0].getNumber()
            description = annotations[0].getDescription()
            newmut.createAnnotation(annotationName=annotName, annotationValue=value, annotationSource=source,
                                    annotationDataType=datatype, annotationDescription=description,
                                    tags=tags, number=number)
        return newmut

    def _combine_with_indels(self, output):
        """add the indels to output and sort by start position"""
        output += self.indel_queue or []
        self.indel_queue = []
        output.sort(key=lambda x: int(x.start))
        return

    def get_combined_mutations(self):
        """
        :return: a generator yielding mutations, adjacent SNPs,DNPs, and ONPs will be merged together.
        """
        # assumes mutations are sorted by start position and then sample
        #if they're not, it won't find DNPs
        last_chr = -1
        last_start = -1
        for mut in self.mutations:
            output = []
            #if we're on a new chromosome, dump all mutations, then add the new one to the queue
            if mut.chr != last_chr:
                output = self._dump_all()
                self._add(mut)
            #if we're at the same start position, add the new mutation to the queue
            elif mut.start == last_start:
                self._add(mut)
            #if we are at a new position on the same chromosome
            elif  self._is_adjacent_to_any_xnp(mut):
                #  if we are adjacent/overlapping to one of our existing positions
                #   add the mutation
                if not self.warned_about_order and int(mut.start) < last_start:
                    self.logger.warn("Mutations are not sorted by start position, this may cause unexpected behavior or "
                                     "increased memory requirements.  It is recommended that your sort any files that you"
                                     "using with --infer-onps by position and sample name.")
                    self.warned_about_order = True
                self._add(mut)
            #  if we are not adjacent to any existing queue position,
            #   dump mutations, then add the mutation
            else:
                output = self._dump_all()
                self._add(mut)
            last_chr = mut.chr
            last_start = mut.start

            for mut in output:
                yield mut

        #when we're finished, be sure to dump any last mutations
        output = self._dump_all()
        for mut in output:
            yield mut
