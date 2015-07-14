from oncotator.MutationData import MutationData

class MutationDataFactory(object):
    """
    Stateful class for producing mutation data.  This can be passed into Mutation Creators as a way of defining an initial state
    for the mutations that the InputCreator will eventually render.
    """

    def __init__(self, allow_overwriting=False):
        self._allow_overwriting = allow_overwriting

    def create(self, chr="", start="", end="", ref_allele="", alt_allele="", build=""):
        return MutationData(chr, start, end, ref_allele, alt_allele, build, new_required=not self._allow_overwriting)

    @staticmethod
    def default_create(chr="", start="", end="", ref_allele="", alt_allele="", build="", allow_overwriting=False):
        """Convenience method mostly for unit tests.  Generally, use an instance of the MutationDataFactory and call
         mut_factory.create()
        """
        return MutationData(chr, start, end, ref_allele, alt_allele, build, new_required=not allow_overwriting)

