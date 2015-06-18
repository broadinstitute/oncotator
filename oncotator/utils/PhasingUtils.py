
class PhasingUtils(object):
    """
    Collection of logic for determining the phasing information in mutation(s).
    """

    # These must be defined in the VCF input config file and the VCF output config file.  Typically, with PGT and PID
    PHASING_GT = "phasing_genotype"
    PHASING_ID = "phasing_id"

    @staticmethod
    def is_in_phase(mut1, mut2, unknown_val=True):
        """
        Determine if the two mutations should be considered in phase with each other.
        Currently, this class only supports the M2 convention using phasing_genotype and phasing_id annotations.

        "phasing information" below is defined as having both the phasing_genotype and phasing_id annotations.  If one
         of these annotations is missing and the other exists, then Oncotator behavior is undefined.

        Rules are as follows:
            0.1)  Missing phasing information cannot be in phase with a mutation containing populated phasing
                information.  Therefore, false will always be returned
            0.2)  If both mutations are missing phasing information, then unknown_val will be returned
            1) Check the phasing ID to determine whether there is a match.


        :param mut1: The first mutation (appears before mut2 in genomic space)
        :param mut2: The second mutation (appears after mut1 in genomic space)
        :type mut1 MutationData
        :type mut2 MutationData
        :param unknown_val: if there is no phasing information in the two mutations, return this boolean value.
        :type unknown_val bool
        :return true if the two mutations are in the same phase.  If there is not enough information, return the
            specified unknown_val

        :rtype bool
        """
        does_mut1_have_phasing = PhasingUtils.has_phasing_information(mut1)
        does_mut2_have_phasing = PhasingUtils.has_phasing_information(mut2)

        # Rule 0.1
        if does_mut1_have_phasing ^ does_mut2_have_phasing:
            return False

        # Rule 0.2
        if not does_mut1_have_phasing and not does_mut2_have_phasing:
            return unknown_val

        if mut1[PhasingUtils.PHASING_ID] == mut2[PhasingUtils.PHASING_ID]:
            return True

        return False


    @staticmethod
    def has_phasing_information(mut):
        """
        Check that there is phasing information for this mutation.

        This currently just checks that there is both a phasing_genotype and phasing_id annotation.

        :param mut:
        :type mut MutationsData
        :return:
        """
        return PhasingUtils.PHASING_GT in mut.keys() and PhasingUtils.PHASING_ID in mut.keys()