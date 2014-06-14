# LICENSE_GOES_HERE
from VcfOutputConfigTableCreator import VcfOutputConfigTableCreator
from VcfInputConfigTableCreator import VcfInputConfigTableCreator


class ConfigTableCreatorFactory():

    @staticmethod
    def getConfigTableCreatorInstance(configTableType):
        if configTableType in ("input_vcf",):
            return VcfInputConfigTableCreator()
        elif configTableType in ("output_vcf",):
            return VcfOutputConfigTableCreator()