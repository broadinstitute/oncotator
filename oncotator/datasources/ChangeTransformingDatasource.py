from oncotator.datasources.Datasource import Datasource

class ChangeTransformingDatasource(Datasource):
    def __init__(self, src_file, title='', version=None, use_binary=True):
        self.title = title
        self.version = version
        self.src_file = src_file

    def annotate_mutation(self, mutation):
        raise NotImplementedError("Not implemented ever")