class DataConfig:
    def __init__(self):
        self._bootstrapThreshold = 0
        self._lsThreshold = 60
        self._windowSize = 200
        self._stepSize = 100
        self._dataNames = ['ALLSKY_SFC_SW_DWN_newick', 'T2M_newick', 'QV2M_newick', 'PRECTOTCORR_newick',
                           'WS10M_newick']
        self._referenceGeneFile = '../datasets/small_seq.fasta'
        self._fileName = '../datasets/geo.csv'
        self._specimen = 'id'
        self._names = ['id', 'ALLSKY_SFC_SW_DWN', 'T2M', 'QV2M', 'PRECTOTCORR', 'WS10M']
        self._bootstrapList = []
        self._data = []
        self._bootstrapAmount = 100
        self._referenceGeneDir = '../datasets/'
        self._makeDebugFiles = True

    # getters
    def get_bootstrapThreshold(self):
        return self._bootstrapThreshold

    def get_lsThreshold(self):
        return self._lsThreshold

    def get_windowSize(self):
        return self._windowSize

    def get_stepSize(self):
        return self._stepSize

    def get_dataNames(self):
        return self._dataNames

    def get_referenceGeneFile(self):
        return self._referenceGeneFile

    def get_fileName(self):
        return self._fileName

    def get_specimen(self):
        return self._specimen

    def get_names(self):
        return self._names

    def get_bootstrapList(self):
        return self._bootstrapList

    def get_data(self):
        return self._data

    def get_bootstrapAmount(self):
        return self._bootstrapAmount

    def get_referenceGeneDir(self):
        return self._referenceGeneDir

    def get_makeDebugFiles(self):
        return self._makeDebugFiles

    # setters
    def set_bootstrapThreshold(self, value):
        self._bootstrapThreshold = value

    def set_lsThreshold(self, value):
        self._lsThreshold = value

    def set_windowSize(self, value):
        self._windowSize = value

    def set_stepSize(self, value):
        self._stepSize = value

    def set_dataNames(self, value):
        self._dataNames = value

    def set_referenceGeneFile(self, value):
        self._referenceGeneFile = value

    def set_fileName(self, value):
        self._fileName = value

    def set_specimen(self, value):
        self._specimen = value

    def set_names(self, value):
        self._names = value

    def set_bootstrapList(self, value):
        self._bootstrapList = value

    def set_data(self, value):
        self._data = value

    def set_bootstrapAmount(self, value):
        self._bootstrapAmount = value

    def set_referenceGeneDir(self, value):
        self._referenceGeneDir = value

    def set_makeDebugFiles(self, value):
        self._makeDebugFiles = value
