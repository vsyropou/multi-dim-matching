

class MultiDimMatching:
    def __init__(self, **kwargs):
        SourceDataSet = kwargs.pop('SourceDataSet', None)
        TargetDataSet = kwargs.pop('TargetDataSet', None)
        Weighted      = kwargs.pop('Weighted', {'Source':None, 'Target':None})
        Variables     = kwargs.pop('Variables', None)
        Verbose       = kwargs.pop('Verbose', False)

        nBins   = kwargs.pop( 'NumBins', 1000 )
        mean    = kwargs.pop( 'Mean', 0. )
        sigma   = kwargs.pop( 'Sigma', 1. )
        outFile = kwargs.pop( 'OutFileName', 'out.root' )
        rootLib = kwargs.pop( 'UseRootLibs', False )

        # check libraries availability
        if rootLib:
            try:
                from ROOT import TMath
            except ImportError:('Cannot locate ROOT module library in $PYTHONPATH')
        else:
            try:
                from scipy.special import erfinv
            except ImportError:('Cannot find locate scipy module in $PYTHONPATH')

        # input checks
        print 'P2VV-INFO: Initialising Distributions matching for variables %s'%Variables.values()
        from ROOT import TTree
        assert type(SourceDataSet) == type(TargetDataSet) == TTree, 'P2VV-INFO: Input datasets are not of type %s'%TTree
        assert Variables, 'P2VV-INFO: Specify variables.'

        # read trees into python lists
        from P2VV import RooFitDecorators
        from P2VV.Load import P2VVLibrary

        # heling stuff
        def append_row(data,row):
            for key, var in Variables.iteritems():
                data[key] += [ row[key] ]

        if Weighted['Source']:
            _readSourceRow = lambda ev,var: (getattr(ev,var),getattr(ev,Weighted['Source']))
        else:
            _readSourceRow = lambda ev,var: getattr(ev,var)

        if Weighted['Target']:
            _readTargetRow = lambda ev,var: (getattr(ev,var),getattr(ev,Weighted['Target']))
        else:
            _readTargetRow = lambda ev,var: getattr(ev,var)

        sourceData = { var : [] for var in Variables.keys() }
        targetData = { var : [] for var in Variables.keys() }

        print 'P2VV-INFO: Reading source branches'
        for ev in SourceDataSet:
            row = { key : _readSourceRow(ev,var) for key, var in Variables.iteritems() }
            append_row(sourceData,row)

        print 'P2VV-INFO: Reading target branches'
        for ev in TargetDataSet:
            row = { key : _readTargetRow(ev,var) for key, var in Variables.iteritems() }
            append_row(targetData,row)

        self._sourceData = sourceData
        self._targetData = targetData
        self._weighted   = Weighted

        self._SourceDataSet = SourceDataSet
        self._TargetDataSet = TargetDataSet

        self._variables = Variables.keys()
        self._dimention = len(self._variables)

        self._verbose = Verbose

        self._nBins   = nBins
        self._mean    = mean
        self._sigma   = sigma
        self._outFile = outFile
        self._rootLib = rootLib

#    @classmethod
    # def fromOtherContainerType(self, SourceDataSet, TargetDataSet, Variables, FlatPrecision = 500):
    #     assert Variables, 'P2VV-INFO: Speciafy variables.'

    def makeFlat(self, NumBins):

        from MultiDimMatching import UniformTransformation

        isWeighted = lambda distr: self._weighted[distr]

        print 'P2VV-INFO: Building flat transformations.'
        self._sourceDataFlatTransformations = { var : UniformTransformation( self._sourceData[var], Weights=isWeighted('Source'), Bins = NumBins ) for var in self._variables }
        self._targetDataFlatTransformations = { var : UniformTransformation( self._targetData[var], Weights=isWeighted('Target'), Bins = NumBins ) for var in self._variables }

        print 'P2VV-INFO: Transforming data to flat.'
        source_data = lambda var: self._sourceData[var] if not isWeighted('Source') else map(lambda x: x[0], self._sourceData[var])
        target_data = lambda var: self._targetData[var] if not isWeighted('Target') else map(lambda x: x[0], self._targetData[var])


        self._sourceDataFlat = { var : self._sourceDataFlatTransformations[var].makeFlat(source_data(var)) for var in self._variables }
        self._targetDataFlat = { var : self._targetDataFlatTransformations[var].makeFlat(target_data(var)) for var in self._variables }


    def makeInvFlat(self):
        print 'P2VV-INFO: Transforming correlated flat source data to target data.'

        self._sourceDataTar = { var : self._targetDataFlatTransformations[var].inverseFlat(self._sourceDataFlatTar[var]) for var in self._variables }

    def gaussTransformation(self, Mean, Sigma):
        if self._rootLib:
            from ROOT import TMath
            erfinv = TMath.ErfInverse
            sqrt   = TMath.Sqrt
        else:
            from scipy.special import erfinv
            from math import sqrt

        return lambda x: Mean + sqrt(2) * Sigma * erfinv(2 * x - 1)

    def gaussInvTransformation(self, Mean, Sigma):
        if self._rootLib:
            from ROOT import TMath
            erfc = TMath.Erfc
            sqrt   = TMath.Sqrt
        else:
            from scipy.special import erfc
            from math import sqrt

        return lambda x: erfc( -(x-Mean) / (sqrt(2)*Sigma) ) / 2.

    def makeGauss(self, mean, sigma):

        gauss = self.gaussTransformation(mean, sigma)

        print 'P2VV-INFO: Transforming data to normal.'
        self._sourceDataNormal = { var : map(gauss,self._sourceDataFlat[var]) for var in self._variables }
        self._targetDataNormal = { var : map(gauss,self._targetDataFlat[var]) for var in self._variables }

    def makeInvGauss(self, mean, sigma):

        invGauss = self.gaussInvTransformation(mean, sigma)

        print 'P2VV-INFO: Transforming correlated normal source data to flat target data.'
        self._sourceDataFlatTar = { var : map(invGauss,self._sourceDataNormalTarCorr[var]) for var in self._variables }

    def computeVariance(self):
        self._source_size = float( len( self._sourceData[self._variables[0]] ) )
        self._target_size = float( len( self._targetData[self._variables[0]] ) )

        # compute mean
        print 'P2VV-INFO: Computing mean'
        self._sourceDataMean = { var : sum(self._sourceDataNormal[var]) / self._source_size for var in self._variables }
        self._targetDataMean = { var : sum(self._targetDataNormal[var]) / self._target_size for var in self._variables }

        # compute variance
        print 'P2VV-INFO: Computing covariance'
        self._sourceDataVariance = { var : {} for var in self._variables }
        self._targetDataVariance = { var : {} for var in self._variables }
        for var1 in self._variables:
            for var2 in self._variables:
                if var1 in self._sourceDataVariance[var2].keys() and var1 in self._targetDataVariance[var2].keys():
                    self._sourceDataVariance[var1][var2] = self._sourceDataVariance[var2][var1]
                    self._targetDataVariance[var1][var2] = self._targetDataVariance[var2][var1]

                var1_values_src, var1_mean_src = self._sourceDataNormal[var1], self._sourceDataMean[var1]
                var2_values_src, var2_mean_src = self._sourceDataNormal[var2], self._sourceDataMean[var2]
                var1_values_tar, var1_mean_tar = self._targetDataNormal[var1], self._targetDataMean[var1]
                var2_values_tar, var2_mean_tar = self._targetDataNormal[var2], self._targetDataMean[var2]

                source_covariance = lambda tpl: (tpl[0]-var1_mean_src) * (tpl[1]-var2_mean_src)
                target_covariance = lambda tpl: (tpl[0]-var1_mean_tar) * (tpl[1]-var2_mean_tar)

                self._sourceDataVariance[var1][var2] = map(source_covariance, zip(var1_values_src,var2_values_src))
                self._targetDataVariance[var1][var2] = map(target_covariance, zip(var1_values_tar,var2_values_tar))

    def averageCovariance(self):
        from ROOT import TMatrixD
        dim = self._dimention
        self._covMatrix = { 'source':TMatrixD(dim,dim), 'target':TMatrixD(dim,dim) }

        self._sourceDataVariance_Av = { var : {} for var in self._variables }
        self._targetDataVariance_Av = { var : {} for var in self._variables }
        for idx1, var1 in enumerate(self._variables):
            for idx2, var2 in enumerate(self._variables):

                source_cov = sum(self._sourceDataVariance[var1][var2]) / self._source_size
                target_cov = sum(self._targetDataVariance[var1][var2]) / self._target_size

                self._covMatrix['source'][idx1][idx2] = source_cov
                self._covMatrix['target'][idx1][idx2] = target_cov

                self._sourceDataVariance_Av[var1][var2] = source_cov
                self._targetDataVariance_Av[var1][var2] = target_cov

    def diagVals(self, Row, TransformationMatrix, TempContainer ):
        for idx, var in enumerate(self._variables):
            TempContainer[idx] = Row[var]
        diagVls = TransformationMatrix * TempContainer

        return { key : val for key, val in zip(self._variables,diagVls) }

    def corrVals(self, Row, InvTransformationMatrix, TempContainer ):
        for idx, var in enumerate(self._variables):
            TempContainer[idx] = Row[var]
        corrVls = InvTransformationMatrix * TempContainer

        return { key : val for key, val in zip(self._variables,corrVls) }

    def decorelate(self):

        self.averageCovariance()

        # diagonalization
        from ROOT import TVectorD, TMatrixD
        dataSize     = { 'source':int(self._source_size), 'target':int(self._target_size) }
        eigenVectors = { 'source' : TVectorD(self._dimention), 'target' : TVectorD(self._dimention)}
        dataNormal   = { 'source': self._sourceDataNormal, 'target':self._targetDataNormal }
        self._invTransMatrix   = {}
        self._transMatrix      = {}
        self._dataNormalUncorr = dict( source = { var : [] for var in self._variables },
                                       target = { var : [] for var in self._variables }
                                       )

        vals = TVectorD(self._dimention) # dummy container
        for specie, eigenVec in eigenVectors.iteritems():

            self._invTransMatrix[specie] = self._covMatrix[specie].EigenVectors(eigenVec)
            self._transMatrix[specie]    = TMatrixD(self._invTransMatrix[specie]).Invert()

            print 'P2VV-INFO: Decorelating %s data.'%specie

            for idx in range(dataSize[specie]):
                # TODO: This is not efficient looping. Use sql or numpy or something else.
                row = {}
                for key in self._variables:
                    row[key] = dataNormal[specie][key][idx]

                diag_vals = self.diagVals(row,self._transMatrix[specie],vals)
                for key in self._variables:
                    self._dataNormalUncorr[specie][key] += [diag_vals[key]]

    def correlateSourceToTarget(self):
        from ROOT import TVectorD

        vals = TVectorD( self._dimention )
        print 'P2VV-INFO: Transforming uncorrelated normal source data to correlated normal target data.'
        self._sourceDataNormalTarCorr = { var : [] for var in self._variables }
        for idx in range(int(self._source_size)):
            # TODO: This is not efficient looping. Use sql or numpy or something else.
            row = {}
            for key in self._variables:
                row[key] = self._dataNormalUncorr['source'][key][idx]

            corr_vals = self.corrVals(row,self._invTransMatrix['target'],vals)
            for key in self._variables:
                self._sourceDataNormalTarCorr[key] += [corr_vals[key]]

    def writeToTree(self,fileName):
        print 'P2VV-INFO: Writting branches to input source tree.'

        from ROOT import TFile, addSWeightToTree
        outFile = TFile.Open(fileName, 'recreate')

        # helping stuff
        inData      = dict( source = self._SourceDataSet,
                            target = self._TargetDataSet )
        outBranches = dict( source = {'flat':self._sourceDataFlat, 'normal':self._sourceDataNormal, 'uncorr':self._dataNormalUncorr['source']},
                            target = {'flat':self._targetDataFlat, 'normal':self._targetDataNormal, 'uncorr':self._dataNormalUncorr['target']} )

        if self._verbose:
            outBranches['source'].update( {'corr_norm':self._sourceDataNormalTarCorr, 'corr_flat':self._sourceDataFlatTar, 'matched':self._sourceDataTar})

        from array import array
        _ar = lambda l: array('d',l)
        for specie, input_data in inData.iteritems():

            directory = outFile.mkdir(specie)
            directory.cd()
            outTree = input_data.CloneTree()

            for stage, branches in outBranches[specie].iteritems():
                for branch_name, branch in branches.iteritems():
                    addSWeightToTree(_ar(branch), len(branch), outTree, '%s_%s'%(branch_name,stage))

            directory.cd()
            outTree.Write()
        outFile.Close()

    def match(self, **kwargs):

        self.makeFlat(self._nBins)
        self.makeGauss(self._mean, self._sigma)
        self.computeVariance()
        self.decorelate()
        self.correlateSourceToTarget()
        self.makeInvGauss(self._mean, self._sigma)
        self.makeInvFlat()

        if '.root' in self._outFile:
            self.writeToTree(self._outFile)
        else:
            print 'Do not know how to write to the provided file type.'

class UniformTransformation:
    def __init__(self, var, **kwargs):
        weights = kwargs.pop('Weights', False)
        nBins   = kwargs.pop('Bins', 500)
        posDef  = kwargs.pop('PosDefData', False)

        assert type(var) == list, 'Input variable array type must be %s'%list
        assert type(nBins) == int or type(nBins) == list, 'Bins type must be either int or array of bin edges'
        if weights:
            assert type(var[0]) == tuple and len(var[0]) == 2, \
            'You are trying to pass a weighted distriubution in the wrong format. Correct format is a list of (variable,weight).'

        # assending order of input
        if not weights:
            var = map(lambda x: (x,1), var)
        var.sort( key = lambda x: x[0] )

        from numpy import array as nparray
        self._var   = nparray(var)
        self._nBins = nBins
        self._oFlow = 0
        self._uFlow = 0
        self._isPos = posDef

        # build cumulative distribution
        self.buildCumulative()

    def buildCumulative(self):
        from numpy import histogram, cumsum, array as nparray

        print 'Building cumulative distribution function'
        data    = map(lambda x: x[0], self._var)
        weights = map(lambda x: x[1], self._var)
        hist, bins = histogram( data, weights = weights, bins = self._nBins, density = False )

        sumW = float( sum(weights) )
        cdf  = nparray( map(lambda c: float(c) / sumW, cumsum(hist)) )

        # Assign cdf values to bin centers
        b_width = lambda i: float(bins[i+1] - bins[i]) / 2.
        self._cdf_x = { i + 1 : int(bins[i] + b_width(i)) for i in xrange(len(bins)-1) }
        self._cdf_y = { i + 1 : c for i,c in enumerate(cdf) }

        # include fake superflows
        from math import ceil
        self._cdf_x.update( { 0 : int(min(data)) } )
        self._cdf_y.update( { 0 : 1. / sumW } )

        self._cdf_x.update( { max(self._cdf_x.keys()) : max(ceil(max(data)), max(self._cdf_x.values()) ) } )
        self._cdf_y.update( { max(self._cdf_y.keys()) : 1. } )

        # cdf accessors
        self._x1 = lambda i: self._cdf_x[i-1]
        self._y1 = lambda i: self._cdf_y[i-1]
        self._x2 = lambda i: self._cdf_x[i]
        self._y2 = lambda i: self._cdf_y[i]

    def _interpolate(self, tpl):
        x1, x2 = tpl[0], tpl[1]
        y1, y2 = tpl[2], tpl[3]
        value = tpl[4]

        try:
            dydx = float(y2-y1) / float(x2-x1)
        except ZeroDivisionError:
            if y2-y1 == 0:
                dydx = 0.
            else:
                raise ValueError('Derivative is not defined')

        return y1 + dydx * (value-x1)

    def _digitize(self,data,bins,accessors):

        x1,x2 = accessors['x1'],accessors['x2']
        y1,y2 = accessors['y1'],accessors['y2']

        from numpy import digitize
        indices = digitize(data,bins)
        try:
            points  = [ ( x1(i), x2(i), y1(i), y2(i), val) for i, val in zip(indices,data) ]
        except:
            import pdb; pdb.set_trace()
            raise ValueError('Digitization failed. There is a problem with the accesors')

        check = lambda tpl: tpl[0] <= tpl[4] and tpl[1] >= tpl[4] and tpl[2] <= tpl[3]
        assert len(filter(check,points)) == len(data), 'Bad digitization'

        return points

    def makeFlat(self,data):

        assert hasattr(data,'__getitem__'), 'Input data is not iterable'

        accessors = { 'x1':self._x1, 'x2':self._x2, 'y1':self._y1, 'y2':self._y2 }

        points  = self._digitize(data, self._cdf_x.values(), accessors)
        mapped  = map(self._interpolate, points)

        # sanity check
        mapping_check  = lambda x: x >= 0 and x<=1
        assert len(filter(mapping_check, mapped)) == len(data), 'Not all data are mapped properly'

        return mapped

    def inverseFlat(self,data):

        assert hasattr(data,'__getitem__'), 'Input data is not iterable'
        assert min(data) >= 0 and max(data) <= 1, 'Data are not bounded within [0,1]'

        # clean data from superflows (relevant when inverseFlat is applied to a different set)
        l_bound, h_bound  = self._cdf_y[1], self._cdf_y[max(self._cdf_y.keys())-1]
        overF   = filter(lambda x: x > h_bound, data)
        undrF   = filter(lambda x: x < l_bound, data)
        nrmData = filter(lambda x: x >= l_bound and x <= h_bound, data)

        accessors = { 'x1':self._y1, 'x2':self._y2, 'y1':self._x1, 'y2':self._x2 }

        points  = self._digitize(nrmData, self._cdf_y.values(), accessors)
        invrtd  = map(self._interpolate, points)

        # handle superflows by extrapolating
        bin_edge = 1
        while True:
            tpl = [self._cdf_y[0], self._cdf_y[bin_edge], self._cdf_x[0], self._cdf_x[bin_edge]]
            try:
                invrtd += [ self._interpolate(tuple(tpl + [val])) for val in undrF ]
                break
            except ValueError:
                bin_edge += 1

        maxKey = max(self._cdf_y.keys())
        sbtrct = 1
        while True:
            tpl = [self._cdf_y[maxKey-sbtrct], self._cdf_y[maxKey], self._cdf_x[maxKey-sbtrct], self._cdf_x[maxKey]]
            try:
                invrtd += [ self._interpolate(tuple(tpl + [val])) for val in overF ]
                break
            except ValueError:
                sbtrct +=1

        if len(overF):
            print 'Data contain %s overflows.'%len(overF)
            self._oFlow = len(overF)
        if len(undrF):
            print 'Data contain %s underflows.'%len(undrF)
            self._uFlow = len(undrF)

        # sanity check
        if self._isPos:
            assert len(filter(lambda x: x >= 0, invrtd)) == len(data), 'Bad inversion data is not pos def.'
        assert len(invrtd) == len(data), 'Not all data can be mapped'

        return invrtd
