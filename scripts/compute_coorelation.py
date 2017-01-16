from P2VV.Load import P2VVLibrary

filePath = 'tempFile'
fileName = 'Bs2JpsiKst'

# v1_name = 'muplus_P_gauss'
# v2_name = 'muminus_P_gauss'

v1_name = 'Kplus_P_gauss'
v2_name = 'piminus_P_gauss'

from ROOT import TFile
data = TFile.Open('%s_%s_%s.root'%(filePath,v1_name.strip('_gauss'),v2_name.strip('_gauss'))).Get(fileName)
size = int(data.GetEntries())

# helping stuff
from P2VV import RooFitDecorators
v1_getter = lambda ev: getattr(ev, v1_name)
v2_getter = lambda ev: getattr(ev, v2_name)

# distributions
v1_distr  = map(v1_getter, data)
v2_distr  = map(v2_getter, data)

# mean
v1_mean  = sum(v1_distr) / float(size)
v2_mean  = sum(v2_distr) / float(size)

# variances
from math import sqrt
v1_variance_raw = map(lambda x: (x-v1_mean)**2, v1_distr)
v2_variance_raw = map(lambda x: (x-v2_mean)**2, v2_distr)

v1_error = sqrt( sum(v1_variance_raw) / float(size) )
v2_error = sqrt( sum(v2_variance_raw) / float(size) )

cov_raw = [ (v1_val-v1_mean) * (v2_val-v2_mean) for v1_val, v2_val in zip(v1_distr, v2_distr) ]

cov  = sum(cov_raw) / float(size)
corr = cov / v1_error * v2_error

# decorelate
from ROOT import TMatrixD, TVectorD
covMat = TMatrixD(2,2)
covMat[0][0] = v1_error * v1_error
covMat[1][1] = v2_error * v2_error
covMat[0][1] = cov
covMat[1][0] = cov

eigenVec = TVectorD(2)
invTrans = covMat.EigenVectors(eigenVec)
trans    = TMatrixD(invTrans).Invert()

vals = TVectorD(2)
def diagVals(v1, v2):
    vals[0] = v1
    vals[1] = v2
    diagVls = trans * vals
    return tuple([diagVls[i] for i in range(len(diagVls))])


uncorr_distr = [ diagVals(v1_val, v2_val) for (v1_val, v2_val) in zip(v1_distr,v2_distr) ]
v1_distr_uncorr = map(lambda item: item[0], uncorr_distr)
v2_distr_uncorr = map(lambda item: item[1], uncorr_distr)

# append to tree
from ROOT import addSWeightToTree
from array import array

_ar = lambda l: array('d',l)

fout = TFile.Open('out_%s_%s.root'%(v1_name,v2_name),'recreate')
tout = data.CloneTree()

addSWeightToTree(_ar(v1_distr_uncorr), size, tout, '%s_uncorr'%v1_name)
addSWeightToTree(_ar(v2_distr_uncorr), size, tout, '%s_uncorr'%v2_name)

fout.cd()
tout.Write()
fout.Close()
