import argparse
parser = argparse.ArgumentParser()
parser.add_argument('numBins', type=int)
parser.add_argument('--inverse', action='store_true')
parser.add_argument('--high_stat', action='store_true')
opts = parser.parse_args()

# inputs
path   = '/project/bfys/vsyropou/data/Bs2JpsiKst'
name   = 'Bs2JpsiKst'
wgtNam = '' #'Bs_sWeight'

fileName_trans = path + '/MC/Bs_monte_carlo/P2VVDataSet_Bs2JpsiKst_2011_Sim08_Reco14_negKaons_896_931_KpiBin_250615.root'
fileName_build = path + '/RealData/forIterProcedure/P2VVDataSet_Bs2JpsiKst_2011_Reco14_negKaons_bin3KpiBins_250615_BsWeights.root'

# read data
from ROOT import TFile
inFile_bld     = TFile.Open(fileName_build)
inFile_trs     = TFile.Open(fileName_trans)

from P2VV import RooFitDecorators
temp_file = TFile.Open('temp_file.root','recreate')
inTree_bld = inFile_bld.Get(name).buildTree(WeightName = wgtNam)
inTree_trs = inFile_trs.Get(name).buildTree()

# inData_bld = [ (getattr(ev,'Kplus_P'),getattr(ev,wgtNam)) for ev in inTree_bld ]
inData_bld = [ (getattr(ev,'Kplus_P'), 1.) for ev in inTree_bld ]
inData_raw_bld = map(lambda x: x[0], inData_bld)
inData_wgt_bld = map(lambda x: x[1], inData_bld)

inData_trs = [ (getattr(ev,'Kplus_P'), 1.) for ev in inTree_trs ]
inData_raw_trs = map(lambda x: x[0], inData_trs)
inData_wgt_trs = map(lambda x: x[1], inData_trs)

# make binning
print 'Making binning'
class BinWidth:
    def __init__(self, dataset, approx_nbins, scale = 1.2):

        self._first   = int(min(dataset))
        self._mean    = sum(dataset) / len(dataset)
        self._scale   = scale
        self._aprx_bn = approx_nbins

    def __call__(self,value):
        quantile = int(value) / int(self._mean)
        nBins    = float(self._aprx_bn) / float(self._scale**(quantile + 1))
        return (quantile + 1) * float(self._mean - self._first) / self._aprx_bn

bin_width = BinWidth(inData_raw_bld,opts.numBins)

bins = [int(min(inData_raw_bld))]
max_inData_raw_bld = max(inData_raw_bld)
while bins[-1] < max_inData_raw_bld:
    bins += [ round( bins[-1] + bin_width(bins[-1]),2) ]

# build trnsformations
data = inData_raw_trs if opts.high_stat else inData_raw_bld

from MultiDimMatching import UniformTransformation
uni_trans = UniformTransformation(data, Weights = wgtNam, Bins = bins, PosDefData = True)

from P2VV.Utilities.MCReweighting import UniFunc
old_trans = UniFunc(data, nbinsmax = len(bins))

# apply transformations
print 'Transforming data'


import time
start = (time.time(), time.clock())
new_data = uni_trans.makeFlat(data)
print 'New data: %s'%(start[0] - time.time())
print 'New data: %s'%(start[1] - time.clock())

start = (time.time(), time.clock())
old_data = [ old_trans(event) for event in data ]
print 'Old data: %s'%(start[0] - time.time())
print 'Old data: %s'%(start[1] - time.clock())

# apply invert transfomation
if opts.inverse:
    print 'Inverting'
    start = (time.time(), time.clock())
    new_data_inv = uni_trans.inverseFlat(new_data)
    print 'New data: %s'%(start[0] - time.time())
    print 'New data: %s'%(start[1] - time.clock())

    start = (time.time(), time.clock())
    old_data_inv = [ old_trans.inverse(event) for event in old_data ]
    print 'Old data: %s'%(start[0] - time.time())
    print 'Old data: %s'%(start[1] - time.clock())


# create a uniform dataset for comparision
from ROOT import RooRandom
uni_def_data = [ RooRandom.uniform() for n in  xrange(len(new_data)) ]
uni_def_data.sort()

## Fill histograms and test
from ROOT import TCanvas, TMath, TH1F
from array import array

l2ar = lambda x: array('d',x)
fill = lambda ev, histogram: histogram.Fill( ev )
canv = lambda name: TCanvas(name,name,600,400)

ks_test = lambda d1, d2: TMath.KolmogorovTest(len(d1),l2ar(d1),len(d2),l2ar(d2),'D')

# create histos
h_new = TH1F('test_new', 'test_new', 100, 0., 1.)
h_old = TH1F('test_old', 'test_old', 100, 0., 1.)
h_ref = TH1F('unif_def', 'unif_def', 100, 0., 1.)

h_init_inv = TH1F('test_init_inv', 'test_init_inv', 100, min(data), max(data))
h_new_inv  = TH1F('test_new_inv', 'test_new_inv', 100, min(data), max(data))
h_old_inv  = TH1F('test_old_inv', 'test_old_inv', 100, min(data), max(data))

print 'Filling histograms'
data.sort()
for dataset, histogram, info in zip( [new_data, old_data, uni_def_data ],
                                     [h_new, h_old, h_ref],
                                     ['uni_new', 'uni_old', None]):

                                     [ fill(ev,histogram) for ev in dataset ]

                                     if info:
                                         dataset.sort()
                                         print '%s: %s'%(info,ks_test(dataset,uni_def_data))

if opts.inverse:
    for dataset, empty_hist, info in zip( [new_data_inv, old_data_inv, data],
                                          [h_new_inv, h_old_inv, h_init_inv],
                                          ['new_inv', 'old_inv', None]
                                          ):

                                          [ fill(event,empty_hist) for event in dataset ]

                                          if info:
                                              dataset.sort()
                                              print '%s: %s'%(info,ks_test(dataset,data))

# Draw
stash = []
canv_temp = canv('temp')
for hn, ho, canvas in zip( [h_new, h_init_inv, h_init_inv],
                           [h_old, h_old_inv, h_new_inv],
                           [canv('trans'),canv('inv_old'),canv('inv_new')] ):
    if not opts.inverse and 'inv' in canvas.GetName():
        del canvas
        continue

    hn.SetLineColor(2)
    hn.SetMarkerColor(2)
    canv_temp.cd()
    ho.Draw('err')

    hmin,hmax = ho.GetMinimum() * 0.6, ho.GetMaximum() * 1.2
    for h in [hn,ho]:
        h.SetAxisRange(hmin,hmax,'Y')
        h.SetStats(0)

    canvas.cd()
    ho.Draw()
    hn.Draw('same err')
    stash += [canvas]

import os
os.system('rm -f %s'%temp_file.GetName())
