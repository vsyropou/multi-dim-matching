import argparse
parser = argparse.ArgumentParser()
parser.add_argument('numBins', type=int)
parser.add_argument('--weights', action='store_true')
opts = parser.parse_args()

# inputs
path   = '/project/bfys/vsyropou/data/Bs2JpsiKst'
name   = 'Bs2JpsiKst'
wgtNam = 'Bs_sWeight'

fileName_src = path + '/MC/Bs_monte_carlo/P2VVDataSet_Bs2JpsiKst_2011_Sim08_Reco14_negKaons_896_931_KpiBin_250615.root'
fileName_trt = path + '/RealData/forIterProcedure/P2VVDataSet_Bs2JpsiKst_2011_Reco14_negKaons_bin3KpiBins_250615_BsWeights.root'

# read data
from ROOT import TFile
inFile_trt     = TFile.Open(fileName_trt)
inFile_src     = TFile.Open(fileName_src)

from P2VV import RooFitDecorators
temp_file = TFile.Open('temp_file.root','recreate')
inTree_trt = inFile_trt.Get(name).buildTree(WeightName = wgtNam)
inTree_src = inFile_src.Get(name).buildTree()

inData_trt = [ (getattr(ev,'Kplus_P'),getattr(ev,wgtNam)) for ev in inTree_trt ]
inData_trt_raw = map(lambda x: x[0], inData_trt)
if not opts.weights:
    inData_trt = map(lambda x: x[0], inData_trt)

inData_src = [ getattr(ev,'Kplus_P') for ev in inTree_src ]

# make binning
print 'Making binning'
class BinWidth:
    def __init__(self, dataset, approx_nbins, scale = 1.2):
        #TODO: extend for weighted datasets

        self._first   = int(min(dataset))
        self._mean    = sum(dataset) / len(dataset)
        self._scale   = scale
        self._aprx_bn = approx_nbins

    def __call__(self,value):
        quantile = int(value) / int(self._mean)
        nBins    = float(self._aprx_bn) / float(self._scale**(quantile + 1))
        return (quantile + 1) * float(self._mean - self._first) / self._aprx_bn

bin_width_src = BinWidth(inData_src,opts.numBins)
bin_width_trt = BinWidth(inData_trt_raw,opts.numBins)

bins_src = [int(min(inData_src))]
bins_trt = [int(min(inData_trt_raw))]

max_inData_src = max(inData_src)
while bins_src[-1] < max_inData_src:
    bins_src += [ round( bins_src[-1] + bin_width_src(bins_src[-1]),2) ]

max_inData_trt = max(inData_trt_raw)
while bins_trt[-1] < max_inData_trt:
    bins_trt += [ round( bins_trt[-1] + bin_width_trt(bins_trt[-1]),2) ]

# build trnsformations
from MultiDimMatching import UniformTransformation
source_trans = UniformTransformation(inData_src, Weights = False,  Bins = bins_src, PosDefData = True)
target_trans = UniformTransformation(inData_trt, Weights = wgtNam if opts.weights else False, Bins = bins_trt, PosDefData = True)

from P2VV.Utilities.MCReweighting import UniFunc
source_trans_old = UniFunc(inData_src, nbinsmax = len(bins_src))
target_trans_old = UniFunc(inData_trt, nbinsmax = len(bins_trt))



# apply transformations
print 'Transforming data my func'
source_data_flat = source_trans.makeFlat(inData_src)
target_data_flat = target_trans.makeFlat(inData_trt_raw)

print 'Transforming data diegos func'
source_data_flat_old = [ source_trans_old(event) for event in inData_src ]
target_data_flat_old = [ target_trans_old(event) for event in inData_trt_raw ]

print 'Matching my func'
target_data_matched = target_trans.inverseFlat(source_data_flat)

print 'Matching diegos func'
target_data_matched_old = [ target_trans_old.inverse(event) for event in source_data_flat_old ]


## create fill and draw histos
from ROOT import TCanvas, TMath, TH1F
from array import array

# helping stuff
l2ar = lambda x: array('d',x)
fill = lambda ev, histogram: histogram.Fill( ev )
canv = lambda name: TCanvas(name,name,600,400)

ks_test = lambda d1, d2: TMath.KolmogorovTest(len(d1),l2ar(d1),len(d2),l2ar(d2),'D')

# create
Min = min(min(target_data_matched),min(target_data_matched))
Max = max(max(target_data_matched),max(target_data_matched))
nBins = 100
histos = dict( source = TH1F('source',  'source',  nBins, 0., 1.),
               target = TH1F('target',  'target',  nBins, 0., 1.),
               matched = TH1F('matched', 'matched', nBins, Min, Max)
               )

histos_comp = dict( source  = TH1F('source_cmp',  'source_cmp',  nBins, 0., 1.),
                    target  = TH1F('target_cmp',  'target_cmp',  nBins, 0., 1.),
                    matched = TH1F('matched_cmp', 'matched_cmp', nBins, Min, Max)
               )


# fill
print 'Filling histograms'
for data, comp_data, info in zip( [source_data_flat,     target_data_flat,     target_data_matched ],
                                  [source_data_flat_old, target_data_flat_old, target_data_matched_old],
                                  ['source', 'target', 'matched']):
    [ fill(ev,histos[info]) for ev in data ]
    [ fill(ev,histos_comp[info]) for ev in comp_data]
    data.sort()
    comp_data.sort()
    print '%s: %s'%(info,ks_test(data,comp_data))

# Draw
stash = []
canv_temp = canv('temp')
for distr_type in histos.iterkeys():

    h  = histos[distr_type]
    hc = histos_comp[distr_type]

    canv_temp.cd()
    h.Draw('err')
    hmn, hmx = h.GetMinimum() * 0.6, h.GetMaximum() * 1.2
    for hist in [h,hc]:
        hist.SetAxisRange(hmn,hmx,'Y')
        hist.SetStats(0)

    h.SetLineColor(2)
    h.SetMarkerColor(2)

    cn = canv(distr_type)
    h.Draw('err')
    hc.Draw('same')

    stash += [cn]


import os
os.system('rm -f %s'%temp_file.GetName())
