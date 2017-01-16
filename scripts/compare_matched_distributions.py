
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( 'path' )
parser.add_argument( '--sourceSufx', default='matched', type=str )
parser.add_argument( '--nBins', default=100, type=int )
parser.add_argument( '--logY', action = 'store_true')
opts = parser.parse_args()

from ROOT import TFile
f = TFile.Open(opts.path)
input_trees = { 'source' : f.Get('source/Bs2JpsiKst'),
                'target' : f.Get('target/Bs2JpsiKst'),
                'source_%s'%opts.sourceSufx : f.Get('source/Bs2JpsiKst')
                }

bin_def = dict( piminus_P = (opts.nBins,0,8e4),
                Kplus_P   = (opts.nBins,0,15e4)
                )

# make histos
from ROOT import TH1F
from math import sqrt
histos = {}
for var, bins in bin_def.iteritems():
    histos[var] = {}
    for tree_name, tree in input_trees.iteritems():

        h_name = '%s_%s'%(tree_name,var)
        hist = TH1F(h_name, h_name, *bins)

        sufx = '_%s'%opts.sourceSufx if opts.sourceSufx in tree_name else ''
        tree.Draw('%s>>%s'%(var + sufx,h_name),'','goff')

        scale = 1. /float(hist.GetEntries())

        for i in range( hist.GetNbinsX() ):
            N = hist.GetBinContent(i+1)
            hist.SetBinContent(i+1, N*scale)
            hist.SetBinError(i+1, sqrt(N) * scale)

        histos[var][tree_name] = hist

# draw histos
from ROOT import TCanvas, gPad
c = TCanvas('comp','comp',600,400)
c.Divide(2,2)
idx = 1
for tree_name in [ ['source','target'], ['source_matched','target'] ]:
    for _idx, var in enumerate(bin_def.keys()):
        # _idx += 1
        c.cd( idx + _idx )

        histos[var][tree_name[0]].SetLineColor(2)
        histos[var][tree_name[0]].SetMarkerColor(2)

        histos[var][tree_name[0]].SetStats(0)
        histos[var][tree_name[1]].SetStats(0)

        histos[var][tree_name[0]].Draw('E')
        histos[var][tree_name[1]].Draw('same E')

        if opts.logY: gPad.SetLogy()


    idx += 2
