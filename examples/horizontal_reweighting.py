
dataPath = '/project/bfys/vsyropou/data/Bs2JpsiKst/MC/Bs_monte_carlo/P2VVDataSet_Bs2JpsiKst_3fb_Sim08_Reco14_negKaons_allKpiBins_250615.root'
tempFileName = 'tempFile'

from ROOT import TFile, gROOT
fmc = TFile.Open(dataPath)
dmc = fmc.Get('Bs2JpsiKst')

# v1_name = 'muplus_P'
# v2_name = 'muminus_P'

v1_name = 'Kplus_P'
v2_name = 'piminus_P'

# read distributions into python lists
from P2VV import RooFitDecorators
# gROOT.cd('PyROOT:/')
ftemp = TFile.Open('%s_%s_%s.root'%(tempFileName,v1_name,v2_name),'recreate')
ftemp.cd()
tmc = dmc.buildTree()

lmc_pK  = []
lmc_ppi = []
for event in tmc:
    lmc_pK  += [ getattr(event,v1_name) ]
    lmc_ppi += [ getattr(event,v2_name) ]

# build transformations
from P2VV.Utilities.MCReweighting import UniFunc
gamma_mcpK  = UniFunc(lmc_pK, nbinsmax=4000)
gamma_mcppi = UniFunc(lmc_ppi, nbinsmax=4000)

# morf distributions
lmc_pK_flat  = map(gamma_mcpK, lmc_pK)
lmc_ppi_flat = map(gamma_mcppi, lmc_ppi)

# append to morfed distributions to initial tree
from ROOT import addSWeightToTree
from array import array

_ar = lambda l: array('d',l)
size = len(lmc_pK_flat)

addSWeightToTree(_ar(lmc_pK_flat), size, tmc, '%s_flat'%v1_name)
addSWeightToTree(_ar(lmc_ppi_flat), size, tmc, '%s_flat'%v2_name)

# morf into gaussian
from ROOT import TMath
gauss = lambda x: TMath.Sqrt(2) * TMath.ErfInverse(2 * x - 1)
# inv_gauss = lambda x: TMath.Erfc(-x/TMath.Sqrt(2)) / 2.

lmc_pK_gauss  = map(gauss, lmc_pK_flat)
lmc_ppi_gauss = map(gauss, lmc_ppi_flat)

addSWeightToTree(_ar(lmc_pK_gauss), size, tmc, '%s_gauss'%v1_name)
addSWeightToTree(_ar(lmc_ppi_gauss), size, tmc, '%s_gauss'%v2_name)

# write
ftemp.cd()
tmc.Write()
ftemp.Close()
