import ROOT
import sys
import numpy as np

def averageOverThreshold1p2(numHist, denHist):
  num = 0
  den = 0
  for i in range(numHist.GetNbinsX()+1):
    if (i > 2): 
      num += numHist.GetBinContent(i)
      den += denHist.GetBinContent(i)
  if (den == 0):
    value = 0
  else:
    value = num/den
  return value

def averageOverThreshold1p2Error(numHist, denHist):
  num = 0
  den = 0
  for i in range(numHist.GetNbinsX()+1):
    if (i > 2):
      num += numHist.GetBinContent(i)
      den += denHist.GetBinContent(i)
  if (den == 0):
    value = 0
  else:
    value = np.sqrt((num+den)/den)/den
  return value

def closerMuon(a, b):
  value = 0
  if (a < b):
    value = 1
  else:
    value = 2
  return value


if len(sys.argv) != 3:
  print "USAGE: %s <input file1> <output pdf>"%(sys.argv[0])
  sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]

print "Reading from", inFileName, "and writing to", outFileName

#Flags

regFlag = False
it0Flag = False
it1Flag = False
it2Flag = False

regFlag = True
it0Flag = True
it1Flag = True
it2Flag = True

#Read file

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("test/ntupleTree")

#Define hists

binning = np.array([0, 0.6, 1.2, 1.8, 2.5, 3.5, 5.5, 8.0, 12.0, 30.0])
lengthBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3])
length3DBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3, 20])

numHist = ROOT.TH1D("numHist", "HLT matched tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
numHist.Sumw2()
denHist = ROOT.TH1D("denHist", "Tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
denHist.Sumw2()

numEtaHist = ROOT.TH1D("numEtaHist", "HLT matched tracks;#eta;Counts", 10, -2.5, +2.5)
numEtaHist.Sumw2()
denEtaHist = ROOT.TH1D("denEtaHist", "Tracks;#eta;Counts", 10, -2.5, +2.5)
denEtaHist.Sumw2()
denEtaPtCutHist = ROOT.TH1D("denEtaPtCutHist", "Tracks pT>1.2;#eta;Counts", 10, -2.5, +2.5)
denEtaPtCutHist.Sumw2()

numPhiHist = ROOT.TH1D("numPhiHist", "HLT matched tracks;#phi;Counts", 8, -3.2, +3.2)
numPhiHist.Sumw2()
denPhiHist = ROOT.TH1D("denPhiHist", "Tracks;#phi;Counts", 8, -3.2, +3.2)
denPhiHist.Sumw2()
denPhiPtCutHist = ROOT.TH1D("denPhiPtCutHist", "Tracks pT>1.2;#phi;Counts", 8, -3.2, +3.2)
denPhiPtCutHist.Sumw2()

numLen2DHist = ROOT.TH1D("numLen2DHist", "HLT matched tracks;Bs 2D flight length;Counts", len(lengthBinning)-1, lengthBinning)
numLen2DHist.Sumw2()
denLen2DHist = ROOT.TH1D("denLen2DHist", "Tracks;Bs 2D flight length;Counts", len(lengthBinning)-1, lengthBinning)
denLen2DHist.Sumw2()
denLen2DPtCutHist = ROOT.TH1D("denLen2DPtCutHist", "Tracks pT>1.2;Bs 2D flight length;Counts", len(lengthBinning)-1, lengthBinning)
denLen2DPtCutHist.Sumw2()

numLen3DHist = ROOT.TH1D("numLen3DHist", "HLT matched tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
numLen3DHist.Sumw2()
denLen3DHist = ROOT.TH1D("denLen3DHist", "Tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
denLen3DHist.Sumw2()
denLen3DPtCutHist = ROOT.TH1D("denLen3DPtCutHist", "Tracks pT>1.2;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
denLen3DPtCutHist.Sumw2()

trackDeltaRHist = ROOT.TH1D("trackDeltaR", "#DeltaR hlt vs gen;Track #DeltaR;Counts", 20, 0, 0.01)
muonDeltaRHist = ROOT.TH1D("muonDeltaR", "#DeltaR hlt vs gen;Muon #DeltaR;Counts", 20, 0, 0.01)

myPixHist = ROOT.TH1D("myPixHist", "Pixel tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
myPixHist.Sumw2()
myPixEtaHist = ROOT.TH1D("myPixEtaHist", "Pixel tracks;#eta;Counts", 10, -2.5, +2.5)
myPixEtaHist.Sumw2()
myPixLen3DHist = ROOT.TH1D("myPixLen3DHist", "Pixel tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
myPixLen3DHist.Sumw2()

myPixSizeHist = ROOT.TH1D("myPixSizeHist", "Size of pixelTracks collection", 30, 0, 3000)
myPixDxyHist  = ROOT.TH1D("myPixDxyHist",  "pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
myPixDzHist   = ROOT.TH1D("myPixDzHist",   "pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
myPixDetaHist = ROOT.TH1D("myPixDetaHist", "pixelTracks' eta distance from closest mu cand", 50, 0, 2)
myPixDphiHist = ROOT.TH1D("myPixDphiHist", "pixelTracks' phi distance from closest mu cand", 50, 0, 2)
myPixChi2nHist = ROOT.TH1D("myPixChi2nHist", "pixelTracks' reduced chisquare", 10, 0, 10)
myPixNhitsHist = ROOT.TH1D("myPixNhitsHist", "pixelTracks' number of hits", 3, 2.5, 5.5)

myPixDxyDzHist    = ROOT.TH2D("myPixVtxDist",      "pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
myPixDetaDphiHist = ROOT.TH2D("myPixDetaDphiHist", "pixelTracks' angular distance from closest mu cand", 50, 0, 2, 50, 0, 2)

regPixHist = ROOT.TH1D("regPixHist", "Regional pixel tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
regPixHist.Sumw2()
regPixEtaHist = ROOT.TH1D("regPixEtaHist", "Regional pixel tracks;#eta;Counts", 10, -2.5, +2.5)
regPixEtaHist.Sumw2()
regPixLen3DHist = ROOT.TH1D("regPixLen3DHist", "Regional tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
regPixLen3DHist.Sumw2()

regPixSizeHist = ROOT.TH1D("regPixSizeHist", "Size of regpixelTracks collection", 20, 0, 20)
regPixDxyHist  = ROOT.TH1D("regPixDxyHist",  "reg pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
regPixDzHist   = ROOT.TH1D("regPixDzHist",   "reg pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
regPixDetaHist = ROOT.TH1D("regPixDetaHist", "reg pixelTracks' eta distance from closest mu cand", 50, 0, 2)
regPixDphiHist = ROOT.TH1D("regPixDphiHist", "reg pixelTracks' phi distance from closest mu cand", 50, 0, 2)
regPixChi2nHist = ROOT.TH1D("regPixChi2nHist", "reg pixelTracks' reduced chisquare", 10, 0, 10)

regPixDxyDzHist    = ROOT.TH2D("regPixVtxDist",      "reg pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
regPixDetaDphiHist = ROOT.TH2D("regPixDetaDphiHist", "reg pixelTracks' angular distance from closest mu cand", 50, 0, 2, 50, 0, 2)

ctf0Hist = ROOT.TH1D("ctf0Hist", "Ctf0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf0Hist.Sumw2()
ctf0EtaHist = ROOT.TH1D("ctf0EtaHist", "Ctf0 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf0EtaHist.Sumw2()
ctf0Len3DHist = ROOT.TH1D("ctf0Len3DHist", "Ctf0 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
ctf0Len3DHist.Sumw2()

ctf0SizeHist = ROOT.TH1D("ctf0SizeHist", "Size of ctf Tracks collection;Size;# (norm)", 20, -0.5, 19.5)
ctf0DxyHist  = ROOT.TH1D("ctf0DxyHist",  "Tracks' transverse distance from mumuVtx;#DeltaR [cm];# (norm)", 50, 0, 0.1)
ctf0DzHist   = ROOT.TH1D("ctf0DzHist",   "Tracks' longitudinal distance from mumuVtx;#Deltaz [cm];# (norm)", 50, 0, 0.3)
ctf0DetaHist = ROOT.TH1D("ctf0DetaHist", "Tracks' eta distance from closest mu cand;#Delta#eta;# (norm)", 30, 0, 1.5)
ctf0DphiHist = ROOT.TH1D("ctf0DphiHist", "Tracks' phi distance from closest mu cand;#Delta#phi;# (norm)", 30, 0, 1.5)
ctf0Chi2nHist = ROOT.TH1D("ctf0Chi2nHist", "Tracks' reduced chisquare", 10, 0, 10)

ctf0DxyDzHist    = ROOT.TH2D("ctf0VtxDist",      "ctf0Tracks' distance from mumuVtx", 50, 0, 0.1, 50, 0, 0.3)
ctf0DetaDphiHist = ROOT.TH2D("ctf0DetaDphiHist", "ctf0Tracks' angular distance from closest mu cand", 30, 0, 1.5, 30, 0, 1.5)

iter0Hist = ROOT.TH1D("iter0Hist", "Iter0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter0Hist.Sumw2()
iter0EtaHist = ROOT.TH1D("iter0EtaHist", "Iter0 tracks (;#eta;Counts", 10, -2.5, +2.5)
iter0EtaHist.Sumw2()
iter0Len3DHist = ROOT.TH1D("iter0Len3DHist", "Iter0 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
iter0Len3DHist.Sumw2()

iter0SizeHist = ROOT.TH1D("iter0SizeHist", "Size of iter Tracks collection;size;# (norm)", 20, -0.5, 19.5)

ctf1Hist = ROOT.TH1D("ctf1Hist", "Ctf1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf1Hist.Sumw2()
ctf1EtaHist = ROOT.TH1D("ctf1EtaHist", "Ctf1 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf1EtaHist.Sumw2()
ctf1Len3DHist = ROOT.TH1D("ctf1Len3DHist", "Ctf1 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
ctf1Len3DHist.Sumw2()

ctf1SizeHist = ROOT.TH1D("ctf1SizeHist", "Size of ctf1Tracks collection", 20, -0.5, 19.5)
ctf1DxyHist  = ROOT.TH1D("ctf1DxyHist",  "ctf1Tracks' transverse distance from mumuVtx", 50, 0, 0.1)
ctf1DzHist   = ROOT.TH1D("ctf1DzHist",   "ctf1Tracks' longitudinal distance from mumuVtx", 50, 0, 0.3)
ctf1DetaHist = ROOT.TH1D("ctf1DetaHist", "ctf1Tracks' eta distance from closest mu cand", 30, 0, 1.5)
ctf1DphiHist = ROOT.TH1D("ctf1DphiHist", "ctf1Tracks' phi distance from closest mu cand", 30, 0, 1.5)
ctf1Chi2nHist = ROOT.TH1D("ctf1Chi2nHist", "ctf1Tracks' reduced chisquare", 10, 0, 10)

ctf1DxyDzHist    = ROOT.TH2D("ctf1VtxDist",      "ctf1Tracks' distance from mumuVtx", 50, 0, 0.1, 50, 0, 0.3)
ctf1DetaDphiHist = ROOT.TH2D("ctf1DetaDphiHist", "ctf1Tracks' angular distance from closest mu cand", 30, 0, 1.5, 30, 0, 1.5)

iter1Hist = ROOT.TH1D("iter1Hist", "Iter1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter1Hist.Sumw2()
iter1EtaHist = ROOT.TH1D("iter1EtaHist", "Iter1 tracks (;#eta;Counts", 10, -2.5, +2.5)
iter1EtaHist.Sumw2()
iter1Len3DHist = ROOT.TH1D("iter1Len3DHist", "Iter1 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
iter1Len3DHist.Sumw2()

iter1SizeHist = ROOT.TH1D("iter1SizeHist", "Size of iter1Tracks collection", 20, -0.5, 19.5)

ctf2Hist = ROOT.TH1D("ctf2Hist", "Ctf2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf2Hist.Sumw2()
ctf2EtaHist = ROOT.TH1D("ctf2EtaHist", "Ctf2 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf2EtaHist.Sumw2()
ctf2Len3DHist = ROOT.TH1D("ctf2Len3DHist", "Ctf2 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
ctf2Len3DHist.Sumw2()

ctf2SizeHist = ROOT.TH1D("ctf2SizeHist", "Size of ctf2Tracks collection", 20, -0.5, 19.5)
ctf2DxyHist  = ROOT.TH1D("ctf2DxyHist",  "ctf2Tracks' transverse distance from mumuVtx", 50, 0, 0.1)
ctf2DzHist   = ROOT.TH1D("ctf2DzHist",   "ctf2Tracks' longitudinal distance from mumuVtx", 50, 0, 0.3)
ctf2DetaHist = ROOT.TH1D("ctf2DetaHist", "ctf2Tracks' eta distance from closest mu cand", 30, 0, 1.5)
ctf2DphiHist = ROOT.TH1D("ctf2DphiHist", "ctf2Tracks' phi distance from closest mu cand", 30, 0, 1.5)
ctf2Chi2nHist = ROOT.TH1D("ctf2Chi2nHist", "ctf2Tracks' reduced chisquare", 10, 0, 10)

ctf2DxyDzHist    = ROOT.TH2D("ctf2VtxDist",      "ctf2Tracks' distance from mumuVtx", 50, 0, 0.1, 50, 0, 0.3)
ctf2DetaDphiHist = ROOT.TH2D("ctf2DetaDphiHist", "ctf2Tracks' angular distance from closest mu cand", 30, 0, 1.5, 30, 0, 1.5)

iter2Hist = ROOT.TH1D("iter2Hist", "Iter2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter2Hist.Sumw2()
iter2EtaHist = ROOT.TH1D("iter2EtaHist", "Iter2 tracks (;#eta;Counts", 10, -2.5, +2.5)
iter2EtaHist.Sumw2()
iter2Len3DHist = ROOT.TH1D("iter2Len3DHist", "Iter2 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
iter2Len3DHist.Sumw2()

iter2SizeHist = ROOT.TH1D("iter2SizeHist", "Size of iter2Tracks collection", 20, -0.5, 19.5)

fakePtHist = ROOT.TH1D("fakePtHist", "Unmatched trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
fakePtHist.Sumw2()
fakeEtaHist = ROOT.TH1D("fakeEtaHist", "Unmatched trigger tracks;#eta;Counts", 10, -2.5, 2.5)
fakeEtaHist.Sumw2()

trigPtHist = ROOT.TH1D("trigPtHist", "All trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
trigPtHist.Sumw2()
trigEtaHist = ROOT.TH1D("trigEtaHist", "All trigger tracks;#eta;Counts", 10, -2.5, 2.5)
trigEtaHist.Sumw2()

#Fill hists

for entryNum in range(0, tree.GetEntries()):

  tree.GetEntry(entryNum)

  #Get variables

  hltJpsiTrkMatched = getattr(tree, "gpuMatch_Track")
  hltJpsiMatched = getattr(tree, "muonsMatched")
  hltTrack1Matched = getattr(tree, "gpuMatch_Track1")
  hltTrack2Matched = getattr(tree, "gpuMatch_Track2")

  hltTrigPtVec  = getattr(tree, "gpuTriggerTrack_pt")
  hltTrigEtaVec = getattr(tree, "gpuTriggerTrack_eta")
  hltFakePtVec  = getattr(tree, "gpuTriggerFake_pt")
  hltFakeEtaVec = getattr(tree, "gpuTriggerFake_eta")

  myPixTrack1Matched = getattr(tree, "myPixMatch_Track1")
  myPixTrack2Matched = getattr(tree, "myPixMatch_Track2")
  regPixTrack1Matched   = getattr(tree, "regMatch_Track1")
  regPixTrack2Matched   = getattr(tree, "regMatch_Track2")
  ctf0Track1Matched  = getattr(tree, "ctf0Match_Track1")
  ctf0Track2Matched  = getattr(tree, "ctf0Match_Track2")
  iter0Track1Matched = getattr(tree, "iter0Match_Track1")
  iter0Track2Matched = getattr(tree, "iter0Match_Track2")
  ctf1Track1Matched  = getattr(tree, "ctf1Match_Track1")
  ctf1Track2Matched  = getattr(tree, "ctf1Match_Track2")
  iter1Track1Matched = getattr(tree, "iter1Match_Track1")
  iter1Track2Matched = getattr(tree, "iter1Match_Track2")
  ctf2Track1Matched  = getattr(tree, "ctf2Match_Track1")
  ctf2Track2Matched  = getattr(tree, "ctf2Match_Track2")
  iter2Track1Matched = getattr(tree, "iter2Match_Track1")
  iter2Track2Matched = getattr(tree, "iter2Match_Track2")

  genTrack1Pt = getattr(tree, "track1_pT")
  genTrack2Pt = getattr(tree, "track2_pT")
  genTrack1eta = getattr(tree, "track1_eta")
  genTrack2eta = getattr(tree, "track2_eta")
  genTrack1phi = getattr(tree, "track1_phi")
  genTrack2phi = getattr(tree, "track2_phi")
  genTrack1Len2D = getattr(tree, "bsFlightLength2D")
  genTrack2Len2D = getattr(tree, "bsFlightLength2D")
  genTrack1Len3D = getattr(tree, "bsFlightLength3D")
  genTrack2Len3D = getattr(tree, "bsFlightLength3D")

  myPixTrack1dxy   = getattr(tree, "myPix_Track1_dxy")
  myPixTrack1dz    = getattr(tree, "myPix_Track1_dz")
  myPixTrack1dEta1 = getattr(tree, "myPix_Track1_dEtaMu1")
  myPixTrack1dPhi1 = getattr(tree, "myPix_Track1_dPhiMu1")
  myPixTrack1dEta2 = getattr(tree, "myPix_Track1_dEtaMu2")
  myPixTrack1dPhi2 = getattr(tree, "myPix_Track1_dPhiMu2")
  myPixTrack1chi2n = getattr(tree, "myPix_Track1_normChi2")
  myPixTrack1nHits = getattr(tree, "myPix_Track1_nHits")

  myPixTrack2dxy   = getattr(tree, "myPix_Track2_dxy")
  myPixTrack2dz    = getattr(tree, "myPix_Track2_dz")
  myPixTrack2dEta1 = getattr(tree, "myPix_Track2_dEtaMu1")
  myPixTrack2dPhi1 = getattr(tree, "myPix_Track2_dPhiMu1")
  myPixTrack2dEta2 = getattr(tree, "myPix_Track2_dEtaMu2")
  myPixTrack2dPhi2 = getattr(tree, "myPix_Track2_dPhiMu2")
  myPixTrack2chi2n = getattr(tree, "myPix_Track2_normChi2")
  myPixTrack2nHits = getattr(tree, "myPix_Track2_nHits")

  dRTrack1Mu1sq = myPixTrack1dEta1**2 + myPixTrack1dPhi1**2
  dRTrack1Mu2sq = myPixTrack1dEta2**2 + myPixTrack1dPhi2**2
  myPixTrack1MuFlag = closerMuon(dRTrack1Mu1sq, dRTrack1Mu2sq)

  dRTrack2Mu1sq = myPixTrack2dEta1**2 + myPixTrack2dPhi1**2
  dRTrack2Mu2sq = myPixTrack2dEta2**2 + myPixTrack2dPhi2**2
  myPixTrack2MuFlag = closerMuon(dRTrack2Mu1sq, dRTrack2Mu2sq)

  myPixSize = getattr(tree, "myPix_collSize")

  regPixTrack1dxy   = getattr(tree, "reg_Track1_dxy")
  regPixTrack1dz    = getattr(tree, "reg_Track1_dz")
  regPixTrack1dEta1 = getattr(tree, "reg_Track1_dEtaMu1")
  regPixTrack1dPhi1 = getattr(tree, "reg_Track1_dPhiMu1")
  regPixTrack1dEta2 = getattr(tree, "reg_Track1_dEtaMu2")
  regPixTrack1dPhi2 = getattr(tree, "reg_Track1_dPhiMu2")
  regPixTrack1chi2n = getattr(tree, "reg_Track1_normChi2")

  regPixTrack2dxy   = getattr(tree, "reg_Track2_dxy")
  regPixTrack2dz    = getattr(tree, "reg_Track2_dz")
  regPixTrack2dEta1 = getattr(tree, "reg_Track2_dEtaMu1")
  regPixTrack2dPhi1 = getattr(tree, "reg_Track2_dPhiMu1")
  regPixTrack2dEta2 = getattr(tree, "reg_Track2_dEtaMu2")
  regPixTrack2dPhi2 = getattr(tree, "reg_Track2_dPhiMu2")
  regPixTrack2chi2n = getattr(tree, "reg_Track2_normChi2")

  dRTrack1Mu1sq = regPixTrack1dEta1**2 + regPixTrack1dPhi1**2
  dRTrack1Mu2sq = regPixTrack1dEta2**2 + regPixTrack1dPhi2**2
  regPixTrack1MuFlag = closerMuon(dRTrack1Mu1sq, dRTrack1Mu2sq)

  dRTrack2Mu1sq = regPixTrack2dEta1**2 + regPixTrack2dPhi1**2
  dRTrack2Mu2sq = regPixTrack2dEta2**2 + regPixTrack2dPhi2**2
  regPixTrack2MuFlag = closerMuon(dRTrack2Mu1sq, dRTrack1Mu2sq)

  regPixSize = getattr(tree, "reg_collSize")

  ctf0Track1dxy   = getattr(tree, "ctf0_Track1_dxy")
  ctf0Track1dz    = getattr(tree, "ctf0_Track1_dz")
  ctf0Track1dEta1 = getattr(tree, "ctf0_Track1_dEtaMu1")
  ctf0Track1dPhi1 = getattr(tree, "ctf0_Track1_dPhiMu1")
  ctf0Track1dEta2 = getattr(tree, "ctf0_Track1_dEtaMu2")
  ctf0Track1dPhi2 = getattr(tree, "ctf0_Track1_dPhiMu2")
  ctf0Track1chi2n = getattr(tree, "ctf0_Track1_normChi2")

  ctf0Track2dxy   = getattr(tree, "ctf0_Track2_dxy")
  ctf0Track2dz    = getattr(tree, "ctf0_Track2_dz")
  ctf0Track2dEta1 = getattr(tree, "ctf0_Track2_dEtaMu1")
  ctf0Track2dPhi1 = getattr(tree, "ctf0_Track2_dPhiMu1")
  ctf0Track2dEta2 = getattr(tree, "ctf0_Track2_dEtaMu2")
  ctf0Track2dPhi2 = getattr(tree, "ctf0_Track2_dPhiMu2")
  ctf0Track2chi2n = getattr(tree, "ctf0_Track2_normChi2")

  dRTrack1Mu1sq = ctf0Track1dEta1**2 + ctf0Track1dPhi1**2
  dRTrack1Mu2sq = ctf0Track1dEta2**2 + ctf0Track1dPhi2**2
  ctf0Track1MuFlag = closerMuon(dRTrack1Mu1sq, dRTrack1Mu2sq)

  dRTrack2Mu1sq = ctf0Track2dEta1**2 + ctf0Track2dPhi1**2
  dRTrack2Mu2sq = ctf0Track2dEta2**2 + ctf0Track2dPhi2**2
  ctf0Track2MuFlag = closerMuon(dRTrack2Mu1sq, dRTrack1Mu2sq)
  
  ctf0Size = getattr(tree, "ctf0_collSize")

  iter0Size = getattr(tree, "iter0_collSize")

  ctf1Track1dxy   = getattr(tree, "ctf1_Track1_dxy")
  ctf1Track1dz    = getattr(tree, "ctf1_Track1_dz")
  ctf1Track1dEta1 = getattr(tree, "ctf1_Track1_dEtaMu1")
  ctf1Track1dPhi1 = getattr(tree, "ctf1_Track1_dPhiMu1")
  ctf1Track1dEta2 = getattr(tree, "ctf1_Track1_dEtaMu2")
  ctf1Track1dPhi2 = getattr(tree, "ctf1_Track1_dPhiMu2")
  ctf1Track1chi2n = getattr(tree, "ctf1_Track1_normChi2")

  ctf1Track2dxy   = getattr(tree, "ctf1_Track2_dxy")
  ctf1Track2dz    = getattr(tree, "ctf1_Track2_dz")
  ctf1Track2dEta1 = getattr(tree, "ctf1_Track2_dEtaMu1")
  ctf1Track2dPhi1 = getattr(tree, "ctf1_Track2_dPhiMu1")
  ctf1Track2dEta2 = getattr(tree, "ctf1_Track2_dEtaMu2")
  ctf1Track2dPhi2 = getattr(tree, "ctf1_Track2_dPhiMu2")
  ctf1Track2chi2n = getattr(tree, "ctf1_Track2_normChi2")

  dRTrack1Mu1sq = ctf1Track1dEta1**2 + ctf1Track1dPhi1**2
  dRTrack1Mu2sq = ctf1Track1dEta2**2 + ctf1Track1dPhi2**2
  ctf1Track1MuFlag = closerMuon(dRTrack1Mu1sq, dRTrack1Mu2sq)

  dRTrack2Mu1sq = ctf1Track2dEta1**2 + ctf1Track2dPhi1**2
  dRTrack2Mu2sq = ctf1Track2dEta2**2 + ctf1Track2dPhi2**2
  ctf1Track2MuFlag = closerMuon(dRTrack2Mu1sq, dRTrack1Mu2sq)

  ctf1Size = getattr(tree, "ctf1_collSize")

  iter1Size = getattr(tree, "iter1_collSize") - iter0Size

  ctf2Track1dxy   = getattr(tree, "ctf2_Track1_dxy")
  ctf2Track1dz    = getattr(tree, "ctf2_Track1_dz")
  ctf2Track1dEta1 = getattr(tree, "ctf2_Track1_dEtaMu1")
  ctf2Track1dPhi1 = getattr(tree, "ctf2_Track1_dPhiMu1")
  ctf2Track1dEta2 = getattr(tree, "ctf2_Track1_dEtaMu2")
  ctf2Track1dPhi2 = getattr(tree, "ctf2_Track1_dPhiMu2")
  ctf2Track1chi2n = getattr(tree, "ctf2_Track1_normChi2")

  ctf2Track2dxy   = getattr(tree, "ctf2_Track2_dxy")
  ctf2Track2dz    = getattr(tree, "ctf2_Track2_dz")
  ctf2Track2dEta1 = getattr(tree, "ctf2_Track2_dEtaMu1")
  ctf2Track2dPhi1 = getattr(tree, "ctf2_Track2_dPhiMu1")
  ctf2Track2dEta2 = getattr(tree, "ctf2_Track2_dEtaMu2")
  ctf2Track2dPhi2 = getattr(tree, "ctf2_Track2_dPhiMu2")
  ctf2Track2chi2n = getattr(tree, "ctf2_Track2_normChi2")

  dRTrack1Mu1sq = ctf2Track1dEta1**2 + ctf2Track1dPhi1**2
  dRTrack1Mu2sq = ctf2Track1dEta2**2 + ctf2Track1dPhi2**2
  ctf2Track1MuFlag = closerMuon(dRTrack1Mu1sq, dRTrack1Mu2sq)

  dRTrack2Mu1sq = ctf2Track2dEta1**2 + ctf2Track2dPhi1**2
  dRTrack2Mu2sq = ctf2Track2dEta2**2 + ctf2Track2dPhi2**2
  ctf2Track2MuFlag = closerMuon(dRTrack2Mu1sq, dRTrack1Mu2sq)

  ctf2Size = getattr(tree, "ctf2_collSize")

  iter2Size = getattr(tree, "iter2_collSize") - (iter0Size + iter1Size)

  #Efficiencies

  if (hltJpsiMatched):
    denHist.Fill(genTrack1Pt)
    denHist.Fill(genTrack2Pt)
    denEtaHist.Fill(genTrack1eta)
    denPhiHist.Fill(genTrack1phi)
    denLen2DHist.Fill(genTrack1Len2D)
    denLen3DHist.Fill(genTrack1Len3D)
    denEtaHist.Fill(genTrack2eta)
    denPhiHist.Fill(genTrack2phi)
    denLen2DHist.Fill(genTrack2Len2D)
    denLen3DHist.Fill(genTrack2Len3D)
    if (genTrack1Pt > 1.2):
      denEtaPtCutHist.Fill(genTrack1eta)
      denPhiPtCutHist.Fill(genTrack1phi)
      denLen2DPtCutHist.Fill(genTrack1Len2D)
      denLen3DPtCutHist.Fill(genTrack1Len3D)
    if (genTrack2Pt > 1.2):
      denEtaPtCutHist.Fill(genTrack2eta)
      denPhiPtCutHist.Fill(genTrack2phi)
      denLen2DPtCutHist.Fill(genTrack2Len2D)
      denLen3DPtCutHist.Fill(genTrack2Len3D)

  if (hltJpsiTrkMatched):
    if (hltTrack1Matched):
      numHist.Fill(genTrack1Pt)
      if (genTrack1Pt > 1.2):
        numEtaHist.Fill(genTrack1eta)
        numPhiHist.Fill(genTrack1phi)
        numLen2DHist.Fill(genTrack1Len2D)
        numLen3DHist.Fill(genTrack1Len3D)
    if (hltTrack2Matched):
      numHist.Fill(genTrack2Pt)
      if (genTrack2Pt > 1.2):
        numEtaHist.Fill(genTrack2eta)
        numPhiHist.Fill(genTrack2phi)
        numLen2DHist.Fill(genTrack2Len2D)
        numLen3DHist.Fill(genTrack2Len3D)

  #Truth-matching

  track1DeltaR = getattr(tree, "hlt_Track1_deltaR")
  track2DeltaR = getattr(tree, "hlt_Track2_deltaR")

  if (hltJpsiTrkMatched):
    if (hltTrack1Matched): 
      trackDeltaRHist.Fill(track1DeltaR)
    if (hltTrack2Matched): 
      trackDeltaRHist.Fill(track2DeltaR)

  muon1DeltaR = getattr(tree, "Mu1_deltaR")
  muon2DeltaR = getattr(tree, "Mu2_deltaR")

  if (hltJpsiMatched):
    muonDeltaRHist.Fill(muon1DeltaR)
    muonDeltaRHist.Fill(muon2DeltaR)

  #New pixelTracks

  if (hltJpsiMatched):
    myPixSizeHist.Fill(myPixSize)
    if (myPixTrack1Matched):
      myPixHist.Fill(genTrack1Pt)
      myPixEtaHist.Fill(genTrack1eta)
      myPixLen3DHist.Fill(genTrack1Len3D)
      myPixDxyHist.Fill(myPixTrack1dxy)
      myPixDzHist.Fill(myPixTrack1dz)
      myPixChi2nHist.Fill(myPixTrack1chi2n)
      myPixDxyDzHist.Fill(myPixTrack1dxy, myPixTrack1dz)
      myPixNhitsHist.Fill(myPixTrack1nHits)
      if (myPixTrack1MuFlag == 1):
        myPixDetaHist.Fill(myPixTrack1dEta1)
        myPixDphiHist.Fill(myPixTrack1dPhi1)
        myPixDetaDphiHist.Fill(myPixTrack1dEta1, myPixTrack1dPhi1)
      if (myPixTrack1MuFlag == 2):
        myPixDetaHist.Fill(myPixTrack1dEta2)
        myPixDphiHist.Fill(myPixTrack1dPhi2)
        myPixDetaDphiHist.Fill(myPixTrack1dEta2, myPixTrack1dPhi2)
    if (myPixTrack2Matched):
      myPixHist.Fill(genTrack2Pt)
      myPixEtaHist.Fill(genTrack2eta)
      myPixLen3DHist.Fill(genTrack2Len3D)
      myPixDxyHist.Fill(myPixTrack2dxy)
      myPixDzHist.Fill(myPixTrack2dz)
      myPixChi2nHist.Fill(myPixTrack2chi2n)
      myPixDxyDzHist.Fill(myPixTrack2dxy, myPixTrack2dz)
      myPixNhitsHist.Fill(myPixTrack2nHits)
      if (myPixTrack2MuFlag == 1):
        myPixDetaHist.Fill(myPixTrack2dEta1)
        myPixDphiHist.Fill(myPixTrack2dPhi1)
        myPixDetaDphiHist.Fill(myPixTrack2dEta1, myPixTrack2dPhi1)
      if (myPixTrack2MuFlag == 2):
        myPixDetaHist.Fill(myPixTrack2dEta2)
        myPixDphiHist.Fill(myPixTrack2dPhi2)
        myPixDetaDphiHist.Fill(myPixTrack2dEta2, myPixTrack2dPhi2)

  if (hltJpsiMatched):
    regPixSizeHist.Fill(regPixSize)
    if (regPixTrack1Matched):
      regPixHist.Fill(genTrack1Pt)
      regPixEtaHist.Fill(genTrack1eta)
      regPixLen3DHist.Fill(genTrack1Len3D)
      regPixDxyHist.Fill(regPixTrack1dxy)
      regPixDzHist.Fill(regPixTrack1dz)
      regPixChi2nHist.Fill(regPixTrack1chi2n)
      regPixDxyDzHist.Fill(regPixTrack1dxy, regPixTrack1dz)
      if (regPixTrack1MuFlag == 1):
        regPixDetaHist.Fill(regPixTrack1dEta1)
        regPixDphiHist.Fill(regPixTrack1dPhi1)
        regPixDetaDphiHist.Fill(regPixTrack1dEta1, regPixTrack1dPhi1)
      if (regPixTrack1MuFlag == 2):
        regPixDetaHist.Fill(regPixTrack1dEta2)
        regPixDphiHist.Fill(regPixTrack1dPhi2)
        regPixDetaDphiHist.Fill(regPixTrack1dEta2, regPixTrack1dPhi2)
    if (regPixTrack2Matched):
      regPixHist.Fill(genTrack2Pt)
      regPixEtaHist.Fill(genTrack2eta)
      regPixLen3DHist.Fill(genTrack2Len3D)
      regPixDxyHist.Fill(regPixTrack2dxy)
      regPixDzHist.Fill(regPixTrack2dz)
      regPixChi2nHist.Fill(regPixTrack2chi2n)
      regPixDxyDzHist.Fill(regPixTrack2dxy, regPixTrack2dz)
      if (regPixTrack2MuFlag == 1):
        regPixDetaHist.Fill(regPixTrack2dEta1)
        regPixDphiHist.Fill(regPixTrack2dPhi1)
        regPixDetaDphiHist.Fill(regPixTrack2dEta1, regPixTrack2dPhi1)
      if (regPixTrack2MuFlag == 2):
        regPixDetaHist.Fill(regPixTrack2dEta2)
        regPixDphiHist.Fill(regPixTrack2dPhi2)
        regPixDetaDphiHist.Fill(regPixTrack2dEta2, regPixTrack2dPhi2)

  if (hltJpsiMatched):
    ctf0SizeHist.Fill(ctf0Size)
    if (ctf0Track1Matched):
      ctf0Hist.Fill(genTrack1Pt)
      ctf0EtaHist.Fill(genTrack1eta)
      ctf0Len3DHist.Fill(genTrack1Len3D)
      ctf0DxyHist.Fill(ctf0Track1dxy)
      ctf0DzHist.Fill(ctf0Track1dz)
      ctf0Chi2nHist.Fill(ctf0Track1chi2n)
      ctf0DxyDzHist.Fill(ctf0Track1dxy, ctf0Track1dz)
      if (ctf0Track1MuFlag == 1):
        ctf0DetaHist.Fill(ctf0Track1dEta1)
        ctf0DphiHist.Fill(ctf0Track1dPhi1)
        ctf0DetaDphiHist.Fill(ctf0Track1dEta1, ctf0Track1dPhi1)
      if (ctf0Track1MuFlag == 2):
        ctf0DetaHist.Fill(ctf0Track1dEta2)
        ctf0DphiHist.Fill(ctf0Track1dPhi2)
        ctf0DetaDphiHist.Fill(ctf0Track1dEta2, ctf0Track1dPhi2)
    if (ctf0Track2Matched):
      ctf0Hist.Fill(genTrack2Pt)
      ctf0EtaHist.Fill(genTrack2eta)
      ctf0Len3DHist.Fill(genTrack2Len3D)
      ctf0DxyHist.Fill(ctf0Track2dxy)
      ctf0DzHist.Fill(ctf0Track2dz)
      ctf0Chi2nHist.Fill(ctf0Track2chi2n)
      ctf0DxyDzHist.Fill(ctf0Track2dxy, ctf0Track2dz)
      if (ctf0Track2MuFlag == 1):
        ctf0DetaHist.Fill(ctf0Track2dEta1)
        ctf0DphiHist.Fill(ctf0Track2dPhi1)
        ctf0DetaDphiHist.Fill(ctf0Track2dEta1, ctf0Track2dPhi1)
      if (ctf0Track2MuFlag == 2):
        ctf0DetaHist.Fill(ctf0Track2dEta2)
        ctf0DphiHist.Fill(ctf0Track2dPhi2)
        ctf0DetaDphiHist.Fill(ctf0Track2dEta2, ctf0Track2dPhi2)

    iter0SizeHist.Fill(iter0Size)
    if (iter0Track1Matched):
      iter0Hist.Fill(genTrack1Pt)
      iter0EtaHist.Fill(genTrack1eta)
      iter0Len3DHist.Fill(genTrack1Len3D)
    if (iter0Track2Matched):
      iter0Hist.Fill(genTrack2Pt)
      iter0EtaHist.Fill(genTrack2eta)
      iter0Len3DHist.Fill(genTrack2Len3D)

  if (hltJpsiMatched):
    ctf1SizeHist.Fill(ctf1Size)
    if (ctf1Track1Matched):
      ctf1Hist.Fill(genTrack1Pt)
      ctf1EtaHist.Fill(genTrack1eta)
      ctf1Len3DHist.Fill(genTrack1Len3D)
      ctf1DxyHist.Fill(ctf1Track1dxy)
      ctf1DzHist.Fill(ctf1Track1dz)
      ctf1Chi2nHist.Fill(ctf1Track1chi2n)
      ctf1DxyDzHist.Fill(ctf1Track1dxy, ctf1Track1dz)
      if (ctf1Track1MuFlag == 1):
        ctf1DetaHist.Fill(ctf1Track1dEta1)
        ctf1DphiHist.Fill(ctf1Track1dPhi1)
        ctf1DetaDphiHist.Fill(ctf1Track1dEta1, ctf1Track1dPhi1)
      if (ctf1Track1MuFlag == 2):
        ctf1DetaHist.Fill(ctf1Track1dEta2)
        ctf1DphiHist.Fill(ctf1Track1dPhi2)
        ctf1DetaDphiHist.Fill(ctf1Track1dEta2, ctf1Track1dPhi2)
    if (ctf1Track2Matched):
      ctf1Hist.Fill(genTrack2Pt)
      ctf1EtaHist.Fill(genTrack2eta)
      ctf1Len3DHist.Fill(genTrack2Len3D)
      ctf1DxyHist.Fill(ctf1Track2dxy)
      ctf1DzHist.Fill(ctf1Track2dz)
      ctf1Chi2nHist.Fill(ctf1Track2chi2n)
      ctf1DxyDzHist.Fill(ctf1Track2dxy, ctf1Track2dz)
      if (ctf1Track2MuFlag == 1):
        ctf1DetaHist.Fill(ctf1Track2dEta1)
        ctf1DphiHist.Fill(ctf1Track2dPhi1)
        ctf1DetaDphiHist.Fill(ctf1Track2dEta1, ctf1Track2dPhi1)
      if (ctf1Track2MuFlag == 2):
        ctf1DetaHist.Fill(ctf1Track2dEta2)
        ctf1DphiHist.Fill(ctf1Track2dPhi2)
        ctf1DetaDphiHist.Fill(ctf1Track2dEta2, ctf1Track2dPhi2)

    iter1SizeHist.Fill(iter1Size)
    if (iter1Track1Matched):
      iter1Hist.Fill(genTrack1Pt)
      iter1EtaHist.Fill(genTrack1eta)
      iter1Len3DHist.Fill(genTrack1Len3D)
    if (iter1Track2Matched):
      iter1Hist.Fill(genTrack2Pt)
      iter1EtaHist.Fill(genTrack2eta)
      iter1Len3DHist.Fill(genTrack2Len3D)

  if (hltJpsiMatched):
    ctf2SizeHist.Fill(ctf2Size)
    if (ctf2Track1Matched):
      ctf2Hist.Fill(genTrack1Pt)
      ctf2EtaHist.Fill(genTrack1eta)
      ctf2Len3DHist.Fill(genTrack1Len3D)
      ctf2DxyHist.Fill(ctf2Track1dxy)
      ctf2DzHist.Fill(ctf2Track1dz)
      ctf2Chi2nHist.Fill(ctf2Track1chi2n)
      ctf2DxyDzHist.Fill(ctf2Track1dxy, ctf2Track1dz)
      if (ctf2Track1MuFlag == 1):
        ctf2DetaHist.Fill(ctf2Track1dEta1)
        ctf2DphiHist.Fill(ctf2Track1dPhi1)
        ctf2DetaDphiHist.Fill(ctf2Track1dEta1, ctf2Track1dPhi1)
      if (ctf2Track1MuFlag == 2):
        ctf2DetaHist.Fill(ctf2Track1dEta2)
        ctf2DphiHist.Fill(ctf2Track1dPhi2)
        ctf2DetaDphiHist.Fill(ctf2Track1dEta2, ctf2Track1dPhi2)
    if (ctf2Track2Matched):
      ctf2Hist.Fill(genTrack2Pt)
      ctf2EtaHist.Fill(genTrack2eta)
      ctf2Len3DHist.Fill(genTrack2Len3D)
      ctf2DxyHist.Fill(ctf2Track2dxy)
      ctf2DzHist.Fill(ctf2Track2dz)
      ctf2Chi2nHist.Fill(ctf2Track2chi2n)
      ctf2DxyDzHist.Fill(ctf2Track2dxy, ctf2Track2dz)
      if (ctf2Track2MuFlag == 1):
        ctf2DetaHist.Fill(ctf2Track2dEta1)
        ctf2DphiHist.Fill(ctf2Track2dPhi1)
        ctf2DetaDphiHist.Fill(ctf2Track2dEta1, ctf2Track2dPhi1)
      if (ctf2Track2MuFlag == 2):
        ctf2DetaHist.Fill(ctf2Track2dEta2)
        ctf2DphiHist.Fill(ctf2Track2dPhi2)
        ctf2DetaDphiHist.Fill(ctf2Track2dEta2, ctf2Track2dPhi2)

    iter2SizeHist.Fill(iter2Size)
    if (iter2Track1Matched):
      iter2Hist.Fill(genTrack1Pt)
      iter2EtaHist.Fill(genTrack1eta)
      iter2Len3DHist.Fill(genTrack1Len3D)
    if (iter2Track2Matched):
      iter2Hist.Fill(genTrack2Pt)
      iter2EtaHist.Fill(genTrack2eta)
      iter2Len3DHist.Fill(genTrack2Len3D)

  #Fake rate studies

  if (hltJpsiMatched):
    for trk in hltTrigPtVec:
      trigPtHist.Fill(trk)
    for trk in hltTrigEtaVec:
      trigEtaHist.Fill(trk)
    for trk in hltFakePtVec:
      fakePtHist.Fill(trk)
    for trk in hltFakeEtaVec:
      fakeEtaHist.Fill(trk)

#Compute efficiency

denHist.SetYTitle("Efficiency") #so that effHist copies the right title
denEtaHist.SetYTitle("Efficiency")
denLen2DHist.SetYTitle("Efficiency")
denLen3DHist.SetYTitle("Efficiency")

#trigger

effHist = ROOT.TEfficiency(numHist, denHist)
effHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effHist.SetConfidenceLevel(0.68)
effHist.SetNameTitle("effHist", "Track efficiency")

effEtaHist = ROOT.TEfficiency(numEtaHist, denEtaPtCutHist)
effEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effEtaHist.SetConfidenceLevel(0.68)
effEtaHist.SetNameTitle("effEtaHist", "Track efficiency pT>1.2")

effPhiHist = ROOT.TEfficiency(numPhiHist, denPhiPtCutHist)
effPhiHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effPhiHist.SetConfidenceLevel(0.68)
effPhiHist.SetNameTitle("effPhiHist", "Track efficiency pT>1.2")

effLen2DHist = ROOT.TEfficiency(numLen2DHist, denLen2DPtCutHist)
effLen2DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLen2DHist.SetConfidenceLevel(0.68)
effLen2DHist.SetNameTitle("effLen2DHist", "Track efficiency pT>1.2")

effLen3DHist = ROOT.TEfficiency(numLen3DHist, denLen3DPtCutHist)
effLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLen3DHist.SetConfidenceLevel(0.68)
effLen3DHist.SetNameTitle("effLen3DHist", "Track efficiency pT>1.2")

#new pixel tracks

effMyPixHist = ROOT.TEfficiency(myPixHist, denHist)
effMyPixHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effMyPixHist.SetConfidenceLevel(0.68)
effMyPixHist.SetNameTitle("effMyPixHist", "pixel Tracks efficiency")

effMyPixEtaHist = ROOT.TEfficiency(myPixEtaHist, denEtaHist)
effMyPixEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effMyPixEtaHist.SetConfidenceLevel(0.68)
effMyPixEtaHist.SetNameTitle("effMyPixEtaHist", "pixel Tracks efficiency")

effMyPixLen3DHist = ROOT.TEfficiency(myPixLen3DHist, denLen3DHist)
effMyPixLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effMyPixLen3DHist.SetConfidenceLevel(0.68)
effMyPixLen3DHist.SetNameTitle("effMyPixLen3DHist", "pixel Tracks efficiency")

#regional pixel tracks

effRegPixHist = ROOT.TEfficiency(regPixHist, denHist)
effRegPixHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effRegPixHist.SetConfidenceLevel(0.68)
effRegPixHist.SetNameTitle("effRegPixHist", "reg pixel Tracks efficiency")

effRegPixEtaHist = ROOT.TEfficiency(regPixEtaHist, denEtaHist)
effRegPixEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effRegPixEtaHist.SetConfidenceLevel(0.68)
effRegPixEtaHist.SetNameTitle("effRegPixEtaHist", "reg pixel Tracks efficiency")

effRegPixLen3DHist = ROOT.TEfficiency(regPixLen3DHist, denLen3DHist)
effRegPixLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effRegPixLen3DHist.SetConfidenceLevel(0.68)
effRegPixLen3DHist.SetNameTitle("effRegPixLen3DHist", "reg pixel Tracks efficiency")

#ctf0

effCtf0Hist = ROOT.TEfficiency(ctf0Hist, denHist)
effCtf0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf0Hist.SetConfidenceLevel(0.68)
effCtf0Hist.SetNameTitle("effCtf0Hist", "ctf0 Track efficiency")

effCtf0EtaHist = ROOT.TEfficiency(ctf0EtaHist, denEtaHist)
effCtf0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf0EtaHist.SetConfidenceLevel(0.68)
effCtf0EtaHist.SetNameTitle("effCtf0EtaHist", "ctf0 Track efficiency")

effCtf0Len3DHist = ROOT.TEfficiency(ctf0Len3DHist, denLen3DHist)
effCtf0Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf0Len3DHist.SetConfidenceLevel(0.68)
effCtf0Len3DHist.SetNameTitle("effCtf0Len3DHist", "pixel Tracks efficiency")

#iter0

effIter0Hist = ROOT.TEfficiency(iter0Hist, denHist)
effIter0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter0Hist.SetConfidenceLevel(0.68)
effIter0Hist.SetNameTitle("effIter0Hist", "iter0 Track efficiency")

effIter0EtaHist = ROOT.TEfficiency(iter0EtaHist, denEtaHist)
effIter0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter0EtaHist.SetConfidenceLevel(0.68)
effIter0EtaHist.SetNameTitle("effIter0EtaHist", "iter0 Track efficiency")

effIter0Len3DHist = ROOT.TEfficiency(iter0Len3DHist, denLen3DHist)
effIter0Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter0Len3DHist.SetConfidenceLevel(0.68)
effIter0Len3DHist.SetNameTitle("effIter0Len3DHist", "iter0 Tracks efficiency")

#ctf1

effCtf1Hist = ROOT.TEfficiency(ctf1Hist, denHist)
effCtf1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf1Hist.SetConfidenceLevel(0.68)
effCtf1Hist.SetNameTitle("effCtf1Hist", "ctf1 Track efficiency")

effCtf1EtaHist = ROOT.TEfficiency(ctf1EtaHist, denEtaHist)
effCtf1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf1EtaHist.SetConfidenceLevel(0.68)
effCtf1EtaHist.SetNameTitle("effCtf1EtaHist", "ctf1 Track efficiency")

effCtf1Len3DHist = ROOT.TEfficiency(ctf1Len3DHist, denLen3DHist)
effCtf1Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf1Len3DHist.SetConfidenceLevel(0.68)
effCtf1Len3DHist.SetNameTitle("effCtf1Len3DHist", "ctf1 Tracks efficiency")

#iter1

effIter1Hist = ROOT.TEfficiency(iter1Hist, denHist)
effIter1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter1Hist.SetConfidenceLevel(0.68)
effIter1Hist.SetNameTitle("effIter1Hist", "iter1 Track efficiency")

effIter1EtaHist = ROOT.TEfficiency(iter1EtaHist, denEtaHist)
effIter1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter1EtaHist.SetConfidenceLevel(0.68)
effIter1EtaHist.SetNameTitle("effIter1EtaHist", "iter1 Track efficiency")

effIter1Len3DHist = ROOT.TEfficiency(iter1Len3DHist, denLen3DHist)
effIter1Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter1Len3DHist.SetConfidenceLevel(0.68)
effIter1Len3DHist.SetNameTitle("effIter1Len3DHist", "iter1 Tracks efficiency")

#ctf2

effCtf2Hist = ROOT.TEfficiency(ctf2Hist, denHist)
effCtf2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf2Hist.SetConfidenceLevel(0.68)
effCtf2Hist.SetNameTitle("effCtf2Hist", "ctf2 Track efficiency")

effCtf2EtaHist = ROOT.TEfficiency(ctf2EtaHist, denEtaHist)
effCtf2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf2EtaHist.SetConfidenceLevel(0.68)
effCtf2EtaHist.SetNameTitle("effCtf2EtaHist", "ctf2 Track efficiency")

effCtf2Len3DHist = ROOT.TEfficiency(ctf2Len3DHist, denLen3DHist)
effCtf2Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf2Len3DHist.SetConfidenceLevel(0.68)
effCtf2Len3DHist.SetNameTitle("effCtf2Len3DHist", "ctf2 Tracks efficiency")

#iter2

effIter2Hist = ROOT.TEfficiency(iter2Hist, denHist)
effIter2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter2Hist.SetConfidenceLevel(0.68)
effIter2Hist.SetNameTitle("effIter2Hist", "iter2 Track efficiency")

effIter2EtaHist = ROOT.TEfficiency(iter2EtaHist, denEtaHist)
effIter2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter2EtaHist.SetConfidenceLevel(0.68)
effIter2EtaHist.SetNameTitle("effIter2EtaHist", "iter2 Track efficiency")

effIter2Len3DHist = ROOT.TEfficiency(iter2Len3DHist, denLen3DHist)
effIter2Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter2Len3DHist.SetConfidenceLevel(0.68)
effIter2Len3DHist.SetNameTitle("effIter2Len3DHist", "iter2 Tracks efficiency")

#fake rate

fakePtHist.Divide(trigPtHist)
fakeEtaHist.Divide(trigEtaHist)

#eff by steps

denHist.SetYTitle("Counts")
denEtaHist.SetYTitle("Counts")
denLen2DHist.SetYTitle("Counts")
denLen3DHist.SetYTitle("Counts")

#Compute average efficiency

print "Average efficiency for", effHist.GetName(), ":", averageOverThreshold1p2(numHist, denHist), "+/-", averageOverThreshold1p2Error(numHist, denHist)

xIter = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
xIterErr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

yEff = np.array([averageOverThreshold1p2(myPixHist, denHist), averageOverThreshold1p2(regPixHist, denHist), averageOverThreshold1p2(iter0Hist, denHist), averageOverThreshold1p2(iter1Hist, denHist), averageOverThreshold1p2(iter2Hist, denHist), averageOverThreshold1p2(numHist, denHist)])

yEffError =  np.array([averageOverThreshold1p2Error(myPixHist, denHist), averageOverThreshold1p2Error(regPixHist, denHist), averageOverThreshold1p2Error(iter0Hist, denHist), averageOverThreshold1p2Error(iter1Hist, denHist), averageOverThreshold1p2Error(iter2Hist, denHist), averageOverThreshold1p2Error(numHist, denHist)])

print(yEff[1]/yEff[0])

effOverItersGraph = ROOT.TGraphErrors(6, xIter, yEff, xIterErr, yEffError)
effOverItersGraph.SetNameTitle("effOverIters", "track efficiency at trigger steps")
effOverItersGraph.GetXaxis().SetTitle("Steps")
effOverItersGraph.GetYaxis().SetTitle("Efficiency")

#Set directories

effHist.SetDirectory(0)
effEtaHist.SetDirectory(0)
effPhiHist.SetDirectory(0)
effLen2DHist.SetDirectory(0)
effLen3DHist.SetDirectory(0)

effMyPixHist.SetDirectory(0)
effMyPixEtaHist.SetDirectory(0)
effMyPixLen3DHist.SetDirectory(0)

myPixSizeHist.SetDirectory(0)
myPixDxyHist.SetDirectory(0)
myPixDzHist.SetDirectory(0)
myPixDetaHist.SetDirectory(0)
myPixDphiHist.SetDirectory(0)
myPixChi2nHist.SetDirectory(0)
myPixDxyDzHist.SetDirectory(0)
myPixDetaDphiHist.SetDirectory(0)
myPixNhitsHist.SetDirectory(0)

effRegPixHist.SetDirectory(0)
effRegPixEtaHist.SetDirectory(0)
effRegPixLen3DHist.SetDirectory(0)

regPixSizeHist.SetDirectory(0)
regPixDxyHist.SetDirectory(0)
regPixDzHist.SetDirectory(0)
regPixDetaHist.SetDirectory(0)
regPixDphiHist.SetDirectory(0)
regPixChi2nHist.SetDirectory(0)
regPixDxyDzHist.SetDirectory(0)
regPixDetaDphiHist.SetDirectory(0)

effCtf0Hist.SetDirectory(0)
effCtf0EtaHist.SetDirectory(0)
effCtf0Len3DHist.SetDirectory(0)
effIter0Hist.SetDirectory(0)
effIter0EtaHist.SetDirectory(0)
effIter0Len3DHist.SetDirectory(0)

ctf0SizeHist.SetDirectory(0)
ctf0DxyHist.SetDirectory(0)
ctf0DzHist.SetDirectory(0)
ctf0DetaHist.SetDirectory(0)
ctf0DphiHist.SetDirectory(0)
ctf0Chi2nHist.SetDirectory(0)
ctf0DxyDzHist.SetDirectory(0)
ctf0DetaDphiHist.SetDirectory(0)
iter0SizeHist.SetDirectory(0)

effCtf1Hist.SetDirectory(0)
effCtf1EtaHist.SetDirectory(0)
effCtf1Len3DHist.SetDirectory(0)
effIter1Hist.SetDirectory(0)
effIter1EtaHist.SetDirectory(0)
effIter1Len3DHist.SetDirectory(0)

ctf1SizeHist.SetDirectory(0)
ctf1DxyHist.SetDirectory(0)
ctf1DzHist.SetDirectory(0)
ctf1DetaHist.SetDirectory(0)
ctf1DphiHist.SetDirectory(0)
ctf1Chi2nHist.SetDirectory(0)
ctf1DxyDzHist.SetDirectory(0)
ctf1DetaDphiHist.SetDirectory(0)
iter1SizeHist.SetDirectory(0)

effCtf2Hist.SetDirectory(0)
effCtf2EtaHist.SetDirectory(0)
effCtf2Len3DHist.SetDirectory(0)
effIter2Hist.SetDirectory(0)
effIter2EtaHist.SetDirectory(0)
effIter2Len3DHist.SetDirectory(0)

ctf2SizeHist.SetDirectory(0)
ctf2DxyHist.SetDirectory(0)
ctf2DzHist.SetDirectory(0)
ctf2DetaHist.SetDirectory(0)
ctf2DphiHist.SetDirectory(0)
ctf2Chi2nHist.SetDirectory(0)
ctf2DxyDzHist.SetDirectory(0)
ctf2DetaDphiHist.SetDirectory(0)
iter2SizeHist.SetDirectory(0)

muonDeltaRHist.SetDirectory(0)
trackDeltaRHist.SetDirectory(0)

fakePtHist.SetDirectory(0)
fakeEtaHist.SetDirectory(0)

inFile.Close()

#set graphics

emptyHist = ROOT.TH1D("emptyHist", "title;p_{T} [GeV/c];Efficiency", 10, 0, 30)
emptyEtaHist = ROOT.TH1D("emptyEtaHist", "title;#eta;Efficiency", 10, -2.5, 2.5)
emptyPhiHist = ROOT.TH1D("emptyPhiHist", "title;#phi;Efficiency", 8, -3.2, +3.2)
emptyLen2DHist = ROOT.TH1D("emptyLen2DHist", "title;Bs 2D flight length [cm];Efficiency", len(lengthBinning)-1, lengthBinning)
emptyLen3DHist = ROOT.TH1D("emptyLen3DHist", "title;Bs 3D flight length [cm];Efficiency", len(length3DBinning)-1, length3DBinning)
emptyStepHist = ROOT.TH1D("emptyStepHist", "title;;Efficiency", 6, 0, 6)
emptySizeHist = ROOT.TH1D("emptySizeHist", "title;N tracks;# (norm.)", 20, -0.5, 19.5)

emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyEtaHist.GetYaxis().SetRangeUser(0, 1.05)
emptyPhiHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen2DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen3DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyStepHist.GetYaxis().SetRangeUser(0., 1.05)
emptySizeHist.GetYaxis().SetRangeUser(0., 0.8)

emptyHist.SetStats(ROOT.kFALSE)
emptyEtaHist.SetStats(ROOT.kFALSE)
emptyPhiHist.SetStats(ROOT.kFALSE)
emptyLen2DHist.SetStats(ROOT.kFALSE)
emptyLen3DHist.SetStats(ROOT.kFALSE)
emptyStepHist.SetStats(ROOT.kFALSE)
emptySizeHist.SetStats(ROOT.kFALSE)

idx = 0
binLabels = ["pixeltracks", "regional", "Iter0", "Iter1", "Iter2", "After cuts"]
for label in binLabels:
  emptyStepHist.GetXaxis().SetBinLabel(idx+1, label)
  idx += 1

emptyStepHist.SetLabelSize(0.04)

line = ROOT.TLine(1.2, 0, 1.2, 1.05)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

effHist.SetLineColor(1)
effEtaHist.SetLineColor(1)
effPhiHist.SetLineColor(1)
effLen2DHist.SetLineColor(1)
effLen3DHist.SetLineColor(1)

effMyPixHist.SetLineColor(7)
effMyPixEtaHist.SetLineColor(7)
effMyPixLen3DHist.SetLineColor(7)
effRegPixHist.SetLineColor(6)
effRegPixEtaHist.SetLineColor(6)
effRegPixLen3DHist.SetLineColor(6)
effCtf0Hist.SetLineColor(2)
effCtf0EtaHist.SetLineColor(2)
effCtf0Len3DHist.SetLineColor(2)
effIter0Hist.SetLineColor(3)
effIter0EtaHist.SetLineColor(3)
effIter0Len3DHist.SetLineColor(3)
effCtf1Hist.SetLineColor(4)
effCtf1EtaHist.SetLineColor(4)
effCtf1Len3DHist.SetLineColor(4)
effIter1Hist.SetLineColor(5)
effIter1EtaHist.SetLineColor(5)
effIter1Len3DHist.SetLineColor(5)
effCtf2Hist.SetLineColor(8)
effCtf2EtaHist.SetLineColor(8)
effCtf2Len3DHist.SetLineColor(8)
effIter2Hist.SetLineColor(9)
effIter2EtaHist.SetLineColor(9)
effIter2Len3DHist.SetLineColor(9)

fakePtHist.SetLineColor(1)
fakeEtaHist.SetLineColor(1)

effHist.SetLineWidth(2)
effEtaHist.SetLineWidth(2)
effPhiHist.SetLineWidth(2)
effLen2DHist.SetLineWidth(2)
effLen3DHist.SetLineWidth(2)

effMyPixHist.SetLineWidth(2)
effMyPixEtaHist.SetLineWidth(2)
effMyPixLen3DHist.SetLineWidth(2)
effRegPixHist.SetLineWidth(2)
effRegPixEtaHist.SetLineWidth(2)
effRegPixLen3DHist.SetLineWidth(2)
effCtf0Hist.SetLineWidth(2)
effCtf0EtaHist.SetLineWidth(2)
effCtf0Len3DHist.SetLineWidth(2)
effIter0Hist.SetLineWidth(2)
effIter0EtaHist.SetLineWidth(2)
effIter0Len3DHist.SetLineWidth(2)
effCtf1Hist.SetLineWidth(2)
effCtf1EtaHist.SetLineWidth(2)
effCtf1Len3DHist.SetLineWidth(2)
effIter1Hist.SetLineWidth(2)
effIter1EtaHist.SetLineWidth(2)
effIter1Len3DHist.SetLineWidth(2)
effCtf2Hist.SetLineWidth(2)
effCtf2EtaHist.SetLineWidth(2)
effCtf2Len3DHist.SetLineWidth(2)
effIter2Hist.SetLineWidth(2)
effIter2EtaHist.SetLineWidth(2)
effIter2Len3DHist.SetLineWidth(2)

fakePtHist.SetLineWidth(2)
fakeEtaHist.SetLineWidth(2)

effOverItersGraph.SetLineColor(1)
effOverItersGraph.SetMarkerStyle(2)

myPixSizeHist.GetYaxis().SetTitle("Counts")
myPixDxyHist.GetYaxis().SetTitle("Counts")
myPixDzHist.GetYaxis().SetTitle("Counts")
myPixDetaHist.GetYaxis().SetTitle("Counts")
myPixDphiHist.GetYaxis().SetTitle("Counts")
myPixChi2nHist.GetYaxis().SetTitle("Counts")
myPixDxyDzHist.GetYaxis().SetTitle("#Delta z [cm]")
myPixDetaDphiHist.GetYaxis().SetTitle("#Delta #phi")
myPixNhitsHist.GetYaxis().SetTitle("Counts")

myPixSizeHist.GetXaxis().SetTitle("Track collection's size")
myPixDxyHist.GetXaxis().SetTitle("#Delta xy [cm]")
myPixDzHist.GetXaxis().SetTitle("#Delta z [cm]")
myPixDetaHist.GetXaxis().SetTitle("#Delta #eta")
myPixDphiHist.GetXaxis().SetTitle("#Delta #phi")
myPixChi2nHist.GetXaxis().SetTitle("#chi^{2}")
myPixDxyDzHist.GetXaxis().SetTitle("#Delta xy [cm]")
myPixDetaDphiHist.GetXaxis().SetTitle("#Delta #eta")
myPixNhitsHist.GetXaxis().SetTitle("nHits")

myPixSizeHist.SetLineColor(1)
myPixDxyHist.SetLineColor(1)
myPixDzHist.SetLineColor(1)
myPixDetaHist.SetLineColor(1)
myPixDphiHist.SetLineColor(1)
myPixNhitsHist.SetLineColor(1)

regPixSizeHist.GetYaxis().SetTitle("Counts")
regPixDxyHist.GetYaxis().SetTitle("Counts")
regPixDzHist.GetYaxis().SetTitle("Counts")
regPixDetaHist.GetYaxis().SetTitle("Counts")
regPixDphiHist.GetYaxis().SetTitle("Counts")
regPixChi2nHist.GetYaxis().SetTitle("Counts")
regPixDxyDzHist.GetYaxis().SetTitle("#Delta z [cm]")
regPixDetaDphiHist.GetYaxis().SetTitle("#Delta #phi")

regPixSizeHist.GetXaxis().SetTitle("Track collection's size")
regPixDxyHist.GetXaxis().SetTitle("#Delta xy [cm]")
regPixDzHist.GetXaxis().SetTitle("#Delta z [cm]")
regPixDetaHist.GetXaxis().SetTitle("#Delta #eta")
regPixDphiHist.GetXaxis().SetTitle("#Delta #phi")
regPixChi2nHist.GetXaxis().SetTitle("#chi^{2}")
regPixDxyDzHist.GetXaxis().SetTitle("#Delta xy [cm]")
regPixDetaDphiHist.GetXaxis().SetTitle("#Delta #eta")

fakePtHist.GetYaxis().SetTitle("Fake rate")
fakeEtaHist.GetYaxis().SetTitle("Fake rate")

regPixSizeHist.SetLineColor(2)
regPixDxyHist.SetLineColor(2)
regPixDzHist.SetLineColor(2)
regPixDetaHist.SetLineColor(2)
regPixDphiHist.SetLineColor(2)

ctf0SizeHist.SetLineColor(1)
ctf0DxyHist.SetLineColor(1)
ctf0DzHist.SetLineColor(1)
ctf0DetaHist.SetLineColor(1)
ctf0DphiHist.SetLineColor(1)
ctf1SizeHist.SetLineColor(2)
ctf1DxyHist.SetLineColor(2)
ctf1DzHist.SetLineColor(2)
ctf1DetaHist.SetLineColor(2)
ctf1DphiHist.SetLineColor(2)
ctf2SizeHist.SetLineColor(3)
ctf2DxyHist.SetLineColor(3)
ctf2DzHist.SetLineColor(3)
ctf2DetaHist.SetLineColor(3)
ctf2DphiHist.SetLineColor(3)

iter0SizeHist.SetLineColor(1)
iter1SizeHist.SetLineColor(2)
iter2SizeHist.SetLineColor(3)

regPixSizeHist.SetFillColor(2)
regPixDxyHist.SetFillColor(2)
regPixDzHist.SetFillColor(2)
regPixDetaHist.SetFillColor(2)
regPixDphiHist.SetFillColor(2)

ctf0SizeHist.SetLineWidth(2)
ctf0DxyHist.SetLineWidth(2)
ctf0DzHist.SetLineWidth(2)
ctf0DetaHist.SetLineWidth(2)
ctf0DphiHist.SetLineWidth(2)
ctf1SizeHist.SetLineWidth(2)
ctf1DxyHist.SetLineWidth(2)
ctf1DzHist.SetLineWidth(2)
ctf1DetaHist.SetLineWidth(2)
ctf1DphiHist.SetLineWidth(2)
ctf2SizeHist.SetLineWidth(2)
ctf2DxyHist.SetLineWidth(2)
ctf2DzHist.SetLineWidth(2)
ctf2DetaHist.SetLineWidth(2)
ctf2DphiHist.SetLineWidth(2)

iter0SizeHist.SetLineWidth(2)
iter1SizeHist.SetLineWidth(2)
iter2SizeHist.SetLineWidth(2)

effItersLegend = ROOT.TLegend(0.7, 0.15, 0.88, 0.35)
effItersLegend.SetTextSize(0.03)
effItersLegend.AddEntry(effMyPixHist,  "pixelTracks")
if (regFlag): 
  effItersLegend.AddEntry(effRegPixHist, "reg pixTks")
if (it0Flag):
  effItersLegend.AddEntry(effIter0Hist,  "iter0")
if (it1Flag):
  effItersLegend.AddEntry(effIter1Hist,  "iter1")
if (it2Flag):
  effItersLegend.AddEntry(effIter2Hist,  "iter2")
effItersLegend.AddEntry(effHist,       "afterFilter")

effItersEtaLegend = ROOT.TLegend(0.65, 0.15, 0.83, 0.35)
effItersEtaLegend.SetTextSize(0.03)
effItersEtaLegend.AddEntry(effMyPixEtaHist, "pixelTracks")
if (regFlag):
  effItersEtaLegend.AddEntry(effRegPixEtaHist, "reg pixTks")
if (it0Flag):
  effItersEtaLegend.AddEntry(effIter0EtaHist, "iter0")
if (it1Flag):
  effItersEtaLegend.AddEntry(effIter1EtaHist, "iter1")
if (it2Flag):
  effItersEtaLegend.AddEntry(effIter2EtaHist, "iter2")
effItersEtaLegend.AddEntry(effEtaHist, "afterFilter")

effItersLen3DLegend = ROOT.TLegend(0.7, 0.15, 0.88, 0.35)
effItersLen3DLegend.SetTextSize(0.03)
effItersLen3DLegend.AddEntry(effMyPixLen3DHist,  "pixelTracks")
if (regFlag):
  effItersLen3DLegend.AddEntry(effRegPixLen3DHist, "reg pixTks")
if (it0Flag):
  effItersLen3DLegend.AddEntry(effIter0Len3DHist,  "iter0")
if (it1Flag):
  effItersLen3DLegend.AddEntry(effIter1Len3DHist,  "iter1")
if (it2Flag):
  effItersLen3DLegend.AddEntry(effIter2Len3DHist,  "iter2")
effItersLen3DLegend.AddEntry(effLen3DHist,       "afterFilter")

histLegend = ROOT.TLegend(0.65, 0.63, 0.89, 0.73)
histLegend.SetTextSize(0.03)
histLegend.AddEntry(myPixSizeHist, "all pixelTracks")
histLegend.AddEntry(regPixSizeHist, "reg pixelTracks")

ctfLegend = ROOT.TLegend(0.65, 0.63, 0.85, 0.73)
ctfLegend.SetTextSize(0.03)
if (it0Flag):
  ctfLegend.AddEntry(ctf0SizeHist, "ctf0 tracks")
if (it1Flag):
  ctfLegend.AddEntry(ctf1SizeHist, "ctf1 tracks")
if (it2Flag):
  ctfLegend.AddEntry(ctf2SizeHist, "ctf2 tracks")

iterLegend = ROOT.TLegend(0.65, 0.63, 0.85, 0.73)
iterLegend.SetTextSize(0.03)
if (it0Flag):
  iterLegend.AddEntry(iter0SizeHist, "iter0 tracks")
if (it1Flag):
  iterLegend.AddEntry(iter1SizeHist, "iter1 tracks")
if (it2Flag):
  iterLegend.AddEntry(iter2SizeHist, "iter2 tracks")

#define canvas and print

canvas = ROOT.TCanvas("canvas")
canvas.cd()

canvas.Print(outFileName+"[")

emptyHist.SetTitle(effHist.GetTitle())
emptyHist.Draw("")
effHist.Draw("same")
line.Draw("same")
canvas.Print(outFileName)

emptyEtaHist.SetTitle(effEtaHist.GetTitle())
emptyEtaHist.Draw("")
effEtaHist.Draw("same")
canvas.Print(outFileName)

emptyPhiHist.SetTitle(effPhiHist.GetTitle())
emptyPhiHist.Draw("")
effPhiHist.Draw("same")
canvas.Print(outFileName)

canvas.cd(0).SetLogx(1)
emptyLen2DHist.SetTitle(effLen2DHist.GetTitle())
emptyLen2DHist.Draw("")
effLen2DHist.Draw("same")
canvas.Print(outFileName)

emptyLen3DHist.SetTitle(effLen3DHist.GetTitle())
emptyLen3DHist.Draw("")
effLen3DHist.Draw("same")
canvas.Print(outFileName)

emptyLen3DHist.SetTitle(effLen3DHist.GetTitle())
emptyLen3DHist.Draw("")
effLen3DHist.Draw("same")
effMyPixLen3DHist.Draw("same")
if (regFlag):
  effRegPixLen3DHist.Draw("same")
if (it0Flag):
  effIter0Len3DHist.Draw("same")
if (it1Flag):
  effIter1Len3DHist.Draw("same")
if (it2Flag):
  effIter2Len3DHist.Draw("same")
effItersLen3DLegend.Draw("same")
canvas.Print(outFileName)

canvas.cd(0).SetLogx(0)
emptyHist.SetTitle(effHist.GetTitle())
emptyHist.Draw("")
effItersLegend.Draw("same")
effHist.Draw("same")
effMyPixHist.Draw("same")
if (regFlag):
  effRegPixHist.Draw("same")
if (it0Flag):
  effIter0Hist.Draw("same")
if (it1Flag):
  effIter1Hist.Draw("same")
if (it2Flag):
  effIter2Hist.Draw("same")
line.Draw("same")
canvas.Print(outFileName)

emptyEtaHist.SetTitle(effEtaHist.GetTitle())
emptyEtaHist.Draw()
effItersEtaLegend.Draw("same")
effEtaHist.Draw("same")
effMyPixEtaHist.Draw("same")
if (regFlag):
  effRegPixEtaHist.Draw("same")
if (it0Flag):
  effIter0EtaHist.Draw("same")
if (it1Flag):
  effIter1EtaHist.Draw("same")
if (it2Flag):
  effIter2EtaHist.Draw("same")
canvas.Print(outFileName)

emptyStepHist.SetTitle("Track efficiency at HLTrigger steps")
emptyStepHist.Draw()
effOverItersGraph.Draw("pesame")
canvas.Print(outFileName)

myPixDxyHist.Draw()
regPixDxyHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

myPixDzHist.Draw()
regPixDzHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

myPixDetaHist.Draw()
regPixDetaHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

myPixDphiHist.Draw()
regPixDphiHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

myPixChi2nHist.Draw()
regPixChi2nHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

myPixDxyDzHist.Draw("colz")
canvas.Print(outFileName)

myPixDetaDphiHist.Draw("colz")
canvas.Print(outFileName)

myPixSizeHist.Draw()
canvas.Print(outFileName)

myPixNhitsHist.Draw()
canvas.Print(outFileName) 

regPixDxyDzHist.Draw("colz")
canvas.Print(outFileName)

regPixDetaDphiHist.Draw("colz")
canvas.Print(outFileName)

regPixSizeHist.Draw()
canvas.Print(outFileName)

if (it0Flag):
  ctf0DxyHist.Scale(1/ctf0DxyHist.Integral())
  ctf0DxyHist.Draw("")
if (it1Flag):
  ctf1DxyHist.Scale(1/ctf1DxyHist.Integral())
  ctf1DxyHist.Draw("same")
if (it2Flag):
  ctf2DxyHist.Scale(1/ctf2DxyHist.Integral())
  ctf2DxyHist.Draw("same")
ctfLegend.Draw("same")
canvas.Print(outFileName)

if (it0Flag):
  ctf0DzHist.Scale(1/ctf0DzHist.Integral())
  ctf0DzHist.Draw()
if (it1Flag):
  ctf1DzHist.Scale(1/ctf1DzHist.Integral())
  ctf1DzHist.Draw("same")
if (it2Flag):
  ctf2DzHist.Scale(1/ctf2DzHist.Integral())
  ctf2DzHist.Draw("same")
ctfLegend.Draw("same")
canvas.Print(outFileName)

if (it0Flag):
  ctf0DetaHist.Scale(1/ctf0DetaHist.Integral())
  ctf0DetaHist.Draw()
if (it1Flag):
  ctf1DetaHist.Scale(1/ctf1DetaHist.Integral())
  ctf1DetaHist.Draw("same")
if (it2Flag):
  ctf2DetaHist.Scale(1/ctf2DetaHist.Integral())
  ctf2DetaHist.Draw("same")
ctfLegend.Draw("same")
canvas.Print(outFileName)

if (it0Flag):
  ctf0DphiHist.Scale(1/ctf0DphiHist.Integral())
  ctf0DphiHist.Draw()
if (it1Flag):
  ctf1DphiHist.Scale(1/ctf1DphiHist.Integral())
  ctf1DphiHist.Draw("same")
if (it2Flag):
  ctf2DphiHist.Scale(1/ctf2DphiHist.Integral())
  ctf2DphiHist.Draw("same")
ctfLegend.Draw("same")
canvas.Print(outFileName)

emptySizeHist.SetTitle("Input to track fitting")
emptySizeHist.Draw("")
if (it0Flag):
  ctf0SizeHist.Scale(1/ctf0SizeHist.Integral())
  ctf0SizeHist.Draw("same")
if (it1Flag):
  ctf1SizeHist.Scale(1/ctf1SizeHist.Integral())
  ctf1SizeHist.Draw("same")
if (it2Flag):
  ctf2SizeHist.Scale(1/ctf2SizeHist.Integral())
  ctf2SizeHist.Draw("same")
ctfLegend.Draw("same")
canvas.Print(outFileName)

emptySizeHist.SetTitle("Output from track fit&classifier")
emptySizeHist.Draw("")
if (it0Flag):
  iter0SizeHist.Scale(1/iter0SizeHist.Integral())
  iter0SizeHist.Draw("same")
if (it1Flag):
  iter1SizeHist.Scale(1/iter1SizeHist.Integral())
  iter1SizeHist.Draw("same")
if (it2Flag):
  iter2SizeHist.Scale(1/iter2SizeHist.Integral())
  iter2SizeHist.Draw("same")
iterLegend.Draw("same")
canvas.Print(outFileName)

muonDeltaRHist.Draw()
canvas.Print(outFileName)

trackDeltaRHist.Draw()
canvas.Print(outFileName)

fakePtHist.Draw("")
canvas.Print(outFileName)

fakeEtaHist.Draw("")
canvas.Print(outFileName)

canvas.Print(outFileName+"]")

