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

if len(sys.argv) != 3:
  print "USAGE: %s <input root> <output pdf>"%(sys.argv[0])
  sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]

print "Reading from", inFileName, "and writing to", outFileName

#Read file

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("test/ntupleTree")

#Define hists

binning = np.array([0, 0.5, 1.2, 1.8, 2.5, 3.5, 5.5, 8.0, 12.0, 30.0])
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

legPixHist = ROOT.TH1D("legPixHist", "Legacy tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legPixHist.Sumw2()
legPixEtaHist = ROOT.TH1D("legPixEtaHist", "Legacy tracks;#eta;Counts", 10, -2.5, +2.5)
legPixEtaHist.Sumw2()
legPixLen3DHist = ROOT.TH1D("legPixLen3DHist", "legacy tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
legPixLen3DHist.Sumw2()

legPixSizeHist = ROOT.TH1D("legPixSizeHist", "Size of pixelTracks collection", 20, 0, 20)
legPixDxyHist  = ROOT.TH1D("legPixDxyHist",  "pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
legPixDzHist   = ROOT.TH1D("legPixDzHist",   "pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
legPixDetaHist = ROOT.TH1D("legPixDetaHist", "pixelTracks' eta distance from mu cands", 50, 0, 2)
legPixDphiHist = ROOT.TH1D("legPixDphiHist", "pixelTracks' phi distance from mu cands", 50, 0, 2)
legPixChi2nHist = ROOT.TH1D("legPixChi2nHist", "pixelTracks' reduced chisquare", 10, 0, 10)

legPixDxyDzHist    = ROOT.TH2D("legPixVtxDist",      "pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
legPixDetaDphiHist = ROOT.TH2D("legPixDetaDphiHist", "pixelTracks' angular distance from mu cands", 50, 0, 2, 50, 0, 2)

legIt0Hist = ROOT.TH1D("legIt0Hist", "Iter0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt0Hist.Sumw2()
legIt0EtaHist = ROOT.TH1D("legIt0EtaHist", "Iter0 tracks;#eta;Counts", 10, -2.5, +2.5)
legIt0EtaHist.Sumw2()
legIt0Len3DHist = ROOT.TH1D("legIt0Len3DHist", "Iter0 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
legIt0Len3DHist.Sumw2()

legIt1Hist = ROOT.TH1D("legIt1Hist", "Iter1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt1Hist.Sumw2()
legIt1EtaHist = ROOT.TH1D("legIt1EtaHist", "Iter1 tracks (;#eta;Counts", 10, -2.5, +2.5)
legIt1EtaHist.Sumw2()
legIt1Len3DHist = ROOT.TH1D("legIt1Len3DHist", "Iter1 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
legIt1Len3DHist.Sumw2()

legIt2Hist = ROOT.TH1D("legIt2Hist", "Iter2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt2Hist.Sumw2()
legIt2EtaHist = ROOT.TH1D("legIt2EtaHist", "Iter2 tracks;#eta;Counts", 10, -2.5, +2.5)
legIt2EtaHist.Sumw2()
legIt2Len3DHist = ROOT.TH1D("legIt2Len3DHist", "Iter2 tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
legIt2Len3DHist.Sumw2()

leg3ReHist = ROOT.TH1D("leg3ReHist", "3Reco tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
leg3ReHist.Sumw2()
leg3ReEtaHist = ROOT.TH1D("leg3ReEtaHist", "3Reco tracks;#eta;Counts", 10, -2.5, +2.5)
leg3ReEtaHist.Sumw2()
leg3ReLen3DHist = ROOT.TH1D("leg3ReLen3DHist", "3Reco tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
leg3ReLen3DHist.Sumw2()

leg2ReHist = ROOT.TH1D("leg2ReHist", "2Reco tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
leg2ReHist.Sumw2()
leg2ReEtaHist = ROOT.TH1D("leg2ReEtaHist", "2Reco tracks;#eta;Counts", 10, -2.5, +2.5)
leg2ReEtaHist.Sumw2()
leg2ReLen3DHist = ROOT.TH1D("leg2ReLen3DHist", "2Reco tracks;Bs 3D flight length;Counts", len(length3DBinning)-1, length3DBinning)
leg2ReLen3DHist.Sumw2()

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

  hltJpsiTrkMatched = getattr(tree, "hltMatch_Track")
  hltJpsiMatched = getattr(tree, "muonsMatched")
  hltTrack1Matched = getattr(tree, "hltMatch_Track1")
  hltTrack2Matched = getattr(tree, "hltMatch_Track2")

  hltTrigPtVec  = getattr(tree, "hltTriggerTrack_pt")
  hltTrigEtaVec = getattr(tree, "hltTriggerTrack_eta")
  hltFakePtVec  = getattr(tree, "hltTriggerFake_pt")
  hltFakeEtaVec = getattr(tree, "hltTriggerFake_eta")

  legPixTrack1Matched = getattr(tree, "legPixMatch_Track1")
  legPixTrack2Matched = getattr(tree, "legPixMatch_Track2")
  legIt0Track1Matched = getattr(tree, "legI0Match_Track1")
  legIt0Track2Matched = getattr(tree, "legI0Match_Track2")
  legIt1Track1Matched = getattr(tree, "legI1Match_Track1")
  legIt1Track2Matched = getattr(tree, "legI1Match_Track2")
  legIt2Track1Matched = getattr(tree, "legI2Match_Track1")
  legIt2Track2Matched = getattr(tree, "legI2Match_Track2")
  leg3ReTrack1Matched = getattr(tree, "leg3RMatch_Track1")
  leg3ReTrack2Matched = getattr(tree, "leg3RMatch_Track2")
  leg2ReTrack1Matched = getattr(tree, "leg2RMatch_Track1")
  leg2ReTrack2Matched = getattr(tree, "leg2RMatch_Track2")

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

  legPixTrack1dxy   = getattr(tree, "legPix_Track1_dxy")
  legPixTrack1dz    = getattr(tree, "legPix_Track1_dz")
  legPixTrack1dEta1 = getattr(tree, "legPix_Track1_dEtaMu1")
  legPixTrack1dPhi1 = getattr(tree, "legPix_Track1_dPhiMu1")
  legPixTrack1dEta2 = getattr(tree, "legPix_Track1_dEtaMu2")
  legPixTrack1dPhi2 = getattr(tree, "legPix_Track1_dPhiMu2")
  legPixTrack1chi2n = getattr(tree, "legPix_Track1_normChi2")

  legPixTrack2dxy   = getattr(tree, "legPix_Track2_dxy")
  legPixTrack2dz    = getattr(tree, "legPix_Track2_dz")
  legPixTrack2dEta1 = getattr(tree, "legPix_Track2_dEtaMu1")
  legPixTrack2dPhi1 = getattr(tree, "legPix_Track2_dPhiMu1")
  legPixTrack2dEta2 = getattr(tree, "legPix_Track2_dEtaMu2")
  legPixTrack2dPhi2 = getattr(tree, "legPix_Track2_dPhiMu2")
  legPixTrack2chi2n = getattr(tree, "legPix_Track2_normChi2")

  legPixSize = getattr(tree, "legPix_collSize")

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

  #Legacy pixelTracks

  if (hltJpsiMatched):
    legPixSizeHist.Fill(legPixSize)
    if (legPixTrack1Matched):
      legPixHist.Fill(genTrack1Pt)
      legPixEtaHist.Fill(genTrack1eta)
      legPixLen3DHist.Fill(genTrack1Len3D)
      legPixDxyHist.Fill(legPixTrack1dxy)
      legPixDzHist.Fill(legPixTrack1dz)
      legPixDetaHist.Fill(legPixTrack1dEta1)
      legPixDetaHist.Fill(legPixTrack1dEta2)
      legPixDphiHist.Fill(legPixTrack1dPhi1)
      legPixDphiHist.Fill(legPixTrack1dPhi2)
      legPixChi2nHist.Fill(legPixTrack1chi2n)
      legPixDxyDzHist.Fill(legPixTrack1dxy, legPixTrack1dz)
      legPixDetaDphiHist.Fill(legPixTrack1dEta1, legPixTrack1dPhi1)
      legPixDetaDphiHist.Fill(legPixTrack1dEta2, legPixTrack1dPhi2)
    if (legPixTrack2Matched):
      legPixHist.Fill(genTrack2Pt)
      legPixEtaHist.Fill(genTrack2eta)
      legPixLen3DHist.Fill(genTrack2Len3D)
      legPixDxyHist.Fill(legPixTrack2dxy)
      legPixDzHist.Fill(legPixTrack2dz)
      legPixDetaHist.Fill(legPixTrack2dEta1)
      legPixDetaHist.Fill(legPixTrack2dEta2)
      legPixDphiHist.Fill(legPixTrack2dPhi1)
      legPixDphiHist.Fill(legPixTrack2dPhi2)
      legPixChi2nHist.Fill(legPixTrack2chi2n)
      legPixDxyDzHist.Fill(legPixTrack2dxy, legPixTrack2dz)
      legPixDetaDphiHist.Fill(legPixTrack2dEta1, legPixTrack2dPhi1)
      legPixDetaDphiHist.Fill(legPixTrack2dEta2, legPixTrack2dPhi2)

    if (legIt0Track1Matched):
      legIt0Hist.Fill(genTrack1Pt)
      legIt0EtaHist.Fill(genTrack1eta)
      legIt0Len3DHist.Fill(genTrack1Len3D)
    if (legIt0Track2Matched):
      legIt0Hist.Fill(genTrack2Pt)
      legIt0EtaHist.Fill(genTrack2eta)
      legIt0Len3DHist.Fill(genTrack2Len3D)

    if (legIt1Track1Matched):
      legIt1Hist.Fill(genTrack1Pt)
      legIt1EtaHist.Fill(genTrack1eta)
      legIt1Len3DHist.Fill(genTrack1Len3D)
    if (legIt1Track2Matched):
      legIt1Hist.Fill(genTrack2Pt)
      legIt1EtaHist.Fill(genTrack2eta)
      legIt1Len3DHist.Fill(genTrack2Len3D)

    if (legIt2Track1Matched):
      legIt2Hist.Fill(genTrack1Pt)
      legIt2EtaHist.Fill(genTrack1eta)
      legIt2Len3DHist.Fill(genTrack1Len3D)
    if (legIt2Track2Matched):
      legIt2Hist.Fill(genTrack2Pt)
      legIt2EtaHist.Fill(genTrack2eta)
      legIt2Len3DHist.Fill(genTrack2Len3D)

    if (leg3ReTrack1Matched):
      leg3ReHist.Fill(genTrack1Pt)
      leg3ReEtaHist.Fill(genTrack1eta)
      leg3ReLen3DHist.Fill(genTrack1Len3D)
    if (leg3ReTrack2Matched):
      leg3ReHist.Fill(genTrack2Pt)
      leg3ReEtaHist.Fill(genTrack2eta)
      leg3ReLen3DHist.Fill(genTrack2Len3D)

    if (leg2ReTrack1Matched):
      leg2ReHist.Fill(genTrack1Pt)
      leg2ReEtaHist.Fill(genTrack1eta)
      leg2ReLen3DHist.Fill(genTrack1Len3D)
    if (leg2ReTrack2Matched):
      leg2ReHist.Fill(genTrack2Pt)
      leg2ReEtaHist.Fill(genTrack2eta)
      leg2ReLen3DHist.Fill(genTrack2Len3D) 

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

#Legacy pixel tracks

effLegPixHist = ROOT.TEfficiency(legPixHist, denHist)
effLegPixHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegPixHist.SetConfidenceLevel(0.68)
effLegPixHist.SetNameTitle("effLegPixHist", "legacy Tracks efficiency")

effLegPixEtaHist = ROOT.TEfficiency(legPixEtaHist, denEtaHist)
effLegPixEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegPixEtaHist.SetConfidenceLevel(0.68)
effLegPixEtaHist.SetNameTitle("legPixEtaHist", "legacy Tracks efficiency")

effLegPixLen3DHist = ROOT.TEfficiency(legPixLen3DHist, denLen3DHist)
effLegPixLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegPixLen3DHist.SetConfidenceLevel(0.68)
effLegPixLen3DHist.SetNameTitle("effLegPixLen3DHist", "legacy Tracks efficiency")

#iter0

effLegIt0Hist = ROOT.TEfficiency(legIt0Hist, denHist)
effLegIt0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt0Hist.SetConfidenceLevel(0.68)
effLegIt0Hist.SetNameTitle("effLegIt0Hist", "legacy Iter0 Track efficiency")

effLegIt0EtaHist = ROOT.TEfficiency(legIt0EtaHist, denEtaHist)
effLegIt0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt0EtaHist.SetConfidenceLevel(0.68)
effLegIt0EtaHist.SetNameTitle("effLegIt0EtaHist", "legacy Iter0 Track efficiency")

effLegIt0Len3DHist = ROOT.TEfficiency(legIt0Len3DHist, denLen3DHist)
effLegIt0Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt0Len3DHist.SetConfidenceLevel(0.68)
effLegIt0Len3DHist.SetNameTitle("effLegIt0Len3DHist", "legacy Iter0 Tracks efficiency")

#iter1

effLegIt1Hist = ROOT.TEfficiency(legIt1Hist, denHist)
effLegIt1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt1Hist.SetConfidenceLevel(0.68)
effLegIt1Hist.SetNameTitle("effLegIt1Hist", "legacy Iter1 Track efficiency")

effLegIt1EtaHist = ROOT.TEfficiency(legIt1EtaHist, denEtaHist)
effLegIt1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt1EtaHist.SetConfidenceLevel(0.68)
effLegIt1EtaHist.SetNameTitle("effLegIt1EtaHist", "legacy Iter1 Track efficiency")

effLegIt1Len3DHist = ROOT.TEfficiency(legIt1Len3DHist, denLen3DHist)
effLegIt1Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt1Len3DHist.SetConfidenceLevel(0.68)
effLegIt1Len3DHist.SetNameTitle("effLegIt1Len3DHist", "legacy Iter1 Tracks efficiency")

#iter2

effLegIt2Hist = ROOT.TEfficiency(legIt2Hist, denHist)
effLegIt2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt2Hist.SetConfidenceLevel(0.68)
effLegIt2Hist.SetNameTitle("effLegIt2Hist", "legacy Iter2 Track efficiency")

effLegIt2EtaHist = ROOT.TEfficiency(legIt2EtaHist, denEtaHist)
effLegIt2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt2EtaHist.SetConfidenceLevel(0.68)
effLegIt2EtaHist.SetNameTitle("effLegIt2EtaHist", "legacy Iter2 Track efficiency")

effLegIt2Len3DHist = ROOT.TEfficiency(legIt2Len3DHist, denLen3DHist)
effLegIt2Len3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt2Len3DHist.SetConfidenceLevel(0.68)
effLegIt2Len3DHist.SetNameTitle("effLegIt2Len3DHist", "legacy Iter2 Tracks efficiency")

#triplet recovery

effLeg3ReHist = ROOT.TEfficiency(leg3ReHist, denHist)
effLeg3ReHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg3ReHist.SetConfidenceLevel(0.68)
effLeg3ReHist.SetNameTitle("effLeg3ReHist", "legacy 3-reco Track efficiency")

effLeg3ReEtaHist = ROOT.TEfficiency(leg3ReEtaHist, denEtaHist)
effLeg3ReEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg3ReEtaHist.SetConfidenceLevel(0.68)
effLeg3ReEtaHist.SetNameTitle("effLeg3ReEtaHist", "legacy 3-reco Track efficiency")

effLeg3ReLen3DHist = ROOT.TEfficiency(leg3ReLen3DHist, denLen3DHist)
effLeg3ReLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg3ReLen3DHist.SetConfidenceLevel(0.68)
effLeg3ReLen3DHist.SetNameTitle("effLeg3ReLen3DHist", "legacy 3-reco Tracks efficiency")

#doublet recovery

effLeg2ReHist = ROOT.TEfficiency(leg2ReHist, denHist)
effLeg2ReHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg2ReHist.SetConfidenceLevel(0.68)
effLeg2ReHist.SetNameTitle("effLeg2ReHist", "legacy 2-reco Track efficiency")

effLeg2ReEtaHist = ROOT.TEfficiency(leg2ReEtaHist, denEtaHist)
effLeg2ReEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg2ReEtaHist.SetConfidenceLevel(0.68)
effLeg2ReEtaHist.SetNameTitle("effLeg2ReEtaHist", "legacy 2-reco Track efficiency")

effLeg2ReLen3DHist = ROOT.TEfficiency(leg2ReLen3DHist, denLen3DHist)
effLeg2ReLen3DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg2ReLen3DHist.SetConfidenceLevel(0.68)
effLeg2ReLen3DHist.SetNameTitle("effLeg2ReLen3DHist", "legacy 2-reco Tracks efficiency")

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

xIter = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
xIterErr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

yEff = np.array([averageOverThreshold1p2(legPixHist, denHist), averageOverThreshold1p2(legIt0Hist, denHist), averageOverThreshold1p2(legIt1Hist, denHist), averageOverThreshold1p2(legIt2Hist, denHist), averageOverThreshold1p2(leg3ReHist, denHist), averageOverThreshold1p2(leg2ReHist, denHist), averageOverThreshold1p2(numHist, denHist)])

yEffError =  np.array([averageOverThreshold1p2Error(legPixHist, denHist), averageOverThreshold1p2Error(legIt0Hist, denHist), averageOverThreshold1p2Error(legIt1Hist, denHist), averageOverThreshold1p2Error(legIt2Hist, denHist), averageOverThreshold1p2Error(leg3ReHist, denHist), averageOverThreshold1p2Error(leg2ReHist, denHist), averageOverThreshold1p2Error(numHist, denHist)])

effOverItersGraph = ROOT.TGraphErrors(7, xIter, yEff, xIterErr, yEffError)
effOverItersGraph.SetNameTitle("effOverIters", "track efficiency at trigger steps")
effOverItersGraph.GetXaxis().SetTitle("Steps")
effOverItersGraph.GetYaxis().SetTitle("Efficiency")

#Set directories

effHist.SetDirectory(0)
effEtaHist.SetDirectory(0)
effPhiHist.SetDirectory(0)
effLen2DHist.SetDirectory(0)
effLen3DHist.SetDirectory(0)

effLegPixHist.SetDirectory(0)
effLegPixEtaHist.SetDirectory(0)
effLegPixLen3DHist.SetDirectory(0)

legPixSizeHist.SetDirectory(0)
legPixDxyHist.SetDirectory(0)
legPixDzHist.SetDirectory(0)
legPixDetaHist.SetDirectory(0)
legPixDphiHist.SetDirectory(0)
legPixChi2nHist.SetDirectory(0)
legPixDxyDzHist.SetDirectory(0)
legPixDetaDphiHist.SetDirectory(0)

effLegIt0Hist.SetDirectory(0)
effLegIt0EtaHist.SetDirectory(0)
effLegIt0Len3DHist.SetDirectory(0)
effLegIt1Hist.SetDirectory(0)
effLegIt1EtaHist.SetDirectory(0)
effLegIt1Len3DHist.SetDirectory(0)
effLegIt2Hist.SetDirectory(0)
effLegIt2EtaHist.SetDirectory(0)
effLegIt2Len3DHist.SetDirectory(0)
effLeg3ReHist.SetDirectory(0)
effLeg3ReEtaHist.SetDirectory(0)
effLeg3ReLen3DHist.SetDirectory(0)
effLeg2ReHist.SetDirectory(0)
effLeg2ReEtaHist.SetDirectory(0)
effLeg2ReLen3DHist.SetDirectory(0)

muonDeltaRHist.SetDirectory(0)
trackDeltaRHist.SetDirectory(0)

numLen2DHist.SetDirectory(0)
numLen3DHist.SetDirectory(0)

fakePtHist.SetDirectory(0)
fakeEtaHist.SetDirectory(0)

inFile.Close()

#set graphics

emptyHist = ROOT.TH1D("emptyHist", "title;p_{T} [GeV/c];Efficiency", 10, 0, 30)
emptyEtaHist = ROOT.TH1D("emptyEtaHist", "title;#eta;Efficiency", 10, -2.5, 2.5)
emptyPhiHist = ROOT.TH1D("emptyPhiHist", "title;#phi;Efficiency", 8, -3.2, +3.2) 
emptyLen2DHist = ROOT.TH1D("emptyLen2DHist", "title;Bs 2D flight length [cm];Efficiency", len(lengthBinning)-1, lengthBinning)
emptyLen3DHist = ROOT.TH1D("emptyLen3DHist", "title;Bs 3D flight length [cm];Efficiency", len(length3DBinning)-1, length3DBinning)
emptyStepHist = ROOT.TH1D("emptyStepHist", "title;;Efficiency", 7, 0, 7)

emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyEtaHist.GetYaxis().SetRangeUser(0, 1.05)
emptyPhiHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen2DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen3DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyStepHist.GetYaxis().SetRangeUser(0., 1.05)

emptyHist.SetStats(ROOT.kFALSE)
emptyEtaHist.SetStats(ROOT.kFALSE)
emptyPhiHist.SetStats(ROOT.kFALSE)
emptyLen2DHist.SetStats(ROOT.kFALSE)
emptyLen3DHist.SetStats(ROOT.kFALSE)
emptyStepHist.SetStats(ROOT.kFALSE)

idx = 0
binLabels = ["Legacy pixeltracks", "Iter0", "Iter1", "Iter2", "3letRecJpsiReg", "2letRecJpsiReg", "After cuts"]
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

effLegPixHist.SetLineColor(7)
effLegPixEtaHist.SetLineColor(7)
effLegPixLen3DHist.SetLineColor(7)
effLegIt0Hist.SetLineColor(6)
effLegIt0EtaHist.SetLineColor(6)
effLegIt0Len3DHist.SetLineColor(6)
effLegIt1Hist.SetLineColor(2)
effLegIt1EtaHist.SetLineColor(2)
effLegIt1Len3DHist.SetLineColor(2)
effLegIt2Hist.SetLineColor(3)
effLegIt2EtaHist.SetLineColor(3)
effLegIt2Len3DHist.SetLineColor(3)
effLeg3ReHist.SetLineColor(4)
effLeg3ReEtaHist.SetLineColor(4)
effLeg3ReLen3DHist.SetLineColor(4)
effLeg2ReHist.SetLineColor(5)
effLeg2ReEtaHist.SetLineColor(5)
effLeg2ReLen3DHist.SetLineColor(5)

fakePtHist.SetLineColor(1)
fakeEtaHist.SetLineColor(1)

effHist.SetLineWidth(2)
effEtaHist.SetLineWidth(2)
effPhiHist.SetLineWidth(2)
effLen2DHist.SetLineWidth(2)
effLen3DHist.SetLineWidth(2)

effLegPixHist.SetLineWidth(2)
effLegPixEtaHist.SetLineWidth(2)
effLegPixLen3DHist.SetLineWidth(2)
effLegIt0Hist.SetLineWidth(2)
effLegIt0EtaHist.SetLineWidth(2)
effLegIt0Len3DHist.SetLineWidth(2)
effLegIt1Hist.SetLineWidth(2)
effLegIt1EtaHist.SetLineWidth(2)
effLegIt1Len3DHist.SetLineWidth(2)
effLegIt2Hist.SetLineWidth(2)
effLegIt2EtaHist.SetLineWidth(2)
effLegIt2Len3DHist.SetLineWidth(2)
effLeg3ReHist.SetLineWidth(2)
effLeg3ReEtaHist.SetLineWidth(2)
effLeg3ReLen3DHist.SetLineWidth(2)
effLeg2ReHist.SetLineWidth(2)
effLeg2ReEtaHist.SetLineWidth(2)
effLeg2ReLen3DHist.SetLineWidth(2)

fakePtHist.SetLineWidth(2)
fakeEtaHist.SetLineWidth(2)

effOverItersGraph.SetLineColor(1)
effOverItersGraph.SetMarkerStyle(2)

legPixSizeHist.GetYaxis().SetTitle("Counts")
legPixDxyHist.GetYaxis().SetTitle("Counts")
legPixDzHist.GetYaxis().SetTitle("Counts")
legPixDetaHist.GetYaxis().SetTitle("Counts")
legPixDphiHist.GetYaxis().SetTitle("Counts")
legPixChi2nHist.GetYaxis().SetTitle("Counts")
legPixDxyDzHist.GetYaxis().SetTitle("#Delta z [cm]")
legPixDetaDphiHist.GetYaxis().SetTitle("#Delta #phi")

legPixSizeHist.GetXaxis().SetTitle("Track collection's size")
legPixDxyHist.GetXaxis().SetTitle("#Delta xy [cm]")
legPixDzHist.GetXaxis().SetTitle("#Delta z [cm]")
legPixDetaHist.GetXaxis().SetTitle("#Delta #eta")
legPixDphiHist.GetXaxis().SetTitle("#Delta #phi")
legPixChi2nHist.GetXaxis().SetTitle("#Chi^{2}")
legPixDxyDzHist.GetXaxis().SetTitle("#Delta xy [cm]")
legPixDetaDphiHist.GetXaxis().SetTitle("#Delta #eta")

fakePtHist.GetYaxis().SetTitle("Fake rate")
fakeEtaHist.GetYaxis().SetTitle("Fake rate")

effItersLegend = ROOT.TLegend(0.7, 0.12, 0.88, 0.35)
effItersLegend.SetTextSize(0.03)
effItersLegend.AddEntry(effLegPixHist, "pixelTracks")
effItersLegend.AddEntry(effLegIt0Hist, "iter0")
effItersLegend.AddEntry(effLegIt1Hist, "iter1")
effItersLegend.AddEntry(effLegIt2Hist, "iter2")
effItersLegend.AddEntry(effLeg3ReHist, "3let Reco")
effItersLegend.AddEntry(effLeg2ReHist, "2let Reco")
effItersLegend.AddEntry(effHist, "afterFilter")

effItersEtaLegend = ROOT.TLegend(0.7, 0.12, 0.88, 0.35)
effItersEtaLegend.SetTextSize(0.03)
effItersEtaLegend.AddEntry(effLegPixEtaHist, "pixelTracks")
effItersEtaLegend.AddEntry(effLegIt0EtaHist, "iter0")
effItersEtaLegend.AddEntry(effLegIt1EtaHist, "iter1")
effItersEtaLegend.AddEntry(effLegIt2EtaHist, "iter2")
effItersEtaLegend.AddEntry(effLeg3ReEtaHist, "3let Reco")
effItersEtaLegend.AddEntry(effLeg2ReEtaHist, "2let Reco")
effItersEtaLegend.AddEntry(effEtaHist, "afterFilter")

effItersLen3DLegend = ROOT.TLegend(0.25, 0.15, 0.7, 0.28)
effItersLen3DLegend.SetTextSize(0.03)
effItersLen3DLegend.SetNColumns(2)
effItersLen3DLegend.AddEntry(effLegPixLen3DHist, "pixelTracks")
effItersLen3DLegend.AddEntry(effLegIt0Len3DHist, "iter0")
effItersLen3DLegend.AddEntry(effLegIt1Len3DHist, "iter1")
effItersLen3DLegend.AddEntry(effLegIt2Len3DHist, "iter2")
effItersLen3DLegend.AddEntry(effLeg3ReLen3DHist, "3let Reco")
effItersLen3DLegend.AddEntry(effLeg2ReLen3DHist, "2let Reco")
effItersLen3DLegend.AddEntry(effLen3DHist, "afterFilter")

histLegend = ROOT.TLegend(0.65, 0.68, 0.89, 0.73)
histLegend.SetTextSize(0.03)
histLegend.AddEntry(legPixSizeHist, "legacy pixelTracks")

#Define canvas and print

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
effItersLen3DLegend.Draw("same")
effLen3DHist.Draw("same")
effLegPixLen3DHist.Draw("same")
effLegIt0Len3DHist.Draw("same")
effLegIt1Len3DHist.Draw("same")
effLegIt2Len3DHist.Draw("same")
effLeg3ReLen3DHist.Draw("same")
effLeg2ReLen3DHist.Draw("same")
canvas.Print(outFileName)

canvas.cd(0).SetLogx(0)
emptyHist.SetTitle(effHist.GetTitle())
emptyHist.Draw("")
effItersLegend.Draw("same")
effHist.Draw("same")
effLegPixHist.Draw("same")
effLegIt0Hist.Draw("same")
effLegIt1Hist.Draw("same")
effLegIt2Hist.Draw("same")
effLeg3ReHist.Draw("same")
effLeg2ReHist.Draw("same")
line.Draw("same")
canvas.Print(outFileName)

emptyEtaHist.SetTitle(effEtaHist.GetTitle())
emptyEtaHist.Draw()
effItersEtaLegend.Draw("same")
effEtaHist.Draw("same")
effLegPixEtaHist.Draw("same")
effLegIt0EtaHist.Draw("same")
effLegIt1EtaHist.Draw("same")
effLegIt2EtaHist.Draw("same")
effLeg3ReEtaHist.Draw("same")
effLeg2ReEtaHist.Draw("same")
canvas.Print(outFileName)

emptyStepHist.SetTitle("Track efficiency at HLTrigger steps")
emptyStepHist.Draw()
effOverItersGraph.Draw("pesame")
canvas.Print(outFileName)

legPixSizeHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixDxyHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixDzHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixDetaHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixDphiHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixChi2nHist.Draw()
histLegend.Draw("same")
canvas.Print(outFileName)

legPixDxyDzHist.Draw("colz")
canvas.Print(outFileName)

legPixDetaDphiHist.Draw("colz")
canvas.Print(outFileName)

muonDeltaRHist.Draw()
canvas.Print(outFileName)

trackDeltaRHist.Draw()
canvas.Print(outFileName)

fakePtHist.Draw()
canvas.Print(outFileName)

fakeEtaHist.Draw()
canvas.Print(outFileName)

canvas.Print(outFileName+"]")
