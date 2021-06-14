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
  print "USAGE: %s <input file1> <output file>"%(sys.argv[0])
  sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]

print "Reading from", inFileName, "and writing to", outFileName

#Read file

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("test/ntupleTree")

#Define hists

binning = np.array([0, 0.6, 1.2, 1.8, 2.5, 3.5, 5.5, 8.0, 12.0, 30.0])
lengthBinning = np.array([0, 3e-3, 9e-3, 3e-2, 9e-2, 3e-1, 9e-1, 3, 9])

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

numLen3DHist = ROOT.TH1D("numLen3DHist", "HLT matched tracks;Bs 3D flight length;Counts", len(lengthBinning)-1, lengthBinning)
numLen3DHist.Sumw2()
denLen3DHist = ROOT.TH1D("denLen3DHist", "Tracks;Bs 3D flight length;Counts", len(lengthBinning)-1, lengthBinning)
denLen3DHist.Sumw2()
denLen3DPtCutHist = ROOT.TH1D("denLen3DPtCutHist", "Tracks pT>1.2;Bs 3D flight length;Counts", len(lengthBinning)-1, lengthBinning)
denLen3DPtCutHist.Sumw2()

trackDeltaRHist = ROOT.TH1D("trackDeltaR", "#DeltaR hlt vs gen;Track #DeltaR;Counts", 20, 0, 0.01)
muonDeltaRHist = ROOT.TH1D("muonDeltaR", "#DeltaR hlt vs gen;Muon #DeltaR;Counts", 20, 0, 0.01)

legPixHist = ROOT.TH1D("legPixHist", "Legacy tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legPixHist.Sumw2()
legPixEtaHist = ROOT.TH1D("legPixEtaHist", "Legacy tracks;#eta;Counts", 10, -2.5, +2.5)
legPixEtaHist.Sumw2()

legPixSizeHist = ROOT.TH1D("legPixSizeHist", "Size of pixelTracks collection", 30, 0, 3000)
legPixDxyHist  = ROOT.TH1D("legPixDxyHist",  "pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
legPixDzHist   = ROOT.TH1D("legPixDzHist",   "pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
legPixDetaHist = ROOT.TH1D("legPixDetaHist", "pixelTracks' eta distance from mu cands", 50, 0, 2)
legPixDphiHist = ROOT.TH1D("legPixDphiHist", "pixelTracks' phi distance from mu cands", 50, 0, 2)
legPixChi2nHist = ROOT.TH1D("legPixChi2nHist", "pixelTrakcs' reduced chisquare", len(lengthBinning), lengthBinning)

legPixDxyDzHist    = ROOT.TH2D("legPixVtxDist",      "pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
legPixDetaDphiHist = ROOT.TH2D("legPixDetaDphiHist", "pixelTracks' angular distance from mu cands", 50, 0, 2, 50, 0, 2)

legIt0Hist = ROOT.TH1D("legIt0Hist", "Iter0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt0Hist.Sumw2()
legIt0EtaHist = ROOT.TH1D("legIt0EtaHist", "Iter0 tracks;#eta;Counts", 10, -2.5, +2.5)
legIt0EtaHist.Sumw2()

legIt1Hist = ROOT.TH1D("legIt1Hist", "Iter1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt1Hist.Sumw2()
legIt1EtaHist = ROOT.TH1D("legIt1EtaHist", "Iter1 tracks (;#eta;Counts", 10, -2.5, +2.5)
legIt1EtaHist.Sumw2()

legIt2Hist = ROOT.TH1D("legIt2Hist", "Iter2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
legIt2Hist.Sumw2()
legIt2EtaHist = ROOT.TH1D("legIt2EtaHist", "Iter2 tracks;#eta;Counts", 10, -2.5, +2.5)
legIt2EtaHist.Sumw2()

leg3ReHist = ROOT.TH1D("leg3ReHist", "3Reco tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
leg3ReHist.Sumw2()
leg3ReEtaHist = ROOT.TH1D("leg3ReEtaHist", "3Reco tracks;#eta;Counts", 10, -2.5, +2.5)
leg3ReEtaHist.Sumw2()

leg2ReHist = ROOT.TH1D("leg2ReHist", "2Reco tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
leg2ReHist.Sumw2()
leg2ReEtaHist = ROOT.TH1D("leg2ReEtaHist", "2Reco tracks;#eta;Counts", 10, -2.5, +2.5)
leg2ReEtaHist.Sumw2()

fakePtHist = ROOT.TH1D("fakePtHist", "Unmatched trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
fakePtHist.Sumw2()
fakeEtaHist = ROOT.TH1D("fakeEtaHist", "Unmatched trigger tracks;#eta;Counts", 10, -2.5, 2.5)
fakeEtaHist.Sumw2()

trigPtHist = ROOT.TH1D("trigPtHist", "All trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
trigPtHist.Sumw2()
trigEtaHist = ROOT.TH1D("trigEtaHist", "All trigger tracks;#eta;Counts", 10, -2.5, 2.5)
trigEtaHist.Sumw2()

#Fill hists

sumVal = 0
counts = 0

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

  if (hltJpsiMatched):
    if (legPixTrack1Matched):
      legPixHist.Fill(genTrack1Pt)
      legPixEtaHist.Fill(genTrack1eta)
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
    if (legIt0Track2Matched):
      legIt0Hist.Fill(genTrack2Pt)
      legIt0EtaHist.Fill(genTrack2eta)

    if (legIt1Track1Matched):
      legIt1Hist.Fill(genTrack1Pt)
      legIt1EtaHist.Fill(genTrack1eta)
    if (legIt1Track2Matched):
      legIt1Hist.Fill(genTrack2Pt)
      legIt1EtaHist.Fill(genTrack2eta)

    if (legIt2Track1Matched):
      legIt2Hist.Fill(genTrack1Pt)
      legIt2EtaHist.Fill(genTrack1eta)
    if (legIt2Track2Matched):
      legIt2Hist.Fill(genTrack2Pt)
      legIt2EtaHist.Fill(genTrack2eta)

    if (leg3ReTrack1Matched):
      leg3ReHist.Fill(genTrack1Pt)
      leg3ReEtaHist.Fill(genTrack1eta)
    if (leg3ReTrack2Matched):
      leg3ReHist.Fill(genTrack2Pt)
      leg3ReEtaHist.Fill(genTrack2eta)

    if (leg2ReTrack1Matched):
      leg2ReHist.Fill(genTrack1Pt)
      leg2ReEtaHist.Fill(genTrack1eta)
    if (leg2ReTrack2Matched):
      leg2ReHist.Fill(genTrack2Pt)
      leg2ReEtaHist.Fill(genTrack2eta)

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
effEtaHist.SetNameTitle("effEtaHist", "Track efficiency")

effPhiHist = ROOT.TEfficiency(numPhiHist, denPhiPtCutHist)
effPhiHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effPhiHist.SetConfidenceLevel(0.68)
effPhiHist.SetNameTitle("effPhiHist", "Track efficiency")

effLen2DHist = ROOT.TEfficiency(numLen2DHist, denLen2DPtCutHist)
effLen2DHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLen2DHist.SetConfidenceLevel(0.68)
effLen2DHist.SetNameTitle("effLen2DHist", "Track efficiency")

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
effLegPixEtaHist.SetNameTitle("effLegPixEtaHist", "legacy Tracks efficiency")

#iter0

effLegIt0Hist = ROOT.TEfficiency(legIt0Hist, denHist)
effLegIt0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt0Hist.SetConfidenceLevel(0.68)
effLegIt0Hist.SetNameTitle("effLegIt0Hist", "legacy Iter0 Track efficiency")

effLegIt0EtaHist = ROOT.TEfficiency(legIt0EtaHist, denEtaHist)
effLegIt0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt0EtaHist.SetConfidenceLevel(0.68)
effLegIt0EtaHist.SetNameTitle("effLegIt0EtaHist", "legacy Iter0 Track efficiency")

#iter1

effLegIt1Hist = ROOT.TEfficiency(legIt1Hist, denHist)
effLegIt1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt1Hist.SetConfidenceLevel(0.68)
effLegIt1Hist.SetNameTitle("effLegIt1Hist", "legacy Iter1 Track efficiency")

effLegIt1EtaHist = ROOT.TEfficiency(legIt1EtaHist, denEtaHist)
effLegIt1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt1EtaHist.SetConfidenceLevel(0.68)
effLegIt1EtaHist.SetNameTitle("effLegIt1EtaHist", "legacy Iter1 Track efficiency")

#iter2

effLegIt2Hist = ROOT.TEfficiency(legIt2Hist, denHist)
effLegIt2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt2Hist.SetConfidenceLevel(0.68)
effLegIt2Hist.SetNameTitle("effLegIt2Hist", "legacy Iter2 Track efficiency")

effLegIt2EtaHist = ROOT.TEfficiency(legIt2EtaHist, denEtaHist)
effLegIt2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLegIt2EtaHist.SetConfidenceLevel(0.68)
effLegIt2EtaHist.SetNameTitle("effLegIt2EtaHist", "legacy Iter2 Track efficiency")

#triplet recovery

effLeg3ReHist = ROOT.TEfficiency(leg3ReHist, denHist)
effLeg3ReHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg3ReHist.SetConfidenceLevel(0.68)
effLeg3ReHist.SetNameTitle("effLeg3ReHist", "legacy 3-reco Track efficiency")

effLeg3ReEtaHist = ROOT.TEfficiency(leg3ReEtaHist, denEtaHist)
effLeg3ReEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg3ReEtaHist.SetConfidenceLevel(0.68)
effLeg3ReEtaHist.SetNameTitle("effLeg3ReEtaHist", "legacy 3-reco Track efficiency")

#doublet recovery

effLeg2ReHist = ROOT.TEfficiency(leg2ReHist, denHist)
effLeg2ReHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg2ReHist.SetConfidenceLevel(0.68)
effLeg2ReHist.SetNameTitle("effLeg2ReHist", "legacy 2-reco Track efficiency")

effLeg2ReEtaHist = ROOT.TEfficiency(leg2ReEtaHist, denEtaHist)
effLeg2ReEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effLeg2ReEtaHist.SetConfidenceLevel(0.68)
effLeg2ReEtaHist.SetNameTitle("effLeg2ReEtaHist", "legacy 2-reco Track efficiency")

#Fake rate

fakeRatePtHist = ROOT.TEfficiency(fakePtHist, trigPtHist)
fakeRatePtHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
fakeRatePtHist.SetConfidenceLevel(0.68)
fakeRatePtHist.SetNameTitle("fakeRatePtHist", "Fake rate distribution")

fakeRateEtaHist = ROOT.TEfficiency(fakeEtaHist, trigEtaHist)
fakeRateEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
fakeRateEtaHist.SetConfidenceLevel(0.68)
fakeRateEtaHist.SetNameTitle("fakeRateEtaHist", "Fake rate distribution")

num = fakeEtaHist.Integral()
den = trigEtaHist.Integral()

if (den == 0):
  ratio = 0
else:
  ratio = num/den

ratio = averageOverThreshold1p2(fakePtHist, trigPtHist)
errRatio = averageOverThreshold1p2Error(fakePtHist, trigPtHist)
print "fake rate:", ratio, "+/-", errRatio

############################

denHist.SetYTitle("Counts")
denEtaHist.SetYTitle("Counts")
denLen2DHist.SetYTitle("Counts")
denLen3DHist.SetYTitle("Counts")

#Compute average efficiency

print "Average efficiency for", effHist.GetName(), ":", averageOverThreshold1p2(numHist, denHist), "+/-", averageOverThreshold1p2Error(numHist, denHist)

#xIter = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
#xIterErr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

#yEff = np.array([averageOverThreshold1p2(legPixHist, denHist), averageOverThreshold1p2(legIt0Hist, denHist), averageOverThreshold1p2(legIt1Hist, denHist), averageOverThreshold1p2(legIt2Hist, denHist), averageOverThreshold1p2(leg3ReHist, denHist), averageOverThreshold1p2(leg2ReHist, denHist), averageOverThreshold1p2(numHist, denHist)])

#yEffError =  np.array([averageOverThreshold1p2Error(legPixHist, denHist), averageOverThreshold1p2Error(legIt0Hist, denHist), averageOverThreshold1p2Error(legIt1Hist, denHist), averageOverThreshold1p2Error(legIt2Hist, denHist), averageOverThreshold1p2Error(leg3ReHist, denHist), averageOverThreshold1p2Error(leg2ReHist, denHist), averageOverThreshold1p2Error(numHist, denHist)])

xIter = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
xIterErr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

yEff = np.array([averageOverThreshold1p2(legPixHist, denHist), averageOverThreshold1p2(legPixHist, denHist), averageOverThreshold1p2(legIt0Hist, denHist), averageOverThreshold1p2(legIt1Hist, denHist), averageOverThreshold1p2(legIt2Hist, denHist), averageOverThreshold1p2(numHist, denHist)])

yEffError =  np.array([averageOverThreshold1p2Error(legPixHist, denHist), averageOverThreshold1p2Error(legPixHist, denHist), averageOverThreshold1p2Error(legIt0Hist, denHist), averageOverThreshold1p2Error(legIt1Hist, denHist), averageOverThreshold1p2Error(legIt2Hist, denHist), averageOverThreshold1p2Error(numHist, denHist)])

effOverItersGraph = ROOT.TGraphErrors(6, xIter, yEff, xIterErr, yEffError)
effOverItersGraph.SetNameTitle("effOverIters", "track efficiency at trigger steps")
effOverItersGraph.GetXaxis().SetTitle("Steps")
effOverItersGraph.GetYaxis().SetTitle("Efficiency")

#Write hists useful for comparisons

effHist.SetDirectory(0)
effEtaHist.SetDirectory(0)
effPhiHist.SetDirectory(0)
effLen2DHist.SetDirectory(0)
effLen3DHist.SetDirectory(0)

effLegPixHist.SetDirectory(0)
effLegPixEtaHist.SetDirectory(0)

effLegIt0Hist.SetDirectory(0)
effLegIt0EtaHist.SetDirectory(0)
effLegIt1Hist.SetDirectory(0)
effLegIt1EtaHist.SetDirectory(0)
effLegIt2Hist.SetDirectory(0)
effLegIt2EtaHist.SetDirectory(0)
effLeg3ReHist.SetDirectory(0)
effLeg3ReEtaHist.SetDirectory(0)
effLeg2ReHist.SetDirectory(0)
effLeg2ReEtaHist.SetDirectory(0)

fakeRatePtHist.SetDirectory(0)
fakeRateEtaHist.SetDirectory(0)

outHistFile = ROOT.TFile.Open(outFileName, "RECREATE")
outHistFile.cd()

inFile.Close()

effHist.Write()
effEtaHist.Write()
effPhiHist.Write()
effLen2DHist.Write()
effLen3DHist.Write()

effLegPixHist.Write()
effLegPixEtaHist.Write()

effLegIt0Hist.Write()
effLegIt0EtaHist.Write()
effLegIt1Hist.Write()
effLegIt1EtaHist.Write()
effLegIt2Hist.Write()
effLegIt2EtaHist.Write()
effLeg3ReHist.Write()
effLeg3ReEtaHist.Write()
effLeg2ReHist.Write()
effLeg2ReEtaHist.Write()

effOverItersGraph.Write()

fakeRatePtHist.Write()
fakeRateEtaHist.Write()

outHistFile.Close()
