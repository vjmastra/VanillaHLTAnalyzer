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

myPixSizeHist = ROOT.TH1D("myPixSizeHist", "Size of pixelTracks collection", 30, 0, 3000)
myPixDxyHist  = ROOT.TH1D("myPixDxyHist",  "pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
myPixDzHist   = ROOT.TH1D("myPixDzHist",   "pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
myPixDetaHist = ROOT.TH1D("myPixDetaHist", "pixelTracks' eta distance from mu cands", 50, 0, 2)
myPixDphiHist = ROOT.TH1D("myPixDphiHist", "pixelTracks' phi distance from mu cands", 50, 0, 2)
myPixChi2nHist = ROOT.TH1D("myPixChi2nHist", "pixelTrakcs' reduced chisquare", len(lengthBinning), lengthBinning)

myPixDxyDzHist    = ROOT.TH2D("myPixVtxDist",      "pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
myPixDetaDphiHist = ROOT.TH2D("myPixDetaDphiHist", "pixelTracks' angular distance from mu cands", 50, 0, 2, 50, 0, 2)

regPixHist = ROOT.TH1D("regPixHist", "Regional pixel tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
regPixHist.Sumw2()
regPixEtaHist = ROOT.TH1D("regPixEtaHist", "Regional pixel tracks;#eta;Counts", 10, -2.5, +2.5)
regPixEtaHist.Sumw2()

regPixSizeHist = ROOT.TH1D("regPixSizeHist", "Size of regpixelTracks collection", 30, 0, 3000)
regPixDxyHist  = ROOT.TH1D("regPixDxyHist",  "reg pixelTracks' transverse distance from mumuVtx", 100, 0, 0.2)
regPixDzHist   = ROOT.TH1D("regPixDzHist",   "reg pixelTracks' longitudinal distance from mumuVtx", 100, 0, 0.5)
regPixDetaHist = ROOT.TH1D("regPixDetaHist", "reg pixelTracks' eta distance from mu cands", 50, 0, 2)
regPixDphiHist = ROOT.TH1D("regPixDphiHist", "reg pixelTracks' phi distance from mu cands", 50, 0, 2)
regPixChi2nHist = ROOT.TH1D("mregixChi2nHist", "reg pixelTracks' reduced chisquare", len(lengthBinning), lengthBinning)

regPixDxyDzHist    = ROOT.TH2D("regPixVtxDist",      "reg pixelTracks' distance from mumuVtx", 100, 0, 0.2, 100, 0, 0.5)
regPixDetaDphiHist = ROOT.TH2D("regPixDetaDphiHist", "reg pixelTracks' angular distance from mu cands", 50, 0, 2, 50, 0, 2)

ctf0Hist = ROOT.TH1D("ctf0Hist", "Ctf0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf0Hist.Sumw2()
ctf0EtaHist = ROOT.TH1D("ctf0EtaHist", "Ctf0 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf0EtaHist.Sumw2()

iter0Hist = ROOT.TH1D("iter0Hist", "Iter0 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter0Hist.Sumw2()
iter0EtaHist = ROOT.TH1D("iter0EtaHist", "Iter0 tracks (;#eta;Counts", 10, -2.5, +2.5)
iter0EtaHist.Sumw2()

ctf1Hist = ROOT.TH1D("ctf1Hist", "Ctf1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf1Hist.Sumw2()
ctf1EtaHist = ROOT.TH1D("ctf1EtaHist", "Ctf1 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf1EtaHist.Sumw2()

iter1Hist = ROOT.TH1D("iter1Hist", "Iter1 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter1Hist.Sumw2()
iter1EtaHist = ROOT.TH1D("iter1EtaHist", "Iter1 tracks (;#eta;Counts", 10, -2.5, +2.5)
iter1EtaHist.Sumw2()

ctf2Hist = ROOT.TH1D("ctf2Hist", "Ctf2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
ctf2Hist.Sumw2()
ctf2EtaHist = ROOT.TH1D("ctf2EtaHist", "Ctf2 tracks;#eta;Counts", 10, -2.5, +2.5)
ctf2EtaHist.Sumw2()

iter2Hist = ROOT.TH1D("iter2Hist", "Iter2 tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
iter2Hist.Sumw2()
iter2EtaHist = ROOT.TH1D("iter2EtaHist", "Iter2 tracks;#eta;Counts", 10, -2.5, +2.5)
iter2EtaHist.Sumw2()

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

  myPixTrack2dxy   = getattr(tree, "myPix_Track2_dxy")
  myPixTrack2dz    = getattr(tree, "myPix_Track2_dz")
  myPixTrack2dEta1 = getattr(tree, "myPix_Track2_dEtaMu1")
  myPixTrack2dPhi1 = getattr(tree, "myPix_Track2_dPhiMu1")
  myPixTrack2dEta2 = getattr(tree, "myPix_Track2_dEtaMu2")
  myPixTrack2dPhi2 = getattr(tree, "myPix_Track2_dPhiMu2")
  myPixTrack2chi2n = getattr(tree, "myPix_Track2_normChi2")

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

  regPixSize = getattr(tree, "reg_collSize")

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
      myPixDxyHist.Fill(myPixTrack1dxy)
      myPixDzHist.Fill(myPixTrack1dz)
      myPixDetaHist.Fill(myPixTrack1dEta1)
      myPixDetaHist.Fill(myPixTrack1dEta2)
      myPixDphiHist.Fill(myPixTrack1dPhi1)
      myPixDphiHist.Fill(myPixTrack1dPhi2)
      myPixChi2nHist.Fill(myPixTrack1chi2n)
      myPixDxyDzHist.Fill(myPixTrack1dxy, myPixTrack1dz)
      myPixDetaDphiHist.Fill(myPixTrack1dEta1, myPixTrack1dPhi1)
      myPixDetaDphiHist.Fill(myPixTrack1dEta2, myPixTrack1dPhi2)
    if (myPixTrack2Matched):
      myPixHist.Fill(genTrack2Pt)
      myPixEtaHist.Fill(genTrack2eta)
      myPixDxyHist.Fill(myPixTrack2dxy)
      myPixDzHist.Fill(myPixTrack2dz)
      myPixDetaHist.Fill(myPixTrack2dEta1)
      myPixDetaHist.Fill(myPixTrack2dEta2)
      myPixDphiHist.Fill(myPixTrack2dPhi1)
      myPixDphiHist.Fill(myPixTrack2dPhi2)
      myPixChi2nHist.Fill(myPixTrack2chi2n)
      myPixDxyDzHist.Fill(myPixTrack2dxy, myPixTrack2dz)
      myPixDetaDphiHist.Fill(myPixTrack2dEta1, myPixTrack2dPhi1)
      myPixDetaDphiHist.Fill(myPixTrack2dEta2, myPixTrack2dPhi2)

  if (hltJpsiMatched):
    regPixSizeHist.Fill(regPixSize)
    if (regPixTrack1Matched):
      regPixHist.Fill(genTrack1Pt)
      regPixEtaHist.Fill(genTrack1eta)
      regPixDxyHist.Fill(regPixTrack1dxy)
      regPixDzHist.Fill(regPixTrack1dz)
      regPixDetaHist.Fill(regPixTrack1dEta1)
      regPixDetaHist.Fill(regPixTrack1dEta2)
      regPixDphiHist.Fill(regPixTrack1dPhi1)
      regPixDphiHist.Fill(regPixTrack1dPhi2)
      regPixChi2nHist.Fill(regPixTrack1chi2n)
      regPixDxyDzHist.Fill(regPixTrack1dxy, regPixTrack1dz)
      regPixDetaDphiHist.Fill(regPixTrack1dEta1, regPixTrack1dPhi1)
      regPixDetaDphiHist.Fill(regPixTrack1dEta2, regPixTrack1dPhi2)
    if (regPixTrack2Matched):
      regPixHist.Fill(genTrack2Pt)
      regPixEtaHist.Fill(genTrack2eta)
      regPixDxyHist.Fill(regPixTrack2dxy)
      regPixDzHist.Fill(regPixTrack2dz)
      regPixDetaHist.Fill(regPixTrack2dEta1)
      regPixDetaHist.Fill(regPixTrack2dEta2)
      regPixDphiHist.Fill(regPixTrack2dPhi1)
      regPixDphiHist.Fill(regPixTrack2dPhi2)
      regPixChi2nHist.Fill(regPixTrack2chi2n)
      regPixDxyDzHist.Fill(regPixTrack2dxy, regPixTrack2dz)
      regPixDetaDphiHist.Fill(regPixTrack2dEta1, regPixTrack2dPhi1)
      regPixDetaDphiHist.Fill(regPixTrack2dEta2, regPixTrack2dPhi2)

    if (ctf0Track1Matched):
      ctf0Hist.Fill(genTrack1Pt)
      ctf0EtaHist.Fill(genTrack1eta)
    if (ctf0Track2Matched):
      ctf0Hist.Fill(genTrack2Pt)
      ctf0EtaHist.Fill(genTrack2eta)
    if (iter0Track1Matched):
      iter0Hist.Fill(genTrack1Pt)
      iter0EtaHist.Fill(genTrack1eta)
    if (iter0Track2Matched):
      iter0Hist.Fill(genTrack2Pt)
      iter0EtaHist.Fill(genTrack2eta)

    if (ctf1Track1Matched):
      ctf1Hist.Fill(genTrack1Pt)
      ctf1EtaHist.Fill(genTrack1eta)
    if (ctf1Track2Matched):
      ctf1Hist.Fill(genTrack2Pt)
      ctf1EtaHist.Fill(genTrack2eta)
    if (iter1Track1Matched):
      iter1Hist.Fill(genTrack1Pt)
      iter1EtaHist.Fill(genTrack1eta)
    if (iter1Track2Matched):
      iter1Hist.Fill(genTrack2Pt)
      iter1EtaHist.Fill(genTrack2eta)

    if (ctf2Track1Matched):
      ctf2Hist.Fill(genTrack1Pt)
      ctf2EtaHist.Fill(genTrack1eta)
    if (ctf2Track2Matched):
      ctf2Hist.Fill(genTrack2Pt)
      ctf2EtaHist.Fill(genTrack2eta)
    if (iter2Track1Matched):
      iter2Hist.Fill(genTrack1Pt)
      iter2EtaHist.Fill(genTrack1eta)
    if (iter2Track2Matched):
      iter2Hist.Fill(genTrack2Pt)
      iter2EtaHist.Fill(genTrack2eta)

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

#new pixel tracks

effMyPixHist = ROOT.TEfficiency(myPixHist, denHist)
effMyPixHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effMyPixHist.SetConfidenceLevel(0.68)
effMyPixHist.SetNameTitle("effMyPixHist", "pixel Tracks efficiency")

effMyPixEtaHist = ROOT.TEfficiency(myPixEtaHist, denEtaHist)
effMyPixEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effMyPixEtaHist.SetConfidenceLevel(0.68)
effMyPixEtaHist.SetNameTitle("effMyPixEtaHist", "pixel Tracks efficiency")

#regional pixel tracks

effRegPixHist = ROOT.TEfficiency(regPixHist, denHist)
effRegPixHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effRegPixHist.SetConfidenceLevel(0.68)
effRegPixHist.SetNameTitle("effRegPixHist", "reg pixel Tracks efficiency")

effRegPixEtaHist = ROOT.TEfficiency(regPixEtaHist, denEtaHist)
effRegPixEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effRegPixEtaHist.SetConfidenceLevel(0.68)
effRegPixEtaHist.SetNameTitle("effRegPixEtaHist", "reg pixel Tracks efficiency")

#ctf0

effCtf0Hist = ROOT.TEfficiency(ctf0Hist, denHist)
effCtf0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf0Hist.SetConfidenceLevel(0.68)
effCtf0Hist.SetNameTitle("effCtf0Hist", "ctf0 Track efficiency")

effCtf0EtaHist = ROOT.TEfficiency(ctf0EtaHist, denEtaHist)
effCtf0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf0EtaHist.SetConfidenceLevel(0.68)
effCtf0EtaHist.SetNameTitle("effCtf0EtaHist", "ctf0 Track efficiency")

#iter0

effIter0Hist = ROOT.TEfficiency(iter0Hist, denHist)
effIter0Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter0Hist.SetConfidenceLevel(0.68)
effIter0Hist.SetNameTitle("effIter0Hist", "iter0 Track efficiency")

effIter0EtaHist = ROOT.TEfficiency(iter0EtaHist, denEtaHist)
effIter0EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter0EtaHist.SetConfidenceLevel(0.68)
effIter0EtaHist.SetNameTitle("effIter0EtaHist", "iter0 Track efficiency")

#ctf1

effCtf1Hist = ROOT.TEfficiency(ctf1Hist, denHist)
effCtf1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf1Hist.SetConfidenceLevel(0.68)
effCtf1Hist.SetNameTitle("effCtf1Hist", "ctf1 Track efficiency")

effCtf1EtaHist = ROOT.TEfficiency(ctf1EtaHist, denEtaHist)
effCtf1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf1EtaHist.SetConfidenceLevel(0.68)
effCtf1EtaHist.SetNameTitle("effCtf1EtaHist", "ctf1 Track efficiency")

#iter1

effIter1Hist = ROOT.TEfficiency(iter1Hist, denHist)
effIter1Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter1Hist.SetConfidenceLevel(0.68)
effIter1Hist.SetNameTitle("effIter1Hist", "iter1 Track efficiency")

effIter1EtaHist = ROOT.TEfficiency(iter1EtaHist, denEtaHist)
effIter1EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter1EtaHist.SetConfidenceLevel(0.68)
effIter1EtaHist.SetNameTitle("effIter1EtaHist", "iter1 Track efficiency")

#ctf2

effCtf2Hist = ROOT.TEfficiency(ctf2Hist, denHist)
effCtf2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf2Hist.SetConfidenceLevel(0.68)
effCtf2Hist.SetNameTitle("effCtf2Hist", "ctf2 Track efficiency")

effCtf2EtaHist = ROOT.TEfficiency(ctf2EtaHist, denEtaHist)
effCtf2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effCtf2EtaHist.SetConfidenceLevel(0.68)
effCtf2EtaHist.SetNameTitle("effCtf2EtaHist", "ctf2 Track efficiency")

#iter2

effIter2Hist = ROOT.TEfficiency(iter2Hist, denHist)
effIter2Hist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter2Hist.SetConfidenceLevel(0.68)
effIter2Hist.SetNameTitle("effIter2Hist", "iter2 Track efficiency")

effIter2EtaHist = ROOT.TEfficiency(iter2EtaHist, denEtaHist)
effIter2EtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
effIter2EtaHist.SetConfidenceLevel(0.68)
effIter2EtaHist.SetNameTitle("effIter2EtaHist", "iter2 Track efficiency")

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

xIter = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
xIterErr = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

yEff = np.array([averageOverThreshold1p2(myPixHist, denHist), averageOverThreshold1p2(regPixHist, denHist), averageOverThreshold1p2(iter0Hist, denHist), averageOverThreshold1p2(iter1Hist, denHist), averageOverThreshold1p2(iter2Hist, denHist), averageOverThreshold1p2(numHist, denHist)])

yEffError =  np.array([averageOverThreshold1p2Error(myPixHist, denHist), averageOverThreshold1p2Error(regPixHist, denHist), averageOverThreshold1p2Error(iter0Hist, denHist), averageOverThreshold1p2Error(iter1Hist, denHist), averageOverThreshold1p2Error(iter2Hist, denHist), averageOverThreshold1p2Error(numHist, denHist)])

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

effMyPixHist.SetDirectory(0)
effMyPixEtaHist.SetDirectory(0)
effRegPixHist.SetDirectory(0)
effRegPixEtaHist.SetDirectory(0)

effCtf0Hist.SetDirectory(0)
effCtf0EtaHist.SetDirectory(0)
effIter0Hist.SetDirectory(0)
effIter0EtaHist.SetDirectory(0)
effCtf1Hist.SetDirectory(0)
effCtf1EtaHist.SetDirectory(0)
effIter1Hist.SetDirectory(0)
effIter1EtaHist.SetDirectory(0)
effCtf2Hist.SetDirectory(0)
effCtf2EtaHist.SetDirectory(0)
effIter2Hist.SetDirectory(0)
effIter2EtaHist.SetDirectory(0)

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

effMyPixHist.Write()
effMyPixEtaHist.Write()
effRegPixHist.Write()
effRegPixEtaHist.Write()

effCtf0Hist.Write()
effCtf0EtaHist.Write()
effIter0Hist.Write()
effIter0EtaHist.Write()
effCtf1Hist.Write()
effCtf1EtaHist.Write()
effIter1Hist.Write()
effIter1EtaHist.Write()
effCtf2Hist.Write()
effCtf2EtaHist.Write()
effIter2Hist.Write()
effIter2EtaHist.Write()

effOverItersGraph.Write()

fakeRatePtHist.Write()
fakeRateEtaHist.Write()

outHistFile.Close()
