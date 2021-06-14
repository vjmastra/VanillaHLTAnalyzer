import ROOT
import sys
import numpy as np

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

#Read file

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("test/ntupleTree")

#Define hists

binning         = np.array([0, 0.6, 1.2, 1.8, 2.5, 3.5, 5.5, 8.0, 12.0, 30.0])
lengthBinning   = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3])
length3DBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3, 20])
d0sigBinning    = np.array([0, 1., 2., 5., 10., 20., 35., 50., 65., 80., 100.])

#After filter

fakePtHist = ROOT.TH1D("fakePtHist", "Unmatched trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
fakePtHist.Sumw2()
fakeEtaHist = ROOT.TH1D("fakeEtaHist", "Unmatched trigger tracks;#eta;Counts", 10, -2.5, 2.5)
fakeEtaHist.Sumw2()

trigPtHist = ROOT.TH1D("trigPtHist", "All trigger tracks;p_{T} [GeV/c];Counts", len(binning)-1, binning)
trigPtHist.Sumw2()
trigEtaHist = ROOT.TH1D("trigEtaHist", "All trigger tracks;#eta;Counts", 10, -2.5, 2.5)
trigEtaHist.Sumw2()

#Before filter

preSigPtHist  = ROOT.TH1D("preSigPtHist",  "Matched tracks before filter;p_{T} [GeV/c];Counts", len(binning)-1, binning)
preSigEtaHist = ROOT.TH1D("preSigEtaHist", "Matched tracks before filter;#eta;Counts", 10, -2.5, 2.5)
preSigPhiHist = ROOT.TH1D("preSigPhiHist", "Matched tracks before filter;#phi;Counts", 18, -3.2, 3.2)

preSigDxyHist   = ROOT.TH1D("preSigDxyHist",   "preSigTracks' transverse distance from mumuVtx;dxy [cm];Counts", 50, 0, 0.1)
preSigDzHist    = ROOT.TH1D("preSigDzHist",    "preSigTracks' longitudinal distance from mumuVtx;dz [cm];Counts", 50, 0, 0.3)
preSigDetaHist  = ROOT.TH1D("preSigDetaHist",  "preSigTracks' eta distance from closest mu cand;#Delta#eta;Counts", 30, 0, 1.5)
preSigDphiHist  = ROOT.TH1D("preSigDphiHist",  "preSigTracks' phi distance from closest mu cand;#Delta#phi;Counts", 30, 0, 1.5)
preSigD0Hist    = ROOT.TH1D("preSigD0Hist",    "preSigTracks' transverse distance from beamspot;d0 [cm];Counts", 20, 0, 0.5)
preSigD0sigHist = ROOT.TH1D("preSigD0sigHist", "preSigTracks' transverse distance signif from bs;d0 significance;Counts", len(d0sigBinning)-1, d0sigBinning)

preBkgPtHist  = ROOT.TH1D("preBkgPtHist",  "HLT tracks before filter;p_{T} [GeV/c];Counts (norm.)", len(binning)-1, binning)
preBkgEtaHist = ROOT.TH1D("preBkgEtaHist", "HLT tracks before filter;#eta;Counts (norm.)", 10, -2.5, 2.5)
preBkgPhiHist = ROOT.TH1D("preBkgPhiHist", "HLT tracks before filter;#phi;Counts (norm.)", 18, -3.2, 3.2)

preBkgDxyHist   = ROOT.TH1D("preBkgDxyHist",   "HLT tracks before filter transverse distance from #mu#muVtx;d_{xy} [cm];Counts (norm.)", 50, 0, 0.1)
preBkgDzHist    = ROOT.TH1D("preBkgDzHist",    "HLT tracks before filter longitudinal distance from #mu#muVtx;d_{z} [cm];Counts (norm.)", 50, 0, 0.3)
preBkgDetaHist  = ROOT.TH1D("preBkgDetaHist",  "HLT tracks before filter #eta distance from closest #mu cand;#Delta#eta;Counts (norm.)", 30, 0, 1.5)
preBkgDphiHist  = ROOT.TH1D("preBkgDphiHist",  "HLT tracks before fitler #phi distance from closest #mu cand;#Delta#phi;Counts (norm.)", 30, 0, 1.5)
preBkgD0Hist    = ROOT.TH1D("preBkgD0Hist",    "HLT tracks before filter transverse distance from beamspot;d_{0} [cm];Counts (norm.)", 20, 0, 0.5)
preBkgD0sigHist = ROOT.TH1D("preBkgD0sigHist", "HLT tracks before filter transverse distance signif from BS;d_{0}/#sigma(d_{0});Counts (norm.)", len(d0sigBinning)-1, d0sigBinning)

#Fill hists

for entryNum in range(0, tree.GetEntries()):

  tree.GetEntry(entryNum)

  #Get variables

  hltJpsiMatched = getattr(tree, "muonsMatched")
  hltTrack1Matched = getattr(tree, "gpuMatch_Track1")
  hltTrack2Matched = getattr(tree, "gpuMatch_Track2")
  hltJpsiTrkMatched = getattr(tree, "gpuMatch_Track") 

  hltTrigPtVec  = getattr(tree, "gpuTriggerTrack_pt")
  hltTrigEtaVec = getattr(tree, "gpuTriggerTrack_eta")
  hltFakePtVec  = getattr(tree, "gpuTriggerFake_pt")
  hltFakeEtaVec = getattr(tree, "gpuTriggerFake_eta")

  prehltTrack1Matched = getattr(tree, "iter2Match_Track1")
  prehltTrack2Matched = getattr(tree, "iter2Match_Track2")

  prehltPtVec    = getattr(tree, "prefilter_newTks_pT")
  prehltEtaVec   = getattr(tree, "prefilter_newTks_eta")
  prehltPhiVec   = getattr(tree, "prefilter_newTks_phi")
  prehltdEta1Vec = getattr(tree, "prefilter_newTks_dEta1")
  prehltdPhi1Vec = getattr(tree, "prefilter_newTks_dPhi1")
  prehltdEta2Vec = getattr(tree, "prefilter_newTks_dEta2")
  prehltdPhi2Vec = getattr(tree, "prefilter_newTks_dPhi2")
  prehltDxyVec   = getattr(tree, "prefilter_newTks_dxy")
  prehltDzVec    = getattr(tree, "prefilter_newTks_dz")
  prehltD0Vec    = getattr(tree, "prefilter_newTks_d0")
  prehltD0errVec = getattr(tree, "prefilter_newTks_d0err")
  prehltBoolVec  = getattr(tree, "prefilter_newTks_bool")

  #Fake rate studies (new)

  multi = zip(prehltPtVec, prehltEtaVec, prehltPhiVec, prehltdEta1Vec, prehltdPhi1Vec, prehltdEta2Vec, prehltdPhi2Vec, prehltDxyVec, prehltDzVec, prehltD0Vec, prehltD0errVec, prehltBoolVec)
  deta = 0
  dphi = 0

  if (hltJpsiMatched):
  
    for pt, eta, phi, deta1, dphi1, deta2, dphi2, dxy, dz, d0, d0err, flag in multi:

      val1 = deta1**2 + dphi1**2
      val2 = deta2**2 + dphi2**2
      muFlag = closerMuon(val1, val2)
      if (muFlag == 1):
        deta = deta1
        dphi = dphi1
      elif (muFlag == 2): 
        deta = deta2
        dphi = dphi2

      if (d0 < 0):
        d0 = -d0
      d0sig = 0
      if (d0err == 0):
        d0sig = 0
      else:
        d0sig = (d0/d0err)

      if (flag == 0):
        preBkgPtHist.Fill(pt)
        preBkgEtaHist.Fill(eta)
        preBkgPhiHist.Fill(phi)
        preBkgDetaHist.Fill(deta)
        preBkgDphiHist.Fill(dphi)
        preBkgDxyHist.Fill(dxy)
        preBkgDzHist.Fill(dz)
        preBkgD0Hist.Fill(d0)
        preBkgD0sigHist.Fill(d0sig)
      if (flag == 1):
        preSigPtHist.Fill(pt)
        preSigEtaHist.Fill(eta)
        preSigPhiHist.Fill(phi)
        preSigDetaHist.Fill(deta)
        preSigDphiHist.Fill(dphi)
        preSigDxyHist.Fill(dxy)
        preSigDzHist.Fill(dz)
        preSigD0Hist.Fill(d0)
        preSigD0sigHist.Fill(d0sig)   

    if (hltJpsiTrkMatched):
      for trk in hltTrigPtVec:
        trigPtHist.Fill(trk)
      for trk in hltTrigEtaVec:
        trigEtaHist.Fill(trk)
      for trk in hltFakePtVec:
        fakePtHist.Fill(trk)
      for trk in hltFakeEtaVec:
        fakeEtaHist.Fill(trk)

#Fake rate 

fakeRatePtHist = ROOT.TEfficiency(fakePtHist, trigPtHist)
fakeRatePtHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
fakeRatePtHist.SetConfidenceLevel(0.68)
fakeRatePtHist.SetNameTitle("fakeRatePtHist", "Fake rate distribution")

fakeRateEtaHist = ROOT.TEfficiency(fakeEtaHist, trigEtaHist)
fakeRateEtaHist.SetStatisticOption(ROOT.TEfficiency.kFCP)
fakeRateEtaHist.SetConfidenceLevel(0.68)
fakeRateEtaHist.SetNameTitle("fakeRateEtaHist", "Fake rate distribution")

num = fakePtHist.Integral()
den = trigPtHist.Integral()

if (den == 0):
  ratio = 0
else:
  ratio = num/den

print "fake rate:", ratio

#Set directories

fakeRatePtHist.SetDirectory(0)
fakeRateEtaHist.SetDirectory(0)

preBkgPtHist.SetDirectory(0)
preBkgEtaHist.SetDirectory(0)
preBkgPhiHist.SetDirectory(0)
preBkgDetaHist.SetDirectory(0)
preBkgDphiHist.SetDirectory(0)
preBkgDxyHist.SetDirectory(0)
preBkgDzHist.SetDirectory(0)
preBkgD0Hist.SetDirectory(0)
preBkgD0sigHist.SetDirectory(0)

preSigPtHist.SetDirectory(0)
preSigEtaHist.SetDirectory(0)
preSigPhiHist.SetDirectory(0)
preSigDetaHist.SetDirectory(0)
preSigDphiHist.SetDirectory(0)
preSigDxyHist.SetDirectory(0)
preSigDzHist.SetDirectory(0)
preSigD0Hist.SetDirectory(0)
preSigD0sigHist.SetDirectory(0)

inFile.Close()

#set graphics

emptyPtHist = ROOT.TH1D("emptyHist", "title;p_{T} [GeV/c];Fake rate", 10, 0, 30)
emptyEtaHist = ROOT.TH1D("emptyEtaHist", "title;#eta;Fake rate", 10, -2.5, 2.5)

emptyPtHist.GetYaxis().SetRangeUser(0., 1.0)
emptyEtaHist.GetYaxis().SetRangeUser(0, 1.0)

emptyPtHist.SetStats(ROOT.kFALSE)
emptyEtaHist.SetStats(ROOT.kFALSE)

line = ROOT.TLine(1.2, 0, 1.2, 1.0)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

fakeRatePtHist.SetLineColor(1)
fakeRateEtaHist.SetLineColor(1)

preBkgPtHist.SetLineColor(2)
preBkgEtaHist.SetLineColor(2)
preBkgPhiHist.SetLineColor(2)
preBkgDetaHist.SetLineColor(2)
preBkgDphiHist.SetLineColor(2)
preBkgDxyHist.SetLineColor(2)
preBkgDzHist.SetLineColor(2)
preBkgD0Hist.SetLineColor(2)
preBkgD0sigHist.SetLineColor(2)

preSigPtHist.SetLineColor(1)
preSigEtaHist.SetLineColor(1)
preSigPhiHist.SetLineColor(1)
preSigDetaHist.SetLineColor(1)
preSigDphiHist.SetLineColor(1)
preSigDxyHist.SetLineColor(1)
preSigDzHist.SetLineColor(1)
preSigD0Hist.SetLineColor(1)
preSigD0sigHist.SetLineColor(1)

fakeRatePtHist.SetLineWidth(2)
fakeRateEtaHist.SetLineWidth(2)

preBkgPtHist.SetLineWidth(2)
preBkgEtaHist.SetLineWidth(2)
preBkgPhiHist.SetLineWidth(2)
preBkgDetaHist.SetLineWidth(2)
preBkgDphiHist.SetLineWidth(2)
preBkgDxyHist.SetLineWidth(2)
preBkgDzHist.SetLineWidth(2)
preBkgD0Hist.SetLineWidth(2)
preBkgD0sigHist.SetLineWidth(2)

preSigPtHist.SetLineWidth(2)
preSigEtaHist.SetLineWidth(2)
preSigPhiHist.SetLineWidth(2)
preSigDetaHist.SetLineWidth(2)
preSigDphiHist.SetLineWidth(2)
preSigDxyHist.SetLineWidth(2)
preSigDzHist.SetLineWidth(2)
preSigD0Hist.SetLineWidth(2)
preSigD0sigHist.SetLineWidth(2)

preBkgPtHist.SetStats(ROOT.kFALSE)
preBkgEtaHist.SetStats(ROOT.kFALSE)
preBkgPhiHist.SetStats(ROOT.kFALSE)
preBkgDetaHist.SetStats(ROOT.kFALSE)
preBkgDphiHist.SetStats(ROOT.kFALSE)
preBkgDxyHist.SetStats(ROOT.kFALSE)
preBkgDzHist.SetStats(ROOT.kFALSE)
preBkgD0Hist.SetStats(ROOT.kFALSE)
preBkgD0sigHist.SetStats(ROOT.kFALSE)

preBkgPtHist.Scale(   1/preBkgPtHist.Integral())
preBkgEtaHist.Scale(  1/preBkgEtaHist.Integral())
preBkgPhiHist.Scale(  1/preBkgPhiHist.Integral())
preBkgDetaHist.Scale( 1/preBkgDetaHist.Integral())
preBkgDphiHist.Scale( 1/preBkgDphiHist.Integral())
preBkgDxyHist.Scale(  1/preBkgDxyHist.Integral())
preBkgDzHist.Scale(   1/preBkgDzHist.Integral())
preBkgD0Hist.Scale(   1/preBkgD0Hist.Integral())
preBkgD0sigHist.Scale(1/preBkgD0sigHist.Integral())

preSigPtHist.Scale(   1/preSigPtHist.Integral())
preSigEtaHist.Scale(  1/preSigEtaHist.Integral())
preSigPhiHist.Scale(  1/preSigPhiHist.Integral())
preSigDetaHist.Scale( 1/preSigDetaHist.Integral())
preSigDphiHist.Scale( 1/preSigDphiHist.Integral())
preSigDxyHist.Scale(  1/preSigDxyHist.Integral())
preSigDzHist.Scale(   1/preSigDzHist.Integral())
preSigD0Hist.Scale(   1/preSigD0Hist.Integral())
preSigD0sigHist.Scale(1/preSigD0sigHist.Integral())

preBkgPtHist.GetYaxis().SetRangeUser(0, 1)
preBkgEtaHist.GetYaxis().SetRangeUser(0, 1)
preBkgPhiHist.GetYaxis().SetRangeUser(0, 1)
preBkgDetaHist.GetYaxis().SetRangeUser(0, 1)
preBkgDphiHist.GetYaxis().SetRangeUser(0, 1)
preBkgDxyHist.GetYaxis().SetRangeUser(0, 1)
preBkgDzHist.GetYaxis().SetRangeUser(0, 1)
preBkgD0Hist.GetYaxis().SetRangeUser(0, 1)
preBkgD0sigHist.GetYaxis().SetRangeUser(0, 1)

histLegend = ROOT.TLegend(0.65, 0.73, 0.89, 0.83)
histLegend.SetTextSize(0.03)
histLegend.AddEntry(preBkgPtHist, "unmatched tracks")
histLegend.AddEntry(preSigPtHist, "gen-matched tracks")

#define canvas and print

canvas = ROOT.TCanvas("canvas")
canvas.cd()

canvas.Print(outFileName+"[")

emptyPtHist.SetTitle(fakeRatePtHist.GetTitle())
emptyPtHist.Draw("")
line.Draw("same")
fakeRatePtHist.Draw("same")
canvas.Print(outFileName)

emptyEtaHist.SetTitle(fakeRateEtaHist.GetTitle())
emptyEtaHist.Draw("")
fakeRateEtaHist.Draw("same")
canvas.Print(outFileName)

preBkgPtHist.Draw("")
preSigPtHist.Draw("same")
histLegend.Draw("same")
line.Draw("same")
canvas.Print(outFileName)

preBkgEtaHist.Draw("")
preSigEtaHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgPhiHist.Draw("")
preSigPhiHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgDetaHist.Draw("")
preSigDetaHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgDphiHist.Draw("")
preSigDphiHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgDxyHist.Draw("")
preSigDxyHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgDzHist.Draw("")
preSigDzHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgD0Hist.Draw("")
preSigD0Hist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

preBkgD0sigHist.Draw("")
preSigD0sigHist.Draw("same")
histLegend.Draw("same")
canvas.Print(outFileName)

canvas.Print(outFileName+"]")

