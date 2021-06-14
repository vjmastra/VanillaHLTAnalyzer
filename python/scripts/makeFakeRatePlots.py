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

#Read file

inFile = ROOT.TFile.Open(inFileName, "READ")
tree = inFile.Get("test/ntupleTree")

#Define hists

binning = np.array([0, 0.6, 1.2, 1.8, 2.5, 3.5, 5.5, 8.0, 12.0, 30.0])
lengthBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3])
length3DBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3, 20])

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
preSigPtHist.Sumw2()
preSigEtaHist = ROOT.TH1D("preSigEtaHist", "Matched tracks before filter;#eta;Counts", 10, -2.5, 2.5)
preSigEtaHist.Sumw2()
preSigPhiHist = ROOT.TH1D("preSigPhiHist", "Matched tracks before filter;#phi;Counts", 18, -3.2, 3.2)
preSigPhiHist.Sumw2()

preSigDxyHist   = ROOT.TH1D("preSigDxyHist",   "preSigTracks' transverse distance from mumuVtx", 50, 0, 0.1)
preSigDzHist    = ROOT.TH1D("preSigDzHist",    "preSigTracks' longitudinal distance from mumuVtx", 50, 0, 0.3)
preSigDetaHist  = ROOT.TH1D("preSigDetaHist",  "preSigTracks' eta distance from closest mu cand", 30, 0, 1.5)
preSigDphiHist  = ROOT.TH1D("preSigDphiHist",  "preSigTracks' phi distance from closest mu cand", 30, 0, 1.5)
preSigD0Hist    = ROOT.TH1D("preSigD0Hist",    "preSigTracks' transverse distance from beamspot", 30, 0, 1)
preSigD0sigHist = ROOT.TH1D("preSigD0sigHist", "preSigTracks' transverse distance signif from bs", 10, 0, 10)

preBkgPtHist  = ROOT.TH1D("preBkgPtHist",  "unmatched tracks before filter;p_{T} [GeV/c];Counts", len(binning)-1, binning)
preBkgPtHist.Sumw2()
preBkgEtaHist = ROOT.TH1D("preBkgEtaHist", "unmatched tracks before filter;#eta;Counts", 10, -2.5, 2.5)
preBkgEtaHist.Sumw2()
preBkgPhiHist = ROOT.TH1D("preBkgPhiHist", "unmatched tracks before filter;#phi;Counts", 18, -3.2, 3.2)
preBkgPhiHist.Sumw2()

preBkgDxyHist   = ROOT.TH1D("preBkgDxyHist",   "preBkgTracks' transverse distance from mumuVtx", 50, 0, 0.1)
preBkgDzHist    = ROOT.TH1D("preBkgDzHist",    "preBkgTracks' longitudinal distance from mumuVtx", 50, 0, 0.3)
preBkgDetaHist  = ROOT.TH1D("preBkgDetaHist",  "preBkgTracks' eta distance from closest mu cand", 30, 0, 1.5)
preBkgDphiHist  = ROOT.TH1D("preBkgDphiHist",  "preBkgTracks' phi distance from closest mu cand", 30, 0, 1.5)
preBkgD0Hist    = ROOT.TH1D("preBkgD0Hist",    "preBkgTracks' transverse distance from beamspot", 30, 0, 1)
preBkgD0sigHist = ROOT.TH1D("preBkgD0sigHist", "preBkgTracks' transverse distance signif from bs", 10, 0, 10)

#Fill hists

for entryNum in range(0, tree.GetEntries()):

  tree.GetEntry(entryNum)

  #Get variables

  hltJpsiMatched = getattr(tree, "muonsMatched")
  hltTrack1Matched = getattr(tree, "hltMatch_Track1")
  hltTrack2Matched = getattr(tree, "hltMatch_Track2")
  hltJpsiTrkMatched = getattr(tree, "hltMatch_Track") 

  gpuJpsiTrkMatched = getattr(tree, "gpuMatch_Track")
  gpuTrack1Matched = getattr(tree, "gpuMatch_Track1")
  gpuTrack2Matched = getattr(tree, "gpuMatch_Track2")

  hltTrigPtVec  = getattr(tree, "hltTriggerTrack_pt")
  hltTrigEtaVec = getattr(tree, "hltTriggerTrack_eta")
  hltFakePtVec  = getattr(tree, "hltTriggerFake_pt")
  hltFakeEtaVec = getattr(tree, "hltTriggerFake_eta")

  gpuTrigPtVec  = getattr(tree, "gpuTriggerTrack_pt")
  gpuTrigEtaVec = getattr(tree, "gpuTriggerTrack_eta")
  gpuFakePtVec  = getattr(tree, "gpuTriggerFake_pt")
  gpuFakeEtaVec = getattr(tree, "gpuTriggerFake_eta")

  prehltTrack1Matched = getattr(tree, "leg2RMatch_Track1")
  prehltTrack2Matched = getattr(tree, "leg2RMatch_Track2")

  pregpuTrack1Matched = getattr(tree, "iter2Match_Track1")
  pregpuTrack2Matched = getattr(tree, "iter2Match_Track2")

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

  prehltPtVec    = getattr(tree, "prefilter_legTks_pT")
  prehltEtaVec   = getattr(tree, "prefilter_legTks_eta")
  prehltPhiVec   = getattr(tree, "prefilter_legTks_phi")
  prehltdEta1Vec = getattr(tree, "prefilter_legTks_dEta1")
  prehltdPhi1Vec = getattr(tree, "prefilter_legTks_dPhi1")
  prehltdEta2Vec = getattr(tree, "prefilter_legTks_dEta2")
  prehltdPhi2Vec = getattr(tree, "prefilter_legTks_dPhi2")
  prehltDxyVec   = getattr(tree, "prefilter_legTks_dxy")
  prehltDzVec    = getattr(tree, "prefilter_legTks_dz")
  prehltD0Vec    = getattr(tree, "prefilter_legTks_d0")
  prehltD0errVec = getattr(tree, "prefilter_legTks_d0err")
  prehltBoolVec  = getattr(tree, "prefilter_legTks_bool")

  pregpuPtVec    = getattr(tree, "prefilter_newTks_pT")
  pregpuEtaVec   = getattr(tree, "prefilter_newTks_eta")
  pregpuPhiVec   = getattr(tree, "prefilter_newTks_phi")
  pregpudEta1Vec = getattr(tree, "prefilter_newTks_dEta1")
  pregpudPhi1Vec = getattr(tree, "prefilter_newTks_dPhi1")
  pregpudEta2Vec = getattr(tree, "prefilter_newTks_dEta2")
  pregpudPhi2Vec = getattr(tree, "prefilter_newTks_dPhi2")
  pregpuDxyVec   = getattr(tree, "prefilter_newTks_dxy")
  pregpuDzVec    = getattr(tree, "prefilter_newTks_dz")
  pregpuD0Vec    = getattr(tree, "prefilter_newTks_d0")
  pregpuD0errVec = getattr(tree, "prefilter_newTks_d0err")
  pregpuBoolVec  = getattr(tree, "prefilter_newTks_bool")

  #Fake rate studies (legacy)

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

fakePtHist.Divide(trigPtHist)
fakeEtaHist.Divide(trigEtaHist)

#Set directories

fakePtHist.SetDirectory(0)
fakeEtaHist.SetDirectory(0)

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

line = ROOT.TLine(1.2, 0, 1.2, 1.05)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

fakePtHist.SetLineColor(1)
fakeEtaHist.SetLineColor(1)

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

fakePtHist.SetLineWidth(2)
fakeEtaHist.SetLineWidth(2)

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

fakePtHist.GetYaxis().SetTitle("Fake rate")
fakeEtaHist.GetYaxis().SetTitle("Fake rate")

#histLegend = ROOT.TLegend(0.65, 0.63, 0.89, 0.73)
#histLegend.SetTextSize(0.03)
#histLegend.AddEntry(myPixSizeHist, "new pixelTracks")
#histLegend.AddEntry(regPixSizeHist, "reg. pixelTracks")

#define canvas and print

canvas = ROOT.TCanvas("canvas")
canvas.cd()

canvas.Print(outFileName+"[")

fakePtHist.Draw("")
canvas.Print(outFileName)

fakeEtaHist.Draw("")
canvas.Print(outFileName)

preBkgPtHist.Draw("")
preSigPtHist.Draw("same")
canvas.Print(outFileName)

preBkgEtaHist.Draw("")
preSigEtaHist.Draw("same")
canvas.Print(outFileName)

preBkgPhiHist.Draw("")
preSigPhiHist.Draw("same")
canvas.Print(outFileName)

preBkgDetaHist.Draw("")
preSigDetaHist.Draw("same")
canvas.Print(outFileName)

preBkgDphiHist.Draw("")
preSigDphiHist.Draw("same")
canvas.Print(outFileName)

preBkgDxyHist.Draw("")
preSigDxyHist.Draw("same")
canvas.Print(outFileName)

preBkgDzHist.Draw("")
preSigDzHist.Draw("same")
canvas.Print(outFileName)

preBkgD0Hist.Draw("")
preSigD0Hist.Draw("same")
canvas.Print(outFileName)

preBkgD0sigHist.Draw("")
preSigD0sigHist.Draw("same")
canvas.Print(outFileName)

canvas.Print(outFileName+"]")

