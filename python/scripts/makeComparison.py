import ROOT
import sys
import numpy as np

def checkEffHist(effHist, name):
  if not effHist:
    print "Failed to get", name, "histogram"
    sys.exit(1)
  effHist.SetDirectory(0)

if len(sys.argv) != 4:
  print "USAGE: %s <input file1> <input file2> <output file>"%(sys.argv[0])
  sys.exit(1)

inFileName1 = sys.argv[1]
inFileName2 = sys.argv[2]
outFileName = sys.argv[3]

print "Reading from", inFileName1, "and", inFileName2, "and writing to", outFileName

text1 = inFileName1.replace(".root", "")
text1 = text1.replace("efficiencyDirectory/efficiency_", "")
text2 = inFileName2.replace(".root", "")
text2 = text2.replace("efficiencyDirectory/efficiency_", "")

text1 = "Legacy menu&tracks"
text1 = "Run-2 pixelTracks"
text2 = "Patatrack pixeltracks"

#Read histograms

inFile1 = ROOT.TFile.Open(inFileName1, "READ")
inFile2 = ROOT.TFile.Open(inFileName2, "READ")

effHist1 = inFile1.Get("effHist")
effHist2 = inFile2.Get("effHist")
effEtaHist1 = inFile1.Get("effEtaHist")
effEtaHist2 = inFile2.Get("effEtaHist")
effPhiHist1 = inFile1.Get("effPhiHist")
effPhiHist2 = inFile2.Get("effPhiHist")
effLen2DHist1 = inFile1.Get("effLen2DHist")
effLen2DHist2 = inFile2.Get("effLen2DHist")
effLen3DHist1 = inFile1.Get("effLen3DHist")
effLen3DHist2 = inFile2.Get("effLen3DHist")

fakeRateHist1 = inFile1.Get("fakeRatePtHist")
fakeRateHist2 = inFile2.Get("fakeRatePtHist")

effOverIters1 = inFile1.Get("effOverIters")
effOverIters2 = inFile2.Get("effOverIters")

#effLegPixHist1 = inFile1.Get("effLegPixHist")
#effMyPixHist2 = inFile2.Get("effMyPixHist")
#effRegPixHist2 = inFile2.Get("effRegPixHist")

#effLegPixEtaHist1 = inFile1.Get("effLegPixEtaHist")
#effMyPixEtaHist2 = inFile2.Get("effMyPixEtaHist")
#effRegPixEtaHist2 = inFile2.Get("effRegPixEtaHist")

effHist1.SetDirectory(0)
effHist2.SetDirectory(0)
effEtaHist1.SetDirectory(0)
effEtaHist2.SetDirectory(0)
effPhiHist1.SetDirectory(0)
effPhiHist2.SetDirectory(0)
effLen2DHist1.SetDirectory(0)
effLen2DHist2.SetDirectory(0)
effLen3DHist1.SetDirectory(0)
effLen3DHist2.SetDirectory(0)

fakeRateHist1.SetDirectory(0)
fakeRateHist2.SetDirectory(0)

#effLegPixHist1.SetDirectory(0)
#effMyPixHist2.SetDirectory(0)
#effRegPixHist2.SetDirectory(0)

#effLegPixEtaHist1.SetDirectory(0)
#effMyPixEtaHist2.SetDirectory(0)
#effRegPixEtaHist2.SetDirectory(0)

inFile1.Close()
inFile2.Close()

#Drawing options

lengthBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3])
length3DBinning = np.array([1e-2, 5e-2, 1e-1, 3e-1, 9e-1, 3, 20])

emptyHist = ROOT.TH1D("emptyHist", "title;p_{T} [GeV/c];Efficiency", 10, 0, 30)
emptyEtaHist = ROOT.TH1D("emptyEtaHist", "title;#eta;Efficiency", 10, -2.5, 2.5)
emptyPhiHist = ROOT.TH1D("emptyPhiHist", "title;#phi;Efficiency", 8, -3.2, +3.2)
emptyLen2DHist = ROOT.TH1D("emptyLen2DHist", "title;Bs 2D flight length [cm];Efficiency", len(lengthBinning)-1, lengthBinning)
emptyLen3DHist = ROOT.TH1D("emptyLen3DHist", "title;Bs 3D flight length [cm];Efficiency", len(length3DBinning)-1, length3DBinning)
emptyRateHist = ROOT.TH1D("emptyRateHist", "title;p_{T} [GeV/c];Fake rate", 10, 0, 30)
emptyIters = ROOT.TH1D("emptyIters", "title;;Efficiency", 6, 0, 6)

emptyHist.GetYaxis().SetRangeUser(0., 1.05)
emptyEtaHist.GetYaxis().SetRangeUser(0, 1.05)
emptyPhiHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen2DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyLen3DHist.GetYaxis().SetRangeUser(0, 1.05)
emptyRateHist.GetYaxis().SetRangeUser(0, 0.6)
emptyIters.GetYaxis().SetRangeUser(0, 1.05)

emptyHist.SetStats(ROOT.kFALSE)
emptyEtaHist.SetStats(ROOT.kFALSE)
emptyPhiHist.SetStats(ROOT.kFALSE)
emptyLen2DHist.SetStats(ROOT.kFALSE)
emptyLen3DHist.SetStats(ROOT.kFALSE)
emptyRateHist.SetStats(ROOT.kFALSE)
emptyIters.SetStats(ROOT.kFALSE)

idx = 0
binLabels = ["pixeltracks", "reg. pixelTracks", "Iter0", "Iter1", "Iter2", "After filter"]
for label in binLabels:
  emptyIters.GetXaxis().SetBinLabel(idx+1, label)
  idx += 1
emptyIters.SetLabelSize(0.04)

line = ROOT.TLine(1.2, 0, 1.2, 1.05)
line.SetLineColor(ROOT.kBlue)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)

line2 = ROOT.TLine(1.2, 0, 1.2, 0.6)
line2.SetLineColor(ROOT.kBlue)
line2.SetLineWidth(2)
line2.SetLineStyle(ROOT.kDashed)

effHist1.SetLineColor(ROOT.kRed)
effHist2.SetLineColor(ROOT.kBlack)
effEtaHist1.SetLineColor(ROOT.kRed)
effEtaHist2.SetLineColor(ROOT.kBlack)
effPhiHist1.SetLineColor(ROOT.kRed)
effPhiHist2.SetLineColor(ROOT.kBlack)
effLen2DHist1.SetLineColor(ROOT.kRed)
effLen2DHist2.SetLineColor(ROOT.kBlack)
effLen3DHist1.SetLineColor(ROOT.kRed)
effLen3DHist2.SetLineColor(ROOT.kBlack)

fakeRateHist1.SetLineColor(ROOT.kRed)
fakeRateHist2.SetLineColor(ROOT.kBlack)

effOverIters1.SetLineColor(ROOT.kRed)
effOverIters2.SetLineColor(ROOT.kBlack)

#effLegPixHist1.SetLineColor(ROOT.kRed)
#effMyPixHist2.SetLineColor(ROOT.kBlack)
#effRegPixHist2.SetLineColor(ROOT.kBlack)

#effLegPixEtaHist1.SetLineColor(ROOT.kRed)
#effMyPixEtaHist2.SetLineColor(ROOT.kBlack)
#effRegPixEtaHist2.SetLineColor(ROOT.kBlack)

effHist1.SetLineWidth(2)
effHist2.SetLineWidth(2)
effEtaHist1.SetLineWidth(2)
effEtaHist2.SetLineWidth(2)
effPhiHist1.SetLineWidth(2)
effPhiHist2.SetLineWidth(2)
effLen2DHist1.SetLineWidth(2)
effLen2DHist2.SetLineWidth(2)
effLen3DHist1.SetLineWidth(2)
effLen3DHist2.SetLineWidth(2)

fakeRateHist1.SetLineWidth(2)
fakeRateHist2.SetLineWidth(2)

effOverIters1.SetLineWidth(2)
effOverIters2.SetLineWidth(2)

#effLegPixHist1.SetLineWidth(2)
#effMyPixHist2.SetLineWidth(2)
#effRegPixHist2.SetLineWidth(2)

#effLegPixEtaHist1.SetLineWidth(2)
#effMyPixEtaHist2.SetLineWidth(2)
#effRegPixEtaHist2.SetLineWidth(2)

#Legends

effLegend = ROOT.TLegend(0.55, 0.25, 0.85, 0.35)
effLegend.SetTextSize(0.03)
effLegend.AddEntry(effHist1, text1)
effLegend.AddEntry(effHist2, text2) 

effEtaLegend = ROOT.TLegend(0.65, 0.75, 0.95, 0.85)
effEtaLegend.SetTextSize(0.03)
effEtaLegend.AddEntry(effEtaHist1, text1)
effEtaLegend.AddEntry(effEtaHist2, text2)

#Define canvas and print

canvas = ROOT.TCanvas("canvas")
canvas.cd()

canvas.Print(outFileName+"[")

emptyHist.SetTitle(effHist1.GetTitle())
emptyHist.Draw()
effHist1.Draw("same")
effHist2.Draw("same")
line.Draw("same")
effLegend.Draw("same")
canvas.Print(outFileName)

emptyEtaHist.SetTitle(effEtaHist1.GetTitle())
emptyEtaHist.Draw()
effEtaHist1.Draw("same")
effEtaHist2.Draw("same")
effEtaLegend.Draw("same")
canvas.Print(outFileName)

emptyPhiHist.SetTitle(effPhiHist1.GetTitle())
emptyPhiHist.Draw()
effPhiHist1.Draw("same")
effPhiHist2.Draw("same")
effLegend.Draw("same")
canvas.Print(outFileName)

canvas.cd(0).SetLogx(1)
emptyLen2DHist.SetTitle(effLen2DHist1.GetTitle())
emptyLen2DHist.Draw()
effLen2DHist1.Draw("same")
effLen2DHist2.Draw("same")
effLegend.Draw("same")
canvas.Print(outFileName)

emptyLen3DHist.SetTitle(effHist1.GetTitle())
emptyLen3DHist.Draw()
effLen3DHist1.Draw("same")
effLen3DHist2.Draw("same")
effLegend.Draw("same")
canvas.Print(outFileName)

canvas.cd(0).SetLogx(0)
emptyRateHist.SetTitle(fakeRateHist1.GetTitle())
emptyRateHist.Draw()
fakeRateHist1.Draw("same")
fakeRateHist2.Draw("same")
effLegend.Draw("same")
line2.Draw("same")
canvas.Print(outFileName)

emptyIters.SetTitle("Tracking efficiency at HLT steps")
emptyIters.Draw("")
effOverIters1.Draw("pesame")
effOverIters2.Draw("pesame")
effLegend.Draw("same")
canvas.Print(outFileName) 

canvas.Print(outFileName+"]")


