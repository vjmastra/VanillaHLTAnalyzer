import sys
import json
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np

def sequenceFromLabel(label=""):
  """
  Get the logical block (sequence) for a module
  from its name (label)
  """
  sequence = ""
  #modules for legacy trigger paths
  if ("hltPixelTracks" in label and "ForSeedsJpsi" in label and "GPU" not in label):
    sequence = "legacyPixelTracksBuildingSequence"
  elif (label == "hltPixelLayerQuadruplets"):
    sequence = "legacyPixelTracksBuildingSequence"
  elif ("hltIter0DisplacedJpsi" in label and "GPU" not in label):
    sequence = "legacyIter0TrackingSequence"
  elif ("hltIter1DisplacedJpsi" in label and "GPU" not in label):
    sequence = "legacyIter1TrackingSequence"
  elif ("hltIter2DisplacedJpsi" in label) and "GPU" not in label:
    sequence = "legacyIter2TrackingSequence"
  elif ("hltTripletRecovery" in label and "GPU" not in label):
    sequence = "legacyIterTripletRecoverySequence"
  elif ("hltDoubletRecovery" in label and "GPU" not in label):
    sequence = "legacyIterDoubletRecoverySequence"
  elif ("hltMergedTracksDisplacedJpsiReg" in label and "GPU" not in label):
    sequence = "legacyIterDoubletRecoverySequence"
  elif ("hltJpsiTk" in label and "TkTk" not in label and "GPU" not in label):
    sequence = "legacyJpsiTrkFilter"
  #modules for upgraded trigger paths
  elif ("hltOnlineBeamSpotToCUDA" in label):
    sequence = "newPixelTracksBuildingSequence"
  elif ("hltSiPixel" in label):
    sequence = "newPixelTracksBuildingSequence" 
  elif ("hltPixelTracks" in label and "Jpsi" not in label):
    sequence = "newPixelTracksBuildingSequence"
  elif ("PixelVertices" in label):
    sequence = "newPixelTracksBuildingSequence"
  elif ("hltPixelTracksTrackingRegionsForSeedsJpsiGPU" == label):
    sequence = "newSeedsFromRegPixelTracks"
  elif ("hltPixelTracksDisplacedJpsiRegional" == label):
    sequence = "newSeedsFromRegPixelTracks"
  elif ("hltIter0DisplacedJpsiPixelSeedsFromPixelTracksGPU" == label):
    sequence = "newSeedsFromRegPixelTracks"
  elif ("hltIter0DisplacedJpsi" in label and "GPU" in label):
    sequence = "newPixelTracksSequenceIter0"
  elif ("hltIter1DisplacedJpsi" in label and "GPU" in label):
    sequence = "newPixelTracksSequenceIter1"
  elif ("hltIter2DisplacedJpsi" in label and "GPU" in label):
    sequence = "newPixelTracksSequenceIter2"
  elif ("hltJpsiTk" in label and "TkTk" not in label and "GPU" in label):
    sequence = "newJpsiTrkFilter"
  elif ("Muon" in label or "Displacedmumu" in label):
    sequence = "JpsiTriggerPath"
  else:
    sequence = "other"

  return sequence

#Read from json

if len(sys.argv) != 2:
  print "USAGE: %s <input file>"%(sys.argv[0])
  sys.exit(1)

filename = sys.argv[1]
with open(filename) as f:
  data = json.load(f)

df = DataFrame(data['modules'])

#Add useful quantities to DataFrame

df['sequence'] = df['label'].apply(sequenceFromLabel)
df['path'] = df['sequence'].apply(lambda seq: "Legacy" if ("legacy" in seq) else ("New" if ("new" in seq) else ("Jpsi" if ("Jpsi" in seq) else "other")))
df['type'] = df['sequence'].apply(lambda seq: "Building" if ("Building" in seq) else ("Tracking" if ("Iter" in seq or "RegPixelTracks" in seq) else ("Filter" if ("Filter" in seq) else "Other")))
df['time/event'] = df.apply(lambda x: 0 if x['events'] == 0 else x['time_real']/x['events'], axis=1)

newPatatrackBuildingSequenceAveTime = df[df.sequence == "newPixelTracksBuildingSequence"]['time/event'].sum()
newRegionalTrackingSequenceAveTime  = df[df.sequence == "newSeedsFromRegPixelTracks"]['time/event'].sum()

#Variables for plots

pathLabels = ["Jpsi", "Legacy", "New"]
pathValues = [df[df.path==label]['time/event'].sum() for label in pathLabels]

newPatatrackBuildingSequenceAveTime = df[df.sequence == "newPixelTracksBuildingSequence"]['time/event'].sum()

typeLabels = ["Building", "Tracking", "Filter"]
legacyValues = [df[(df['path']=="Legacy") & (df['type']==label)]['time/event'].sum() for label in typeLabels]
newValues = [df[(df['path']=="New") & (df['type']==label)]['time/event'].sum() for label in typeLabels]

newRegionalTrackingSequenceAveTime  = df[df.sequence == "newSeedsFromRegPixelTracks"]['time/event'].sum()

legacySeqLabels = ["legacyPixelTracksBuildingSequence", "legacyIter0TrackingSequence", "legacyIter1TrackingSequence", "legacyIter2TrackingSequence", "legacyIterTripletRecoverySequence", "legacyIterDoubletRecoverySequence", "legacyJpsiTrkFilter"]
legacySeqValues = [df[df['sequence']==seq]['time/event'].sum() for seq in legacySeqLabels]

newSeqLabels = ["newPixelTracksBuildingSequence", "newSeedsFromRegPixelTracks", "newPixelTracksSequenceIter0", "newPixelTracksSequenceIter1", "newPixelTracksSequenceIter2", "newJpsiTrkFilter"]
newSeqValues = [df[df['sequence']==seq]['time/event'].sum() for seq in newSeqLabels]

newPatatrackBuildingSequenceAveTime = df[df.sequence == "newPixelTracksBuildingSequence"]['time/event'].sum()
newRegionalTrackingSequenceAveTime  = df[df.sequence == "newSeedsFromRegPixelTracks"]['time/event'].sum()

#Plots

outputFile1 = "./pathTiming.pdf"
outputFile2 = "./pathTimingStepsComparison.pdf"
outputFile3 = "./legacySequenceTimingBySteps.pdf"
outputFile4 = "./newSequenceTimingBySteps.pdf"

x = np.arange(len(pathLabels))
plt.ylim(1, 1000)
plt.bar(x, pathValues, width=0.3, log=True, label="Tracking sequence", color="orange")
plt.bar(x[-1], newPatatrackBuildingSequenceAveTime, width=0.3, alpha=0.5, hatch="//", color="red", log=True, label="Patatrack pixeltracks' building")
plt.xticks(x, pathLabels)
plt.title("Time spent in HLT paths (avg)")
plt.xlabel("")
plt.ylabel("Time (ms)")
plt.legend()
plt.savefig(outputFile1)
print("printed:\t"+outputFile1)
plt.clf()
plt.cla()

x = np.arange(len(typeLabels))
plt.ylim(0.1, 1000)
plt.bar(x-0.15, legacyValues, width=0.3, log=True, label="Legacy tracking sequence", color="green")
plt.bar(x+0.15, newValues, width=0.3, log=True, label="New tracking sequence", color="orange")
plt.bar(x[0]+0.15, newPatatrackBuildingSequenceAveTime, width=0.3, log=True, alpha=0.5, hatch="//", color="red", label="Patatrack pixeltracks' building", bottom=newRegionalTrackingSequenceAveTime)
plt.xticks(x, typeLabels)
plt.title("Time spent in HLT sequences (avg)")
plt.xlabel("")
plt.ylabel("Time (ms)")
plt.legend()
plt.savefig(outputFile2)
print("printed:\t"+outputFile2)
plt.clf()
plt.cla()

legacySeqLabels2 = ["Building", "Iter0", "Iter1", "Iter2", "3let-Rec.", "2let-Rec", "Filter"]
x = np.arange(len(legacySeqLabels))
plt.ylim(0.1, 100)
plt.bar(x, legacySeqValues, width=0.3, log=True, label="Legacy tracking sequence", color="green")
plt.xticks(x, legacySeqLabels2)
plt.title("Time spent in HLT sequences step-by-step (avg)")
plt.xlabel("")
plt.ylabel("Time (ms)")
plt.legend()
plt.savefig(outputFile3)
print("printed:\t"+outputFile3)
plt.clf()
plt.cla()

newSeqLabels2 = ["Building", "Reg. selection", "Iter0", "Iter1", "Iter2", "JpsiTrk filter"]
x = np.arange(len(newSeqLabels))
plt.ylim(0.1, 100)
plt.bar(x, newSeqValues, width=0.3, log=True, label="New tracking sequence", color="orange")
plt.xticks(x, newSeqLabels2)
plt.title("Time spent in HLT sequences step-by-step (avg)")
plt.xlabel("")
plt.ylabel("Time (ms)")
plt.legend()
plt.savefig(outputFile4)
print("printed:\t"+outputFile4)
plt.clf()
plt.cla()
