Working code in the "buildBsCandidateFromGEN" branch.

This repo is meant to set up some code to run efficiency/fake rate measurements
for BPH HLT trigger paths.

Set up working area:
```sh
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_11_2_1_Patatrack
cd CMSSW_11_2_1_Patatrack/src
cmsenv
```

Clone repository and compile:
```sh
mkdir TriggerPerformance
cd TriggerPerformance
git clone https://github.com/vjmastra/VanillaHLTAnalyzer.git
git checkout buildBsCandidatesFromGEN
cd ..
scram b -j4
```

A test configuration file is available (./python/hltGPU_JpsiTrkMenu.py).
Input files (with corresponding global tag) and other options can be changed in the file.
To run the menu:
```sh
cd python
nohup taskset -c 0-3 cmsRun hltGPU_JpsiTrkMenu.py &> log &
```

The menu calls the analyzer (./plugins/VanillaHLTAnalyzer.\*) to output root files,
which can be analyzed with some scripts (included in the repo)
