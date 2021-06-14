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
git clone git@github.com:vjmastra/VanillaHLTAnalyzer.git
git checkout buildBsCandidatesFromGEN
cd ..
scram b -j4
```
-------------------------------------------------------
Copy configuration file:
```sh
wget https://gist.githubusercontent.com/aboletti/947b4bd0e9629ce6e2d7ad55811be5e7/raw/RunNtuple_miniAOD.py
```
and run on crab or local file.

Read the output files and produce the plots with [this macro](https://gist.github.com/aboletti/a0776d7e35e444b3c452933de52ca39e), adjusting the paths for the input files.
