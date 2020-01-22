Set up working area:
```sh
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_15_patch2
cd CMSSW_10_2_15_patch2/src
cmsenv
```

Clone repository and compile:
```sh
mkdir TriggerPerformance
git clone git@github.com:aboletti/VanillaHLTAnalyzer.git
cd ..
scram b -j2
```

Copy configuration file:
```sh
wget https://gist.githubusercontent.com/aboletti/947b4bd0e9629ce6e2d7ad55811be5e7/raw/RunNtuple_miniAOD.py
```
and run on crab or local file.

Read the output files and produce the plots with [this macro](https://gist.github.com/aboletti/a0776d7e35e444b3c452933de52ca39e), adjusting the paths for the input files.
