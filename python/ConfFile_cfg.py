import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/a/aboletti/public/TriggerPerfPlots/Charmonium-2017C-AOD.root'
    )
)

process.demo = cms.EDAnalyzer('VanillaHLTAnalyzer',

                              triggerResult = cms.untracked.InputTag("TriggerResults::HLT"),
                              triggerSummary = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                              L1Candidates = cms.untracked.InputTag("gmtStage2Digis:Muon:RECO"),

                              offlineVtx = cms.InputTag("offlinePrimaryVertices::RECO"),
                              beamspot = cms.InputTag("offlineBeamSpot::RECO"),

                              OfflineMuonsTag = cms.untracked.InputTag("muons::RECO"),
                              OfflineTkTag = cms.untracked.InputTag("generalTracks::RECO")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("testNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )


process.p = cms.Path(process.demo)
