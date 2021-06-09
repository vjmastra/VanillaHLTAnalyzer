// -*- C++ -*-
//
// Package:    TriggerPerformance/VanillaHLTAnalyzer
// Class:      VanillaHLTAnalyzer
// 
/**\class VanillaHLTAnalyzer VanillaHLTAnalyzer.cc TriggerPerformance/VanillaHLTAnalyzer/plugins/VanillaHLTAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alessio Boletti
//         Created:  Wed, 18 Apr 2018 11:04:38 GMT
//
//

#include "VanillaHLTAnalyzer.h"

//
// constructors and destructor
//
VanillaHLTAnalyzer::VanillaHLTAnalyzer(const edm::ParameterSet& iConfig):

  triggerResultTag_       (iConfig.getUntrackedParameter<edm::InputTag>("triggerResult")), 
  triggerResultToken_     (consumes<edm::TriggerResults>(triggerResultTag_)),

  genParticlesTag_(iConfig.getUntrackedParameter<edm::InputTag>("genParticles")),
  genParticlesToken_(consumes< std::vector<reco::GenParticle> >(genParticlesTag_)),

  triggerObjectsTag_      (iConfig.getUntrackedParameter<edm::InputTag>("triggerObjects")),
  triggerObjectsToken_    (consumes<trigger::TriggerEvent>(triggerObjectsTag_)),

  jpsiVerticesTag_        (iConfig.getUntrackedParameter<edm::InputTag>("jpsiVertices")),
  jpsiVerticesToken_      (consumes<vector<reco::Vertex>>(jpsiVerticesTag_)),

  legacyPixelTracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacyPixelTracks")),
  legacyPixelTracksToken_      (consumes<vector<reco::Track>>(legacyPixelTracksTag_)),
  
  legacyIter0TracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacyIter0Tracks")),
  legacyIter0TracksToken_      (consumes<vector<reco::Track>>(legacyIter0TracksTag_)),

  legacyIter1TracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacyIter1Tracks")),
  legacyIter1TracksToken_      (consumes<vector<reco::Track>>(legacyIter1TracksTag_)),

  legacyIter2TracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacyIter2Tracks")),
  legacyIter2TracksToken_      (consumes<vector<reco::Track>>(legacyIter2TracksTag_)),

  legacy3RecoTracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacy3RecoTracks")),
  legacy3RecoTracksToken_      (consumes<vector<reco::Track>>(legacy3RecoTracksTag_)),

  legacy2RecoTracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("legacy2RecoTracks")),
  legacy2RecoTracksToken_      (consumes<vector<reco::Track>>(legacy2RecoTracksTag_)),

  myPixelTracksTag_        (iConfig.getUntrackedParameter<edm::InputTag>("myPixelTracks")),
  myPixelTracksToken_      (consumes<vector<reco::Track>>(myPixelTracksTag_)),

  regionalPixelTracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("regionalPixelTracks")),
  regionalPixelTracksToken_ (consumes<vector<reco::Track>>(regionalPixelTracksTag_)),

  ctfIter0TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("ctfIter0Tracks")),
  ctfIter0TracksToken_ (consumes<vector<reco::Track>>(ctfIter0TracksTag_)),

  triggerIter0TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("triggerIter0Tracks")),
  triggerIter0TracksToken_ (consumes<vector<reco::Track>>(triggerIter0TracksTag_)),

  ctfIter1TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("ctfIter1Tracks")),
  ctfIter1TracksToken_ (consumes<vector<reco::Track>>(ctfIter1TracksTag_)),

  triggerIter1TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("triggerIter1Tracks")),
  triggerIter1TracksToken_ (consumes<vector<reco::Track>>(triggerIter1TracksTag_)),

  ctfIter2TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("ctfIter2Tracks")),
  ctfIter2TracksToken_ (consumes<vector<reco::Track>>(ctfIter2TracksTag_)),

  triggerIter2TracksTag_   (iConfig.getUntrackedParameter<edm::InputTag>("triggerIter2Tracks")),
  triggerIter2TracksToken_ (consumes<vector<reco::Track>>(triggerIter2TracksTag_)),

  beamSpotTag_            (iConfig.getUntrackedParameter<edm::InputTag>("beamspot")),
  beamSpotToken_          (consumes<reco::BeamSpot>(beamSpotTag_)),

  hltTag_                 (iConfig.getUntrackedParameter<std::string>("hltTag"))
  // hltPrescale_ (iConfig, consumesCollector(), *this)
  // hltPrescale_ (new HLTPrescaleProvider(iConfig, consumesCollector(), *this))

{
  //now do what ever initialization is needed
  // usesResource("TFileService");
  // fGtUtil = new l1t::L1TGlobalUtil(iConfig, consumesCollector(), *this, iConfig.getParameter<edm::InputTag>("l1results"), iConfig.getParameter<edm::InputTag>("l1results"));

}


VanillaHLTAnalyzer::~VanillaHLTAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  // delete hltPrescale_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VanillaHLTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   hltNames.clear();
   hltResults.clear();

   Muon1_pT = 0;
   Muon1_eta = 0;
   Muon1_phi = 0;
   Muon1_cha = 0;
   Muon2_pT = 0;
   Muon2_eta = 0;
   Muon2_phi = 0;
   Muon2_cha = 0;

   track1_pT  = 0;
   track1_eta = 0;
   track1_phi = 0;
   track1_cha = 0;
   track2_pT  = 0;
   track2_eta = 0;
   track2_phi = 0;
   track2_cha = 0;

   bsFlightLength2D = 0;
   bsFlightLength3D = 0;

   passHLT_DoubleMu4_Jpsi_Displaced = 0;
   passHLT_DoubleMu4_JpsiTrk_Displaced = 0;
   passHLT_DoubleMu4_JpsiTrk_Displaced_GPU = 0;

   muonsMatched = 0;

   hltMatch_Track1 = false;
   hltMatch_Track2 = false;
   gpuMatch_Track1 = false;
   gpuMatch_Track2 = false;

   hltTriggerTrack_pt.clear();
   hltTriggerTrack_eta.clear();
   hltTriggerFake_pt.clear();
   hltTriggerFake_eta.clear();
   gpuTriggerTrack_pt.clear();
   gpuTriggerTrack_eta.clear();
   gpuTriggerFake_pt.clear();
   gpuTriggerFake_eta.clear();

   legPixMatch_Track1 = false;
   legPixMatch_Track2 = false;
   legI0Match_Track1 = false;
   legI0Match_Track2 = false;
   legI1Match_Track1 = false;
   legI1Match_Track2 = false;
   legI2Match_Track1 = false;
   legI2Match_Track2 = false;
   leg3RMatch_Track1 = false;
   leg3RMatch_Track2 = false;
   leg2RMatch_Track1 = false;
   leg2RMatch_Track2 = false;

   myPixMatch_Track1 = false;
   myPixMatch_Track2 = false;
   regMatch_Track1 = false;
   regMatch_Track2 = false;
   ctf0Match_Track1 = false;
   ctf0Match_Track2 = false;
   iter0Match_Track1 = false;
   iter0Match_Track2 = false;
   ctf1Match_Track1 = false;
   ctf1Match_Track2 = false;
   iter1Match_Track1 = false;
   iter1Match_Track2 = false;
   ctf2Match_Track1 = false;
   ctf2Match_Track2 = false;
   iter2Match_Track1 = false;
   iter2Match_Track2 = false;

   hltMatch_Track = 0;
   gpuMatch_Track = 0;

   j_vtx_x = -999;
   j_vtx_y = -999;
   j_vtx_z = -999;

   j_vtx_size = 0;   

   hlt_Mu1_pT = -10;
   hlt_Mu1_eta = -10;
   hlt_Mu1_phi = -10;
   hlt_Mu1_cha = 0; //not handled
   hlt_Mu2_pT = -10;
   hlt_Mu2_eta = -10;
   hlt_Mu2_phi = -10;
   hlt_Mu2_cha = 0; //not handled

   hlt_Track1_pT = -10;
   hlt_Track1_eta = -10;
   hlt_Track1_phi = -10;
   hlt_Track2_pT = -10;
   hlt_Track2_eta = -10;
   hlt_Track2_phi = -10;

   hlt_Track1_deltaR = 999;
   hlt_Track2_deltaR = 999; 

   gpu_Track1_pT = -10;
   gpu_Track1_eta = -10;
   gpu_Track1_phi = -10;
   gpu_Track2_pT = -10;
   gpu_Track2_eta = -10;
   gpu_Track2_phi = -10;
 
   gpu_Track1_deltaR = 999;
   gpu_Track2_deltaR = 999;

   Mu1_deltaR = -10;
   Mu2_deltaR = -10;

   legPix_Track1_pT = -10;
   legPix_Track1_eta = -10;
   legPix_Track1_phi = -10;
   legPix_Track1_d0 = -10;
   legPix_Track2_pT = -10;
   legPix_Track2_eta = -10;
   legPix_Track2_phi = -10;
   legPix_Track2_d0 = -10;

   legPix_Track1_deltaR = 999;
   legPix_Track2_deltaR = 999;

   legPix_Track1_dxy  = -1.;
   legPix_Track1_dz   = -1.;
   legPix_Track1_dEtaMu1 = -1.;
   legPix_Track1_dPhiMu1 = -1.;
   legPix_Track1_dEtaMu2 = -1.;
   legPix_Track1_dPhiMu2 = -1.;
   legPix_Track2_dxy  = -1.;
   legPix_Track2_dz   = -1.;
   legPix_Track2_dEtaMu1 = -1.;
   legPix_Track2_dPhiMu1 = -1.;
   legPix_Track2_dEtaMu2 = -1.;
   legPix_Track2_dPhiMu2 = -1.;

   legPix_Track1_nHits = 0;
   legPix_Track1_normChi2 = 0.;
   legPix_Track2_nHits = 0;
   legPix_Track2_normChi2 = 0.;

   legPix_collSize = 0;

   legI0_Track1_pT = -10;
   legI0_Track1_eta = -10;
   legI0_Track1_phi = -10;
   legI0_Track1_d0 = -10;
   legI0_Track2_pT = -10;
   legI0_Track2_eta = -10;
   legI0_Track2_phi = -10;
   legI0_Track2_d0 = -10;

   legI0_collSize = 0;

   legI1_Track1_pT = -10;
   legI1_Track1_eta = -10;
   legI1_Track1_phi = -10;
   legI1_Track1_d0 = -10;
   legI1_Track2_pT = -10;
   legI1_Track2_eta = -10;
   legI1_Track2_phi = -10;
   legI1_Track2_d0 = -10;

   legI1_collSize = 0;

   legI2_Track1_pT = -10;
   legI2_Track1_eta = -10;
   legI2_Track1_phi = -10;
   legI2_Track1_d0 = -10;
   legI2_Track2_pT = -10;
   legI2_Track2_eta = -10;
   legI2_Track2_phi = -10;
   legI2_Track2_d0 = -10;

   legI2_collSize = 0;

   leg3R_Track1_pT = -10;
   leg3R_Track1_eta = -10;
   leg3R_Track1_phi = -10;
   leg3R_Track1_d0 = -10;
   leg3R_Track2_pT = -10;
   leg3R_Track2_eta = -10;
   leg3R_Track2_phi = -10;
   leg3R_Track2_d0 = -10;

   leg2R_Track1_pT = -10;
   leg2R_Track1_eta = -10;
   leg2R_Track1_phi = -10;
   leg2R_Track1_d0 = -10;
   leg2R_Track1_d0err = -10;
   leg2R_Track2_pT = -10;
   leg2R_Track2_eta = -10;
   leg2R_Track2_phi = -10;
   leg2R_Track2_d0 = -10;
   leg2R_Track2_d0err = -10;

   prefilter_legTks_pT.clear();
   prefilter_legTks_eta.clear();
   prefilter_legTks_phi.clear();
   prefilter_legTks_dEta1.clear();
   prefilter_legTks_dPhi1.clear();
   prefilter_legTks_dEta2.clear();
   prefilter_legTks_dPhi2.clear();
   prefilter_legTks_dxy.clear();
   prefilter_legTks_dz.clear();
   prefilter_legTks_d0.clear();
   prefilter_legTks_d0err.clear();
   prefilter_legTks_bool.clear();

   myPix_Track1_pT = -10;
   myPix_Track1_eta = -10;
   myPix_Track1_phi = -10;
   myPix_Track1_d0 = -10;
   myPix_Track2_pT = -10;
   myPix_Track2_eta = -10;
   myPix_Track2_phi = -10;
   myPix_Track2_d0 = -10;

   myPix_Track1_deltaR = 999;
   myPix_Track2_deltaR = 999;

   myPix_Track1_dxy  = -1.;
   myPix_Track1_dz   = -1.;
   myPix_Track1_dEtaMu1 = -1.;
   myPix_Track1_dPhiMu1 = -1.;
   myPix_Track1_dEtaMu2 = -1.;
   myPix_Track1_dPhiMu2 = -1.;
   myPix_Track2_dxy  = -1.;
   myPix_Track2_dz   = -1.;
   myPix_Track2_dEtaMu1 = -1.;
   myPix_Track2_dPhiMu1 = -1.;
   myPix_Track2_dEtaMu2 = -1.;
   myPix_Track2_dPhiMu2 = -1.;

   myPix_Track1_nHits = 0;
   myPix_Track1_normChi2 = 0.;
   myPix_Track2_nHits = 0;
   myPix_Track2_normChi2 = 0.;

   myPix_collSize = 0;

   reg_Track1_pT = -10;
   reg_Track1_eta = -10;
   reg_Track1_phi = -10;
   reg_Track1_d0 = -10;
   reg_Track2_pT = -10;
   reg_Track2_eta = -10;
   reg_Track2_phi = -10;
   reg_Track2_d0 = -10; 

   reg_Track1_deltaR = 999;
   reg_Track2_deltaR = 999;

   reg_Track1_dxy  = -1.;
   reg_Track1_dz   = -1.;
   reg_Track1_dEtaMu1 = -1.;
   reg_Track1_dPhiMu1 = -1.;
   reg_Track1_dEtaMu2 = -1.;
   reg_Track1_dPhiMu2 = -1.;
   reg_Track2_dxy  = -1.;
   reg_Track2_dz   = -1.;
   reg_Track2_dEtaMu1 = -1.;
   reg_Track2_dPhiMu1 = -1.;
   reg_Track2_dEtaMu2 = -1.;
   reg_Track2_dPhiMu2 = -1.;

   reg_Track1_nHits = 0;
   reg_Track1_normChi2 = 0.;
   reg_Track2_nHits = 0;
   reg_Track2_normChi2 = 0.;

   reg_collSize = 0;
   regionalFrac_collSize = 0;

   ctf0_Track1_pT = -10;
   ctf0_Track1_eta = -10;
   ctf0_Track1_phi = -10;
   ctf0_Track1_d0 = -10;
   ctf0_Track2_pT = -10;
   ctf0_Track2_eta = -10;
   ctf0_Track2_phi = -10;
   ctf0_Track2_d0 = -10;

   ctf0_Track1_deltaR = 999;
   ctf0_Track2_deltaR = 999;
 
   ctf0_Track1_dxy  = -1.;
   ctf0_Track1_dz   = -1.;
   ctf0_Track1_dEtaMu1 = -1.;
   ctf0_Track1_dPhiMu1 = -1.;
   ctf0_Track1_dEtaMu2 = -1.;
   ctf0_Track1_dPhiMu2 = -1.;
   ctf0_Track2_dxy  = -1.;
   ctf0_Track2_dz   = -1.;
   ctf0_Track2_dEtaMu1 = -1.;
   ctf0_Track2_dPhiMu1 = -1.;
   ctf0_Track2_dEtaMu2 = -1.;
   ctf0_Track2_dPhiMu2 = -1.;

   ctf0_Track1_nHits = 0;
   ctf0_Track1_normChi2 = 0.;
   ctf0_Track2_nHits = 0;
   ctf0_Track2_normChi2 = 0.;

   ctf0_collSize = 0;

   iter0_Track1_pT = -10;
   iter0_Track1_eta = -10;
   iter0_Track1_phi = -10;
   iter0_Track1_d0 = -10;
   iter0_Track2_pT = -10;
   iter0_Track2_eta = -10;
   iter0_Track2_phi = -10;
   iter0_Track2_d0 = -10;
 
   iter0_collSize = 0;

   ctf1_Track1_pT = -10;
   ctf1_Track1_eta = -10;
   ctf1_Track1_phi = -10;
   ctf1_Track1_d0 = -10;
   ctf1_Track2_pT = -10;
   ctf1_Track2_eta = -10;
   ctf1_Track2_phi = -10;
   ctf1_Track2_d0 = -10; 

   ctf1_Track1_deltaR = 999;
   ctf1_Track2_deltaR = 999;

   ctf1_Track1_dxy  = -1.;
   ctf1_Track1_dz   = -1.;
   ctf1_Track1_dEtaMu1 = -1.;
   ctf1_Track1_dPhiMu1 = -1.;
   ctf1_Track1_dEtaMu2 = -1.;
   ctf1_Track1_dPhiMu2 = -1.;
   ctf1_Track2_dxy  = -1.;
   ctf1_Track2_dz   = -1.;
   ctf1_Track2_dEtaMu1 = -1.;
   ctf1_Track2_dPhiMu1 = -1.;
   ctf1_Track2_dEtaMu2 = -1.;
   ctf1_Track2_dPhiMu2 = -1.;

   ctf1_Track1_nHits = 0;
   ctf1_Track1_normChi2 = 0.;
   ctf1_Track2_nHits = 0;
   ctf1_Track2_normChi2 = 0.;

   ctf1_collSize = 0;

   iter1_Track1_pT = -10;
   iter1_Track1_eta = -10;
   iter1_Track1_phi = -10;
   iter1_Track1_d0 = -10;
   iter1_Track2_pT = -10;
   iter1_Track2_eta = -10;
   iter1_Track2_phi = -10;
   iter1_Track2_d0 = -10;

   iter1_collSize = 0;
 
   ctf2_Track1_pT = -10;
   ctf2_Track1_eta = -10;
   ctf2_Track1_phi = -10;
   ctf2_Track1_d0 = -10;
   ctf2_Track2_pT = -10;
   ctf2_Track2_eta = -10;
   ctf2_Track2_phi = -10;
   ctf2_Track2_d0 = -10;
 
   ctf2_Track1_deltaR = 999;
   ctf2_Track2_deltaR = 999;

   ctf2_Track1_dxy  = -1.;
   ctf2_Track1_dz   = -1.;
   ctf2_Track1_dEtaMu1 = -1.;
   ctf2_Track1_dPhiMu1 = -1.;
   ctf2_Track1_dEtaMu2 = -1.;
   ctf2_Track1_dPhiMu2 = -1.;
   ctf2_Track2_dxy  = -1.;
   ctf2_Track2_dz   = -1.;
   ctf2_Track2_dEtaMu1 = -1.;
   ctf2_Track2_dPhiMu1 = -1.;
   ctf2_Track2_dEtaMu2 = -1.;
   ctf2_Track2_dPhiMu2 = -1.;

   ctf2_Track1_nHits = 0;
   ctf2_Track1_normChi2 = 0.;
   ctf2_Track2_nHits = 0;
   ctf2_Track2_normChi2 = 0.;

   ctf2_collSize = 0;

   iter2_Track1_pT = -10;
   iter2_Track1_eta = -10;
   iter2_Track1_phi = -10;
   iter2_Track1_d0 = -10;
   iter2_Track1_d0err = -10;
   iter2_Track2_pT = -10;
   iter2_Track2_eta = -10;
   iter2_Track2_phi = -10;
   iter2_Track2_d0 = -10;
   iter2_Track2_d0err = -10;

   iter2_collSize = 0;

   prefilter_newTks_pT.clear();
   prefilter_newTks_eta.clear();
   prefilter_newTks_phi.clear();
   prefilter_newTks_dEta1.clear();
   prefilter_newTks_dPhi1.clear();
   prefilter_newTks_dEta2.clear();
   prefilter_newTks_dPhi2.clear();
   prefilter_newTks_dxy.clear();
   prefilter_newTks_dz.clear();
   prefilter_newTks_d0.clear();
   prefilter_newTks_d0err.clear();
   prefilter_newTks_bool.clear();

   dummy = 0;

   edm::Handle<edm::TriggerResults> triggerResults;

   if ( !(iEvent.getByToken(triggerResultToken_, triggerResults) ) ) {
     edm::LogError("") << "Trigger collections not found !!!";
     return;
   }

   const edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResults);
   for (unsigned int itrig=0; itrig < triggerNames_.size(); ++itrig) {
     std::string iName = triggerNames_.triggerName(itrig);
     iName.erase( std::remove_if(iName.end()-2, iName.end(), (int(*)(int))std::isdigit), iName.end() );
     if ( iName.find ("Jpsi") != std::string::npos ) {
       hltNames    .push_back(iName);
       hltResults  .push_back(triggerResults->accept(itrig));
     }
   }

   for (unsigned int i = 0; i < hltNames.size(); i++) {
   //  cout << hltNames[i] << endl;
     if (hltNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v") != std::string::npos)
       hltResults[i] == 1 ? passHLT_DoubleMu4_Jpsi_Displaced = 1 : passHLT_DoubleMu4_Jpsi_Displaced = 0;
     if (hltNames[i].find("HLT_DoubleMu4_JpsiTrk_Displaced_v") != std::string::npos)
       hltResults[i] == 1 ? passHLT_DoubleMu4_JpsiTrk_Displaced = 1 : passHLT_DoubleMu4_JpsiTrk_Displaced = 0;
     if (hltNames[i].find("HLT_DoubleMu4_JpsiTrk_Displaced_GPU_v") != std::string::npos)
       hltResults[i] == 1 ? passHLT_DoubleMu4_JpsiTrk_Displaced_GPU = 1 : passHLT_DoubleMu4_JpsiTrk_Displaced_GPU = 0;
   }    

   //PrimaryVertex
   
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByToken(beamSpotToken_, beamSpotHandle);

   if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;
   else {
     cout << "No beam spot available" << endl;
   }

   float x0 = beamSpot.x0();
   float y0 = beamSpot.y0();
   float z0 = beamSpot.z0();

   //Look for gen particles

   edm::Handle< std::vector<reco::GenParticle> >  genParticles;
   iEvent.getByToken(genParticlesToken_, genParticles);

   float x1 = 0.;
   float y1 = 0.;
   float z1 = 0.;
 
   for(std::vector<reco::GenParticle>::const_iterator genBs=genParticles->begin(); genBs!=genParticles->end(); ++genBs) {

     if ( abs(genBs->pdgId()) != 531 ) continue;

//     printProgeny(*genBs);
//     cout<<endl;

     const reco::Candidate* genMu1 = 0;
     const reco::Candidate* genMu2 = 0;
     const reco::Candidate* genKaon1 = 0;
     const reco::Candidate* genKaon2 = 0;

     bool muonPairFound = false;
     bool kaonPairFound = false;

     for (uint iDaug=0; iDaug<genBs->numberOfDaughters(); ++iDaug) {

       if ( genBs->daughter(iDaug)->pdgId() == 443 ) {

	 //cout<<"Jpsi found!"<<endl;
	 const reco::GenParticle* genJpsi = (reco::GenParticle*)genBs->daughter(iDaug);

	 for (uint iGrandDaug1=0; iGrandDaug1<genJpsi->numberOfDaughters(); ++iGrandDaug1) {
	   if ( abs(genJpsi->daughter(iGrandDaug1)->pdgId()) != 13 ) continue;

	   for (uint iGrandDaug2=0; iGrandDaug2<genJpsi->numberOfDaughters(); ++iGrandDaug2) {
	     if ( iGrandDaug2 == iGrandDaug1 ) continue;
	     if ( abs(genJpsi->daughter(iGrandDaug2)->pdgId()) != 13 ) continue;
	     if ( genJpsi->daughter(iGrandDaug2)->pt() > genJpsi->daughter(iGrandDaug1)->pt() ) continue;
	     if ( genJpsi->daughter(iGrandDaug2)->pdgId()*genJpsi->daughter(iGrandDaug1)->pdgId() > 0 ) continue;

	     //cout<<"Muons from Jpsi decay found!"<<endl;
	     muonPairFound = true;
	     genMu1 = genJpsi->daughter(iGrandDaug1);
	     genMu2 = genJpsi->daughter(iGrandDaug2);
	     break;
	   }
           
	   if (muonPairFound) break;
	 }

       }//end Jpsi if

       if ( genBs->daughter(iDaug)->pdgId() == 333 ) {

	 //cout<<"Phi found!"<<endl;
	 const reco::GenParticle* genPhi = (reco::GenParticle*)genBs->daughter(iDaug);

	 for (uint iGrandDaug1=0; iGrandDaug1<genPhi->numberOfDaughters(); ++iGrandDaug1) {
	   if ( abs(genPhi->daughter(iGrandDaug1)->pdgId()) != 321 ) continue;

	   for (uint iGrandDaug2=0; iGrandDaug2<genPhi->numberOfDaughters(); ++iGrandDaug2) {
	     if ( iGrandDaug2 == iGrandDaug1 ) continue;
	     if ( abs(genPhi->daughter(iGrandDaug2)->pdgId()) != 321 ) continue;
	     if ( genPhi->daughter(iGrandDaug2)->pt() > genPhi->daughter(iGrandDaug1)->pt() ) continue;
	     if ( genPhi->daughter(iGrandDaug2)->pdgId()*genPhi->daughter(iGrandDaug1)->pdgId() > 0 ) continue;

	     //cout<<"Kaons from Psi decay found!"<<endl;
	     kaonPairFound = true;
	     genKaon1 = genPhi->daughter(iGrandDaug1);
	     genKaon2 = genPhi->daughter(iGrandDaug2);
	     break;
	   }

	   if (kaonPairFound) break;
	 }

       }//end phi if
	 
     }//end daug for

     if ( !muonPairFound || !kaonPairFound ) continue;

//     cout<<"Bs accepted"<<endl;

     Muon1_pT  = genMu1->pt();
     Muon1_eta = genMu1->eta();
     Muon1_phi = genMu1->phi();
     Muon1_cha = genMu1->charge();
     Muon2_pT  = genMu2->pt();
     Muon2_eta = genMu2->eta();
     Muon2_phi = genMu2->phi();
     Muon2_cha = genMu2->charge();

     track1_pT  = genKaon1->pt();
     track1_eta = genKaon1->eta();
     track1_phi = genKaon1->phi();
     track1_cha = genKaon1->charge();
     track2_pT  = genKaon2->pt();
     track2_eta = genKaon2->eta();
     track2_phi = genKaon2->phi();
     track2_cha = genKaon2->charge();

     x1 = genBs->vx();
     y1 = genBs->vy();
     z1 = genBs->vz();
     x0 = genKaon1->vx();
     y0 = genKaon1->vy();
     z0 = genKaon1->vz();
     bsFlightLength2D = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
     bsFlightLength3D = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

     //Match with trigger object

     edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
     iEvent.getByToken(triggerObjectsToken_, triggerObjectsSummary);

     trigger::TriggerObjectCollection selectedObjects;

     edm::Handle<vector<reco::Vertex>> jpsiVertexCollection;
     iEvent.getByToken(jpsiVerticesToken_, jpsiVertexCollection);
     math::XYZPoint jpsiVertex;     

     bool matchMuon1 = false;
     bool matchMuon2 = false;
     bool matchHLTMuon1 = false;
     bool matchHLTMuon2 = false;
     bool matchTrack = false;
     bool matchTrackToJpsi = false;

     std::vector<float> hltTriggerObjects_pT;
     std::vector<float> hltTriggerObjects_eta;
     std::vector<float> hltTriggerObjects_phi;

     //DoubleMu4_Jpsi_Displaced
 
     if (passHLT_DoubleMu4_Jpsi_Displaced) {
       if(triggerObjectsSummary.isValid()) {
         size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltDisplacedmumuFilterDoubleMu4Jpsi", "", hltTag_) );
         trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
         if (filterIndex < (*triggerObjectsSummary).sizeFilters()) {
           const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
           for (size_t j = 0; j < keys.size(); j++) {
             trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]]; //0 Mu1 1 Mu2 ...
             hltTriggerObjects_pT.push_back(foundObject.pt()); 
             hltTriggerObjects_eta.push_back(foundObject.eta());
             hltTriggerObjects_phi.push_back(foundObject.phi());
           }
         }
       }//end triggerObjects

       float dRthreshold = 0.01;
       j_vtx_size = jpsiVertexCollection->size();
       for (unsigned int i = 0; i < hltTriggerObjects_pT.size(); i=i+2) { //To handle jpsi multiplicity
         float dR1 = deltaR(Muon1_eta, Muon1_phi, hltTriggerObjects_eta[i], hltTriggerObjects_phi[i]);
         float dR2 = deltaR(Muon2_eta, Muon2_phi, hltTriggerObjects_eta[i+1], hltTriggerObjects_phi[i+1]);
         if (dR1 < dRthreshold && dR2 < dRthreshold) {
           hlt_Mu1_pT = hltTriggerObjects_pT[i];
           hlt_Mu1_eta = hltTriggerObjects_eta[i];
           hlt_Mu1_phi = hltTriggerObjects_phi[i];
           Mu1_deltaR = dR1;
           matchMuon1 = true;
           hlt_Mu2_pT = hltTriggerObjects_pT[i+1];
           hlt_Mu2_eta = hltTriggerObjects_eta[i+1];
           hlt_Mu2_phi = hltTriggerObjects_phi[i+1];
           Mu2_deltaR = dR2;
           matchMuon2 = true;
           int jpsiIndex = (int)i/2;
           jpsiVertex = jpsiVertexCollection->at(jpsiIndex).position();
           j_vtx_x = jpsiVertexCollection->at(jpsiIndex).x();
           j_vtx_y = jpsiVertexCollection->at(jpsiIndex).y();
           j_vtx_z = jpsiVertexCollection->at(jpsiIndex).z();
         }
       }      
     }//end passHLT

     muonsMatched = matchMuon1 && matchMuon2;

     //DoubleMu4_JpsiTrk_Displaced

     hltTriggerObjects_pT.clear();
     hltTriggerObjects_eta.clear();
     hltTriggerObjects_phi.clear();
     matchTrackToJpsi = false;

     if (passHLT_DoubleMu4_JpsiTrk_Displaced) {
       if (triggerObjectsSummary.isValid()) {
         size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltJpsiTkVertexFilter", "", hltTag_) );
         trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
         if (filterIndex < (*triggerObjectsSummary).sizeFilters()) {
           const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
           for (size_t j = 0; j < keys.size(); j++) {
             trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]]; //0 Mu1_A 1 Mu2_A 2 Tk_A 3 Mu1_B 4 Mu2_B 5 Tk_B ... (can be same muon pair and different associated track)
             hltTriggerObjects_pT.push_back(foundObject.pt()); 
             hltTriggerObjects_eta.push_back(foundObject.eta()); 
             hltTriggerObjects_phi.push_back(foundObject.phi());
           }
         }
       }//end triggerObjects

       float dPtRelThreshold = 0.01;
       float dRthreshold = 0.01;
       for (unsigned int i = 0; i < hltTriggerObjects_pT.size(); i=i+3) { //loop on groups Mu1 Mu2 Tk
         //check that Jpsi in hltJpsiTrk is the same from hltJpsi
         matchHLTMuon1 = (deltaPtOverPt(hlt_Mu1_pT, hltTriggerObjects_pT[i+0]) < dPtRelThreshold && deltaR(hlt_Mu1_eta, hlt_Mu1_phi, hltTriggerObjects_eta[i+0], hltTriggerObjects_phi[i+0]) < dRthreshold);
         matchHLTMuon2 = (deltaPtOverPt(hlt_Mu2_pT, hltTriggerObjects_pT[i+1]) < dPtRelThreshold && deltaR(hlt_Mu2_eta, hlt_Mu2_phi, hltTriggerObjects_eta[i+1], hltTriggerObjects_phi[i+1]) < dRthreshold);
         //match trk from hltJpsiTrk to gen track
         float dR1 = deltaR(track1_eta, track1_phi, hltTriggerObjects_eta[i+2], hltTriggerObjects_phi[i+2]);
         float dR2 = deltaR(track2_eta, track2_phi, hltTriggerObjects_eta[i+2], hltTriggerObjects_phi[i+2]);
         if (dR1 < dRthreshold) {
           hlt_Track1_pT = hltTriggerObjects_pT[i+2];
           hlt_Track1_eta = hltTriggerObjects_eta[i+2];
           hlt_Track1_phi = hltTriggerObjects_phi[i+2];
           hlt_Track1_deltaR = dR1;
           hltMatch_Track1 = true;
         } 
         else if (dR2 < dRthreshold) {
           hlt_Track2_pT = hltTriggerObjects_pT[i+2];
           hlt_Track2_eta = hltTriggerObjects_eta[i+2];
           hlt_Track2_phi = hltTriggerObjects_phi[i+2];
           hlt_Track2_deltaR = dR2;
           hltMatch_Track2 = true;
         } 
         else {
           hltTriggerFake_pt.push_back(hltTriggerObjects_pT[i+2]);
           hltTriggerFake_eta.push_back(hltTriggerObjects_eta[i+2]);
         }
         hltTriggerTrack_pt.push_back(hltTriggerObjects_pT[i+2]);
         hltTriggerTrack_eta.push_back(hltTriggerObjects_eta[i+2]);
       } 
       matchTrack = hltMatch_Track1 || hltMatch_Track2;
       matchTrackToJpsi = matchHLTMuon1 && matchHLTMuon2 && matchTrack;
     }//end passHLTJpsiTrk
   
     hltMatch_Track = muonsMatched && matchTrackToJpsi;

     //DoubleMu4_JpsiTrk_Displaced_GPU
  
     hltTriggerObjects_pT.clear();
     hltTriggerObjects_eta.clear();
     hltTriggerObjects_phi.clear();
     matchTrackToJpsi = false;

     if (passHLT_DoubleMu4_JpsiTrk_Displaced_GPU) {
       if (triggerObjectsSummary.isValid()) {
         size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltJpsiTkVertexFilterGPU", "", hltTag_) );
         trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
         if (filterIndex < (*triggerObjectsSummary).sizeFilters()) {
           const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
           for (size_t j = 0; j < keys.size(); j++) {
             trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]]; //0 Mu1_A 1 Mu2_A 2 Tk_A 3 Mu1_B 4 Mu2_B 5 Tk_B ... (can be same muon pair and different associated track)
             hltTriggerObjects_pT.push_back(foundObject.pt());
             hltTriggerObjects_eta.push_back(foundObject.eta());
             hltTriggerObjects_phi.push_back(foundObject.phi());
           }
         }
       }//end triggerObjects

       float dPtRelThreshold = 0.01;
       float dRthreshold = 0.01;
       for (unsigned int i = 0; i < hltTriggerObjects_pT.size(); i=i+3) { //loop on groups Mu1 Mu2 Tk

         matchHLTMuon1 = (deltaPtOverPt(hlt_Mu1_pT, hltTriggerObjects_pT[i+0]) < dPtRelThreshold && deltaR(hlt_Mu1_eta, hlt_Mu1_phi, hltTriggerObjects_eta[i+0], hltTriggerObjects_phi[i+0]) < dRthreshold);
         matchHLTMuon2 = (deltaPtOverPt(hlt_Mu2_pT, hltTriggerObjects_pT[i+1]) < dPtRelThreshold && deltaR(hlt_Mu2_eta, hlt_Mu2_phi, hltTriggerObjects_eta[i+1], hltTriggerObjects_phi[i+1]) < dRthreshold);
         float dR1 = deltaR(track1_eta, track1_phi, hltTriggerObjects_eta[i+2], hltTriggerObjects_phi[i+2]);
         float dR2 = deltaR(track2_eta, track2_phi, hltTriggerObjects_eta[i+2], hltTriggerObjects_phi[i+2]);
         if (dR1 < dRthreshold) {
           gpu_Track1_pT = hltTriggerObjects_pT[i+2];
           gpu_Track1_eta = hltTriggerObjects_eta[i+2];
           gpu_Track1_phi = hltTriggerObjects_phi[i+2];
           gpu_Track1_deltaR = dR1;
           gpuMatch_Track1 = true;
         }
         else if (dR2 < dRthreshold) {
           gpu_Track2_pT = hltTriggerObjects_pT[i+2];
           gpu_Track2_eta = hltTriggerObjects_eta[i+2];
           gpu_Track2_phi = hltTriggerObjects_phi[i+2];
           gpu_Track2_deltaR = dR2;
           gpuMatch_Track2 = true;
         }
         else {
           gpuTriggerFake_pt.push_back(hltTriggerObjects_pT[i+2]);
           gpuTriggerFake_eta.push_back(hltTriggerObjects_eta[i+2]);
         }
         gpuTriggerTrack_pt.push_back(hltTriggerObjects_pT[i+2]);
         gpuTriggerTrack_eta.push_back(hltTriggerObjects_eta[i+2]);
       }
       matchTrack = gpuMatch_Track1 || gpuMatch_Track2;
       matchTrackToJpsi = matchHLTMuon1 && matchHLTMuon2 && matchTrack;
     }//end passHLTJpsiTrk_GPU

     gpuMatch_Track = muonsMatched && matchTrackToJpsi;

     //Intemerdiate steps for DoubleMu4_JpsiTrk_Displaced
    
     edm::Handle<vector<reco::Track>> legacyPixelTrackCollection;
     edm::Handle<vector<reco::Track>> legacyIter0TrackCollection;
     edm::Handle<vector<reco::Track>> legacyIter1TrackCollection;
     edm::Handle<vector<reco::Track>> legacyIter2TrackCollection;
     edm::Handle<vector<reco::Track>> legacy3RecoTrackCollection;
     edm::Handle<vector<reco::Track>> legacy2RecoTrackCollection;

     iEvent.getByToken(legacyPixelTracksToken_, legacyPixelTrackCollection);
     iEvent.getByToken(legacyIter0TracksToken_, legacyIter0TrackCollection);
     iEvent.getByToken(legacyIter1TracksToken_, legacyIter1TrackCollection);
     iEvent.getByToken(legacyIter2TracksToken_, legacyIter2TrackCollection);
     iEvent.getByToken(legacy3RecoTracksToken_, legacy3RecoTrackCollection);
     iEvent.getByToken(legacy2RecoTracksToken_, legacy2RecoTrackCollection);

     if (muonsMatched) {

       float dRthreshold = 0.03;

       legPix_collSize = legacyPixelTrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = legacyPixelTrackCollection->begin(); iTrack != legacyPixelTrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           legPix_Track1_pT = iTrack->pt();
           legPix_Track1_eta = iTrack->eta();
           legPix_Track1_phi = iTrack->phi();
           legPix_Track1_d0 = iTrack->dxy();
           legPix_Track1_deltaR = dR1;
           legPix_Track1_dxy = dxy;
           legPix_Track1_dz = dz;
           legPix_Track1_dEtaMu1 = dEtaMu1;
           legPix_Track1_dPhiMu1 = dPhiMu1;
           legPix_Track1_dEtaMu2 = dEtaMu2;
           legPix_Track1_dPhiMu2 = dPhiMu2;
           legPix_Track1_nHits = nHits;
           legPix_Track1_normChi2 = normChi2;
           legPixMatch_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           legPix_Track2_pT = iTrack->pt();
           legPix_Track2_eta = iTrack->eta();
           legPix_Track2_phi = iTrack->phi();
           legPix_Track2_d0 = iTrack->dxy();
           legPix_Track2_deltaR = dR2;
           legPix_Track2_dxy = dxy;
           legPix_Track2_dz = dz;
           legPix_Track2_dEtaMu1 = dEtaMu1;
           legPix_Track2_dPhiMu1 = dPhiMu1;
           legPix_Track2_dEtaMu2 = dEtaMu2;
           legPix_Track2_dPhiMu2 = dPhiMu2;
           legPix_Track2_nHits = nHits;
           legPix_Track2_normChi2 = normChi2;
           legPixMatch_Track2 = true;
         }
       }//end legacyPixelTracks

       dRthreshold = 0.01;

       legI0_collSize = legacyIter0TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = legacyIter0TrackCollection->begin(); iTrack != legacyIter0TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           legI0_Track1_pT = iTrack->pt();
           legI0_Track1_eta = iTrack->eta();
           legI0_Track1_phi = iTrack->phi();
           legI0_Track1_d0 = iTrack->dxy();
           legI0Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           legI0_Track2_pT = iTrack->pt();
           legI0_Track2_eta = iTrack->eta();
           legI0_Track2_phi = iTrack->phi();
           legI0_Track2_d0 = iTrack->dxy();
           legI0Match_Track2 = true;
         }
       }//end Iter0

      legI1_collSize = legacyIter1TrackCollection->size();
      for(vector<reco::Track>::const_iterator iTrack = legacyIter1TrackCollection->begin(); iTrack != legacyIter1TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           legI1_Track1_pT = iTrack->pt();
           legI1_Track1_eta = iTrack->eta();
           legI1_Track1_phi = iTrack->phi();
           legI1_Track1_d0 = iTrack->dxy();
           legI1Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           legI1_Track2_pT = iTrack->pt();
           legI1_Track2_eta = iTrack->eta();
           legI1_Track2_phi = iTrack->phi();
           legI1_Track2_d0 = iTrack->dxy();
           legI1Match_Track2 = true;
         }
       }//end Iter1

      legI2_collSize = legacyIter2TrackCollection->size();
      for(vector<reco::Track>::const_iterator iTrack = legacyIter2TrackCollection->begin(); iTrack != legacyIter2TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           legI2_Track1_pT = iTrack->pt();
           legI2_Track1_eta = iTrack->eta();
           legI2_Track1_phi = iTrack->phi();
           legI2_Track1_d0 = iTrack->dxy();
           legI2Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           legI2_Track2_pT = iTrack->pt();
           legI2_Track2_eta = iTrack->eta();
           legI2_Track2_phi = iTrack->phi();
           legI2_Track2_d0 = iTrack->dxy();
           legI2Match_Track2 = true;
         }
       }//end Iter2

      for(vector<reco::Track>::const_iterator iTrack = legacy3RecoTrackCollection->begin(); iTrack != legacy3RecoTrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           leg3R_Track1_pT = iTrack->pt();
           leg3R_Track1_eta = iTrack->eta();
           leg3R_Track1_phi = iTrack->phi();
           leg3R_Track1_d0 = iTrack->dxy();
           leg3RMatch_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           leg3R_Track2_pT = iTrack->pt();
           leg3R_Track2_eta = iTrack->eta();
           leg3R_Track2_phi = iTrack->phi();
           leg3R_Track2_d0 = iTrack->dxy();
           leg3RMatch_Track2 = true;
         }
       }//end TripletRecovery

      for(vector<reco::Track>::const_iterator iTrack = legacy2RecoTrackCollection->begin(); iTrack != legacy2RecoTrackCollection->end(); iTrack++) {
        float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
        float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
        if (dR1 < dRthreshold) {
          leg2R_Track1_pT = iTrack->pt();
          leg2R_Track1_eta = iTrack->eta();
          leg2R_Track1_phi = iTrack->phi();
          leg2R_Track1_d0 = iTrack->dxy();
          leg2R_Track1_d0err = iTrack->dxyError();
          leg2RMatch_Track1 = true;
        }
        if (dR2 < dRthreshold) {
          leg2R_Track2_pT = iTrack->pt();
          leg2R_Track2_eta = iTrack->eta();
          leg2R_Track2_phi = iTrack->phi();
          leg2R_Track2_d0 = iTrack->dxy();
          leg2R_Track2_d0err = iTrack->dxyError();
          leg2RMatch_Track2 = true;
        }
        //save all track collection
        float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
        float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
        if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
        float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
        float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
        if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
        float dxy = iTrack->dxy(jpsiVertex);
        float dz = fabs( iTrack->dz(jpsiVertex) );
        bool match = dR1 < dRthreshold || dR2 < dRthreshold;
        prefilter_legTks_pT.push_back(iTrack->pt());
        prefilter_legTks_eta.push_back(iTrack->eta());
        prefilter_legTks_phi.push_back(iTrack->phi());
        prefilter_legTks_dEta1.push_back(dEtaMu1);
        prefilter_legTks_dPhi1.push_back(dPhiMu1);
        prefilter_legTks_dEta2.push_back(dEtaMu2);
        prefilter_legTks_dPhi2.push_back(dPhiMu2);
        prefilter_legTks_dxy.push_back(dxy);
        prefilter_legTks_dz.push_back(dz);
        prefilter_legTks_d0.push_back(iTrack->dxy());
        prefilter_legTks_d0err.push_back(iTrack->dxyError());
        prefilter_legTks_bool.push_back(match);
      }//end DoubletRecovery

     }//end muonsMatched for HLT_DoubleMu4_JpsiTrk_Displaced

     //Intermediate steps for HLT_DoubleMu4_JpsiTrk_Displaced_GPU

     edm::Handle<vector<reco::Track>> myPixelTracksCollection;
     edm::Handle<vector<reco::Track>> regionalPixelTrackCollection;
     edm::Handle<vector<reco::Track>> ctfIter0TrackCollection;
     edm::Handle<vector<reco::Track>> triggerIter0TrackCollection;
     edm::Handle<vector<reco::Track>> ctfIter1TrackCollection;
     edm::Handle<vector<reco::Track>> triggerIter1TrackCollection;
     edm::Handle<vector<reco::Track>> ctfIter2TrackCollection;
     edm::Handle<vector<reco::Track>> triggerIter2TrackCollection;

     iEvent.getByToken(myPixelTracksToken_, myPixelTracksCollection);
     iEvent.getByToken(regionalPixelTracksToken_, regionalPixelTrackCollection);
     iEvent.getByToken(ctfIter0TracksToken_, ctfIter0TrackCollection);
     iEvent.getByToken(triggerIter0TracksToken_, triggerIter0TrackCollection);
     iEvent.getByToken(ctfIter1TracksToken_, ctfIter1TrackCollection);
     iEvent.getByToken(triggerIter1TracksToken_, triggerIter1TrackCollection);
     iEvent.getByToken(ctfIter2TracksToken_, ctfIter2TrackCollection);
     iEvent.getByToken(triggerIter2TracksToken_, triggerIter2TrackCollection);

     if (muonsMatched) { //enough to be sure that trackCollections are not empty, but should be checked
       float dRthreshold = 0.03;

       myPix_collSize = myPixelTracksCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = myPixelTracksCollection->begin(); iTrack != myPixelTracksCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           myPix_Track1_pT = iTrack->pt();
           myPix_Track1_eta = iTrack->eta();
           myPix_Track1_phi = iTrack->phi();
           myPix_Track1_d0 = iTrack->dxy();
           myPix_Track1_deltaR = dR1;
           myPix_Track1_dxy = dxy;
           myPix_Track1_dz = dz;
           myPix_Track1_dEtaMu1 = dEtaMu1;
           myPix_Track1_dPhiMu1 = dPhiMu1;
           myPix_Track1_dEtaMu2 = dEtaMu2;
           myPix_Track1_dPhiMu2 = dPhiMu2;
           myPix_Track1_nHits = nHits;
           myPix_Track1_normChi2 = normChi2;
           myPixMatch_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           myPix_Track2_pT = iTrack->pt();
           myPix_Track2_eta = iTrack->eta();
           myPix_Track2_phi = iTrack->phi();
           myPix_Track2_d0 = iTrack->dxy();
           myPix_Track2_deltaR = dR2;
           myPix_Track2_dxy = dxy;
           myPix_Track2_dz = dz;
           myPix_Track2_dEtaMu1 = dEtaMu1;
           myPix_Track2_dPhiMu1 = dPhiMu1;
           myPix_Track2_dEtaMu2 = dEtaMu2;
           myPix_Track2_dPhiMu2 = dPhiMu2;
           myPix_Track2_nHits = nHits;
           myPix_Track2_normChi2 = normChi2;
           myPixMatch_Track2 = true;
         }
       }//end myPixelTracks

       reg_collSize = regionalPixelTrackCollection->size();
       myPix_collSize == 0 ? regionalFrac_collSize = 0 : regionalFrac_collSize = (float)reg_collSize/(float)myPix_collSize;

       for(vector<reco::Track>::const_iterator iTrack = regionalPixelTrackCollection->begin(); iTrack != regionalPixelTrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           reg_Track1_pT = iTrack->pt();
           reg_Track1_eta = iTrack->eta();
           reg_Track1_phi = iTrack->phi();
           reg_Track1_d0 = iTrack->dxy();
           reg_Track1_deltaR = dR1;
           reg_Track1_dxy = dxy;
           reg_Track1_dz = dz;
           reg_Track1_dEtaMu1 = dEtaMu1;
           reg_Track1_dPhiMu1 = dPhiMu1;
           reg_Track1_dEtaMu2 = dEtaMu2;
           reg_Track1_dPhiMu2 = dPhiMu2;
           reg_Track1_nHits = nHits;
           reg_Track1_normChi2 = normChi2;
           regMatch_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           reg_Track2_pT = iTrack->pt();
           reg_Track2_eta = iTrack->eta();
           reg_Track2_phi = iTrack->phi();
           reg_Track2_d0 = iTrack->dxy();
           reg_Track2_deltaR = dR2;
           reg_Track2_dxy = dxy;
           reg_Track2_dz = dz;
           reg_Track2_dEtaMu1 = dEtaMu1;
           reg_Track2_dPhiMu1 = dPhiMu1;
           reg_Track2_dEtaMu2 = dEtaMu2;
           reg_Track2_dPhiMu2 = dPhiMu2;
           reg_Track2_nHits = nHits;
           reg_Track2_normChi2 = normChi2;
           regMatch_Track2 = true;
         }
       }//end regionalPixelTracks
       
       dRthreshold = 0.01; //from now on - high quality tracks

       ctf0_collSize = ctfIter0TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = ctfIter0TrackCollection->begin(); iTrack != ctfIter0TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           ctf0_Track1_pT = iTrack->pt();
           ctf0_Track1_eta = iTrack->eta();
           ctf0_Track1_phi = iTrack->phi();
           ctf0_Track1_d0 = iTrack->dxy();
           ctf0_Track1_deltaR = dR1;
           ctf0_Track1_dxy = dxy;
           ctf0_Track1_dz = dz;
           ctf0_Track1_dEtaMu1 = dEtaMu1;
           ctf0_Track1_dPhiMu1 = dPhiMu1;
           ctf0_Track1_dEtaMu2 = dEtaMu2;
           ctf0_Track1_dPhiMu2 = dPhiMu2;
           ctf0_Track1_nHits = nHits;
           ctf0_Track1_normChi2 = normChi2;
           ctf0Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           ctf0_Track2_pT = iTrack->pt();
           ctf0_Track2_eta = iTrack->eta();
           ctf0_Track2_phi = iTrack->phi();
           ctf0_Track2_d0 = iTrack->dxy();
           ctf0_Track2_deltaR = dR2;
           ctf0_Track2_dxy = dxy;
           ctf0_Track2_dz = dz;
           ctf0_Track2_dEtaMu1 = dEtaMu1;
           ctf0_Track2_dPhiMu1 = dPhiMu1;
           ctf0_Track2_dEtaMu2 = dEtaMu2;
           ctf0_Track2_dPhiMu2 = dPhiMu2;
           ctf0_Track2_nHits = nHits;
           ctf0_Track2_normChi2 = normChi2;
           ctf0Match_Track2 = true;
         }
       }//end ctfIter0

       iter0_collSize = triggerIter0TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = triggerIter0TrackCollection->begin(); iTrack != triggerIter0TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           iter0_Track1_pT = iTrack->pt();
           iter0_Track1_eta = iTrack->eta();
           iter0_Track1_phi = iTrack->phi();
           iter0_Track1_d0 = iTrack->dxy();
           iter0Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           iter0_Track2_pT = iTrack->pt();
           iter0_Track2_eta = iTrack->eta();
           iter0_Track2_phi = iTrack->phi();
           iter0_Track2_d0 = iTrack->dxy();
           iter0Match_Track2 = true;
         }
       }//end Iter0

       ctf1_collSize = ctfIter1TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = ctfIter1TrackCollection->begin(); iTrack != ctfIter1TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           ctf1_Track1_pT = iTrack->pt();
           ctf1_Track1_eta = iTrack->eta();
           ctf1_Track1_phi = iTrack->phi();
           ctf1_Track1_d0 = iTrack->dxy();
           ctf1_Track1_deltaR = dR1;
           ctf1_Track1_dxy = dxy;
           ctf1_Track1_dz = dz;
           ctf1_Track1_dEtaMu1 = dEtaMu1;
           ctf1_Track1_dPhiMu1 = dPhiMu1;
           ctf1_Track1_dEtaMu2 = dEtaMu2;
           ctf1_Track1_dPhiMu2 = dPhiMu2;
           ctf1_Track1_nHits = nHits;
           ctf1_Track1_normChi2 = normChi2;
           ctf1Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           ctf1_Track2_pT = iTrack->pt();
           ctf1_Track2_eta = iTrack->eta();
           ctf1_Track2_phi = iTrack->phi();
           ctf1_Track2_d0 = iTrack->dxy();
           ctf1_Track2_deltaR = dR2;
           ctf1_Track2_dxy = dxy;
           ctf1_Track2_dz = dz;
           ctf1_Track2_dEtaMu1 = dEtaMu1;
           ctf1_Track2_dPhiMu1 = dPhiMu1;
           ctf1_Track2_dEtaMu2 = dEtaMu2;
           ctf1_Track2_dPhiMu2 = dPhiMu2;
           ctf1_Track2_nHits = nHits;
           ctf1_Track2_normChi2 = normChi2;
           ctf1Match_Track2 = true;
         }
       }//end ctfIter1

       iter1_collSize = triggerIter1TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = triggerIter1TrackCollection->begin(); iTrack != triggerIter1TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           iter1_Track1_pT = iTrack->pt();
           iter1_Track1_eta = iTrack->eta();
           iter1_Track1_phi = iTrack->phi();
           iter1_Track1_d0 = iTrack->dxy();
           iter1Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           iter1_Track2_pT = iTrack->pt();
           iter1_Track2_eta = iTrack->eta();
           iter1_Track2_phi = iTrack->phi();
           iter1_Track2_d0 = iTrack->dxy();
           iter1Match_Track2 = true;
         }
       }//end Iter1

       ctf2_collSize = ctfIter2TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = ctfIter2TrackCollection->begin(); iTrack != ctfIter2TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         float dxy = iTrack->dxy(jpsiVertex);
         float dz = fabs( iTrack->dz(jpsiVertex) );
         float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
         float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
         if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
         float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
         float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
         if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
         unsigned int nHits = iTrack->numberOfValidHits();
         float normChi2 = iTrack->normalizedChi2();
         if (dR1 < dRthreshold) {
           ctf2_Track1_pT = iTrack->pt();
           ctf2_Track1_eta = iTrack->eta();
           ctf2_Track1_phi = iTrack->phi();
           ctf2_Track1_d0 = iTrack->dxy();
           ctf2_Track1_deltaR = dR1;
           ctf2_Track1_dxy = dxy;
           ctf2_Track1_dz = dz;
           ctf2_Track1_dEtaMu1 = dEtaMu1;
           ctf2_Track1_dPhiMu1 = dPhiMu1;
           ctf2_Track1_dEtaMu2 = dEtaMu2;
           ctf2_Track1_dPhiMu2 = dPhiMu2;
           ctf2_Track1_nHits = nHits;
           ctf2_Track1_normChi2 = normChi2;
           ctf2Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           ctf2_Track2_pT = iTrack->pt();
           ctf2_Track2_eta = iTrack->eta();
           ctf2_Track2_phi = iTrack->phi();
           ctf2_Track2_d0 = iTrack->dxy();
           ctf2_Track2_deltaR = dR2;
           ctf2_Track2_dxy = dxy;
           ctf2_Track2_dz = dz;
           ctf2_Track2_dEtaMu1 = dEtaMu1;
           ctf2_Track2_dPhiMu1 = dPhiMu1;
           ctf2_Track2_dEtaMu2 = dEtaMu2;
           ctf2_Track2_dPhiMu2 = dPhiMu2;
           ctf2_Track2_nHits = nHits;
           ctf2_Track2_normChi2 = normChi2;
           ctf2Match_Track2 = true;
         }
       }//end ctfIter2

       iter2_collSize = triggerIter2TrackCollection->size();
       for(vector<reco::Track>::const_iterator iTrack = triggerIter2TrackCollection->begin(); iTrack != triggerIter2TrackCollection->end(); iTrack++) {
         float dR1 = deltaR(track1_eta, track1_phi, iTrack->eta(), iTrack->phi());
         float dR2 = deltaR(track2_eta, track2_phi, iTrack->eta(), iTrack->phi());
         if (dR1 < dRthreshold) {
           iter2_Track1_pT = iTrack->pt();
           iter2_Track1_eta = iTrack->eta();
           iter2_Track1_phi = iTrack->phi();
           iter2_Track1_d0 = iTrack->dxy();
           iter2_Track1_d0err = iTrack->dxyError();
           iter2Match_Track1 = true;
         }
         if (dR2 < dRthreshold) {
           iter2_Track2_pT = iTrack->pt();
           iter2_Track2_eta = iTrack->eta();
           iter2_Track2_phi = iTrack->phi();
           iter2_Track2_d0 = iTrack->dxy();
           iter2_Track2_d0err = iTrack->dxyError();
           iter2Match_Track2 = true;
         }
        float dEtaMu1 = fabs( iTrack->eta() - hlt_Mu1_eta );
        float dPhiMu1 = fabs( iTrack->phi() - hlt_Mu1_phi );
        if ( dPhiMu1 > ROOT::Math::Pi() ) dPhiMu1 = fabs( 2*ROOT::Math::Pi() - dPhiMu1 );
        float dEtaMu2 = fabs( iTrack->eta() - hlt_Mu2_eta );
        float dPhiMu2 = fabs( iTrack->phi() - hlt_Mu2_phi );
        if ( dPhiMu2 > ROOT::Math::Pi() ) dPhiMu2 = fabs( 2*ROOT::Math::Pi() - dPhiMu2 );
        float dxy = iTrack->dxy(jpsiVertex);
        float dz = fabs( iTrack->dz(jpsiVertex) );
        bool match = dR1 < dRthreshold || dR2 < dRthreshold;
        prefilter_newTks_pT.push_back(iTrack->pt());
        prefilter_newTks_eta.push_back(iTrack->eta());
        prefilter_newTks_phi.push_back(iTrack->phi());
        prefilter_newTks_dEta1.push_back(dEtaMu1);
        prefilter_newTks_dPhi1.push_back(dPhiMu1);
        prefilter_newTks_dEta2.push_back(dEtaMu2);
        prefilter_newTks_dPhi2.push_back(dPhiMu2);
        prefilter_newTks_dxy.push_back(dxy);
        prefilter_newTks_dz.push_back(dz);
        prefilter_newTks_d0.push_back(iTrack->dxy());
        prefilter_newTks_d0err.push_back(iTrack->dxyError());
        prefilter_newTks_bool.push_back(match);
       }//end Iter2

     }//end muonsMatched

     //Print

     outTree->Fill();

   }//end Bs for

}


float
VanillaHLTAnalyzer::deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dPhi = fabs(phi1 - phi2);
  if ( dPhi > ROOT::Math::Pi() )  dPhi = 2*ROOT::Math::Pi() - dPhi;

  return sqrt( (eta1-eta2)*(eta1-eta2) + dPhi*dPhi );
}

float
VanillaHLTAnalyzer::deltaPtOverPt(float pt1, float pt2)
{
  float deltaPt = fabs(pt1 - pt2);
   float deltaPtOverPt = 0;

  if (pt1 == 0) {
    if (pt2 == 0) {
      deltaPtOverPt = 0;
    } else {
      deltaPtOverPt = deltaPt/pt2;
      }
  } else {
    deltaPtOverPt = deltaPt/pt1;
    }
 
  return deltaPtOverPt;
}

void
VanillaHLTAnalyzer::printProgeny( const reco::GenParticle& part)
{
  cout<<part.pdgId()<<"(";
  for (uint iDaug=0; iDaug<part.numberOfDaughters(); ++iDaug) {
    if (iDaug>0) cout<<" ";
    printProgeny( *((reco::GenParticle*)part.daughter(iDaug)) );
  }
  cout<<")";
}

FreeTrajectoryState 
VanillaHLTAnalyzer::initialFreeState( const reco::Track& tk, const MagneticField* field)
{
  Basic3DVector<float> pos( tk.vertex());
  GlobalPoint gpos( pos);
  Basic3DVector<float> mom( tk.momentum());
  GlobalVector gmom( mom);
  GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err( tk.covariance());
  return FreeTrajectoryState( par, err);
}

void
VanillaHLTAnalyzer::beginRun(const edm::Run & run, const edm::EventSetup & iSetup)
{
  // bool changed(true);
  // if (hltPrescale_.init(run,iSetup,"HLT",changed)) {
  //   // if init returns TRUE, initialisation has succeeded!
  //   if (changed) {
  //     // The HLT config has actually changed wrt the previous Run
  //     std::cout << "Initializing HLTConfigProvider"  << std::endl;
  //   }
  // } 
  // else {
  //   std::cout << " HLT config extraction failure with process name HLT" << std::endl;
  // }
}

// ------------ method called once each job just before starting event loop  ------------
void 
VanillaHLTAnalyzer::beginJob()
{
  outTree = outfile_-> make<TTree>("ntupleTree","ntupleTree");

  outTree->Branch("hltNames"    , &hltNames     );
  outTree->Branch("hltResults"  , &hltResults   );

  outTree->Branch("Muon1_pT" , &Muon1_pT  );
  outTree->Branch("Muon1_eta", &Muon1_eta );
  outTree->Branch("Muon1_phi", &Muon1_phi );
  outTree->Branch("Muon1_cha", &Muon1_cha );
  outTree->Branch("Muon2_pT" , &Muon2_pT  );
  outTree->Branch("Muon2_eta", &Muon2_eta );
  outTree->Branch("Muon2_phi", &Muon2_phi );
  outTree->Branch("Muon2_cha", &Muon2_cha );

  outTree->Branch("track1_pT" , &track1_pT  );
  outTree->Branch("track1_eta", &track1_eta );
  outTree->Branch("track1_phi", &track1_phi );
  outTree->Branch("track1_cha", &track1_cha );
  outTree->Branch("track2_pT" , &track2_pT  );
  outTree->Branch("track2_eta", &track2_eta );
  outTree->Branch("track2_phi", &track2_phi );
  outTree->Branch("track2_cha", &track2_cha );

  outTree->Branch("bsFlightLength2D", &bsFlightLength2D );
  outTree->Branch("bsFlightLength3D", &bsFlightLength3D );

  outTree->Branch("passHLT_DoubleMu4_Jpsi_Displaced", &passHLT_DoubleMu4_Jpsi_Displaced );
  outTree->Branch("passHLT_DoubleMu4_JpsiTrk_Displaced", &passHLT_DoubleMu4_JpsiTrk_Displaced );
  outTree->Branch("passHLT_DoubleMu4_JpsiTrk_Displaced_GPU", &passHLT_DoubleMu4_JpsiTrk_Displaced_GPU );

  outTree->Branch("muonsMatched", &muonsMatched );

  outTree->Branch("hltMatch_Track1", &hltMatch_Track1 );
  outTree->Branch("hltMatch_Track2", &hltMatch_Track2 );
  outTree->Branch("gpuMatch_Track1", &gpuMatch_Track1 );
  outTree->Branch("gpuMatch_Track2", &gpuMatch_Track2 );

  outTree->Branch("hltTriggerTrack_pt", &hltTriggerTrack_pt );
  outTree->Branch("hltTriggerTrack_eta", &hltTriggerTrack_eta );
  outTree->Branch("hltTriggerFake_pt", &hltTriggerFake_pt );
  outTree->Branch("hltTriggerFake_eta", &hltTriggerFake_eta );
  outTree->Branch("gpuTriggerTrack_pt", &gpuTriggerTrack_pt );
  outTree->Branch("gpuTriggerTrack_eta", &gpuTriggerTrack_eta );
  outTree->Branch("gpuTriggerFake_pt", &gpuTriggerFake_pt );
  outTree->Branch("gpuTriggerFake_eta", &gpuTriggerFake_eta );

  outTree->Branch("legPixMatch_Track1", &legPixMatch_Track1 );
  outTree->Branch("legPixMatch_Track2", &legPixMatch_Track2 );
  outTree->Branch("legI0Match_Track1", &legI0Match_Track1 );
  outTree->Branch("legI0Match_Track2", &legI0Match_Track2 );
  outTree->Branch("legI1Match_Track1", &legI1Match_Track1 );
  outTree->Branch("legI1Match_Track2", &legI1Match_Track2 );
  outTree->Branch("legI2Match_Track1", &legI2Match_Track1 );
  outTree->Branch("legI2Match_Track2", &legI2Match_Track2 );
  outTree->Branch("leg3RMatch_Track1", &leg3RMatch_Track1 );
  outTree->Branch("leg3RMatch_Track2", &leg3RMatch_Track2 );
  outTree->Branch("leg2RMatch_Track1", &leg2RMatch_Track1 );
  outTree->Branch("leg2RMatch_Track2", &leg2RMatch_Track2 );

  outTree->Branch("myPixMatch_Track1", &myPixMatch_Track1 );
  outTree->Branch("myPixMatch_Track2", &myPixMatch_Track2 );
  outTree->Branch("regMatch_Track1", &regMatch_Track1 );
  outTree->Branch("regMatch_Track2", &regMatch_Track2 );
  outTree->Branch("ctf0Match_Track1", &ctf0Match_Track1 );
  outTree->Branch("ctf0Match_Track2", &ctf0Match_Track2 );
  outTree->Branch("iter0Match_Track1", &iter0Match_Track1 );
  outTree->Branch("iter0Match_Track2", &iter0Match_Track2 );
  outTree->Branch("ctf1Match_Track1", &ctf1Match_Track1 );
  outTree->Branch("ctf1Match_Track2", &ctf1Match_Track2 );
  outTree->Branch("iter1Match_Track1", &iter1Match_Track1 );
  outTree->Branch("iter1Match_Track2", &iter1Match_Track2 );
  outTree->Branch("ctf2Match_Track1", &ctf2Match_Track1 );
  outTree->Branch("ctf2Match_Track2", &ctf2Match_Track2 );
  outTree->Branch("iter2Match_Track1", &iter2Match_Track1 );
  outTree->Branch("iter2Match_Track2", &iter2Match_Track2 );

  outTree->Branch("hltMatch_Track", &hltMatch_Track );
  outTree->Branch("gpuMatch_Track", &gpuMatch_Track );

  outTree->Branch("j_vtx_x", &j_vtx_x );
  outTree->Branch("j_vtx_y", &j_vtx_y );
  outTree->Branch("j_vtx_z", &j_vtx_z );
  outTree->Branch("j_vtx_size", &j_vtx_size );

  outTree->Branch("hlt_Mu1_pT", &hlt_Mu1_pT );
  outTree->Branch("hlt_Mu1_eta", &hlt_Mu1_eta );
  outTree->Branch("hlt_Mu1_phi", &hlt_Mu1_phi );
  outTree->Branch("hlt_Mu1_cha", &hlt_Mu1_cha );
  outTree->Branch("hlt_Mu2_pT", &hlt_Mu2_pT );
  outTree->Branch("hlt_Mu2_eta", &hlt_Mu2_eta );
  outTree->Branch("hlt_Mu2_phi", &hlt_Mu2_phi );
  outTree->Branch("hlt_Mu2_cha", &hlt_Mu2_cha );

  outTree->Branch("hlt_Track1_pT", &hlt_Track1_pT );
  outTree->Branch("hlt_Track1_eta", &hlt_Track1_eta );
  outTree->Branch("hlt_Track1_phi", &hlt_Track1_phi );
  outTree->Branch("hlt_Track2_pT", &hlt_Track2_pT );
  outTree->Branch("hlt_Track2_eta", &hlt_Track2_eta );
  outTree->Branch("hlt_Track2_phi", &hlt_Track2_phi );

  outTree->Branch("gpu_Track1_pT", &gpu_Track1_pT );
  outTree->Branch("gpu_Track1_eta", &gpu_Track1_eta );
  outTree->Branch("gpu_Track1_phi", &gpu_Track1_phi );
  outTree->Branch("gpu_Track2_pT", &gpu_Track2_pT );
  outTree->Branch("gpu_Track2_eta", &gpu_Track2_eta );
  outTree->Branch("gpu_Track2_phi", &gpu_Track2_phi );

  outTree->Branch("hlt_Track1_deltaR", &hlt_Track1_deltaR );
  outTree->Branch("hlt_Track2_deltaR", &hlt_Track2_deltaR );

  outTree->Branch("gpu_Track1_deltaR", &gpu_Track1_deltaR );
  outTree->Branch("gpu_Track2_deltaR", &gpu_Track2_deltaR );

  outTree->Branch("Mu1_deltaR", &Mu1_deltaR );
  outTree->Branch("Mu2_deltaR", &Mu2_deltaR );

  outTree->Branch("legPix_Track1_pT", &legPix_Track1_pT );
  outTree->Branch("legPix_Track1_eta", &legPix_Track1_eta );
  outTree->Branch("legPix_Track1_phi", &legPix_Track1_phi );
  outTree->Branch("legPix_Track1_d0", &legPix_Track1_d0 ); 
  outTree->Branch("legPix_Track2_pT", &legPix_Track2_pT );
  outTree->Branch("legPix_Track2_eta", &legPix_Track2_eta );
  outTree->Branch("legPix_Track2_phi", &legPix_Track2_phi );
  outTree->Branch("legPix_Track2_d0", &legPix_Track2_d0 );

  outTree->Branch("legPix_Track1_deltaR", &legPix_Track1_deltaR );
  outTree->Branch("legPix_Track2_deltaR", &legPix_Track2_deltaR );

  outTree->Branch("legPix_Track1_dxy", &legPix_Track1_dxy );
  outTree->Branch("legPix_Track1_dz",  &legPix_Track1_dz );
  outTree->Branch("legPix_Track1_dEtaMu1", &legPix_Track1_dEtaMu1 );
  outTree->Branch("legPix_Track1_dPhiMu1", &legPix_Track1_dPhiMu1 );
  outTree->Branch("legPix_Track1_dEtaMu2", &legPix_Track1_dEtaMu2 );
  outTree->Branch("legPix_Track1_dPhiMu2", &legPix_Track1_dPhiMu2 );
  outTree->Branch("legPix_Track2_dxy", &legPix_Track2_dxy );
  outTree->Branch("legPix_Track2_dz",  &legPix_Track2_dz );
  outTree->Branch("legPix_Track2_dEtaMu1", &legPix_Track2_dEtaMu1 );
  outTree->Branch("legPix_Track2_dPhiMu1", &legPix_Track2_dPhiMu1 );
  outTree->Branch("legPix_Track2_dEtaMu2", &legPix_Track2_dEtaMu2 );
  outTree->Branch("legPix_Track2_dPhiMu2", &legPix_Track2_dPhiMu2 );

  outTree->Branch("legPix_Track1_nHits", &legPix_Track1_nHits );
  outTree->Branch("legPix_Track1_normChi2", &legPix_Track1_normChi2 );
  outTree->Branch("legPix_Track2_nHits", &legPix_Track2_nHits );
  outTree->Branch("legPix_Track2_normChi2", &legPix_Track2_normChi2 );

  outTree->Branch("legPix_collSize", &legPix_collSize );

  outTree->Branch("legI0_Track1_pT", &legI0_Track1_pT );
  outTree->Branch("legI0_Track1_eta", &legI0_Track1_eta );
  outTree->Branch("legI0_Track1_phi", &legI0_Track1_phi );
  outTree->Branch("legI0_Track1_d0", &legI0_Track1_d0 );
  outTree->Branch("legI0_Track2_pT", &legI0_Track2_pT );
  outTree->Branch("legI0_Track2_eta", &legI0_Track2_eta );
  outTree->Branch("legI0_Track2_phi", &legI0_Track2_phi );
  outTree->Branch("legI0_Track2_d0", &legI0_Track2_d0 );

  outTree->Branch("legI0_collSize", &legI0_collSize );

  outTree->Branch("legI1_Track1_pT", &legI1_Track1_pT );
  outTree->Branch("legI1_Track1_eta", &legI1_Track1_eta );
  outTree->Branch("legI1_Track1_phi", &legI1_Track1_phi );
  outTree->Branch("legI1_Track1_d0", &legI1_Track1_d0 );
  outTree->Branch("legI1_Track2_pT", &legI1_Track2_pT );
  outTree->Branch("legI1_Track2_eta", &legI1_Track2_eta );
  outTree->Branch("legI1_Track2_phi", &legI1_Track2_phi );
  outTree->Branch("legI1_Track2_d0", &legI1_Track2_d0 );

  outTree->Branch("legI1_collSize", &legI1_collSize );

  outTree->Branch("legI2_Track1_pT", &legI2_Track1_pT );
  outTree->Branch("legI2_Track1_eta", &legI2_Track1_eta );
  outTree->Branch("legI2_Track1_phi", &legI2_Track1_phi );
  outTree->Branch("legI2_Track1_d0", &legI2_Track1_d0 );
  outTree->Branch("legI2_Track2_pT", &legI2_Track2_pT );
  outTree->Branch("legI2_Track2_eta", &legI2_Track2_eta );
  outTree->Branch("legI2_Track2_phi", &legI2_Track2_phi );
  outTree->Branch("legI2_Track2_d0", &legI2_Track2_d0 );

  outTree->Branch("legI2_collSize", &legI2_collSize );

  outTree->Branch("leg3R_Track1_pT", &leg3R_Track1_pT );
  outTree->Branch("leg3R_Track1_eta", &leg3R_Track1_eta );
  outTree->Branch("leg3R_Track1_phi", &leg3R_Track1_phi );
  outTree->Branch("leg3R_Track1_d0", &leg3R_Track1_d0 );
  outTree->Branch("leg3R_Track2_pT", &leg3R_Track2_pT );
  outTree->Branch("leg3R_Track2_eta", &leg3R_Track2_eta );
  outTree->Branch("leg3R_Track2_phi", &leg3R_Track2_phi );
  outTree->Branch("leg3R_Track2_d0", &leg3R_Track2_d0 );

  outTree->Branch("leg2R_Track1_pT", &leg2R_Track1_pT );
  outTree->Branch("leg2R_Track1_eta", &leg2R_Track1_eta );
  outTree->Branch("leg2R_Track1_phi", &leg2R_Track1_phi );
  outTree->Branch("leg2R_Track1_d0", &leg2R_Track1_d0 );
  outTree->Branch("leg2R_Track1_d0err", &leg2R_Track1_d0err );
  outTree->Branch("leg2R_Track2_pT", &leg2R_Track2_pT );
  outTree->Branch("leg2R_Track2_eta", &leg2R_Track2_eta );
  outTree->Branch("leg2R_Track2_phi", &leg2R_Track2_phi );
  outTree->Branch("leg2R_Track2_d0", &leg2R_Track2_d0 );
  outTree->Branch("leg2R_Track2_d0err", &leg2R_Track2_d0err );

  outTree->Branch("prefilter_legTks_pT",    &prefilter_legTks_pT );
  outTree->Branch("prefilter_legTks_eta",   &prefilter_legTks_eta );
  outTree->Branch("prefilter_legTks_phi",   &prefilter_legTks_phi );
  outTree->Branch("prefilter_legTks_dEta1", &prefilter_legTks_dEta1 );
  outTree->Branch("prefilter_legTks_dPhi1", &prefilter_legTks_dPhi1 );
  outTree->Branch("prefilter_legTks_dEta2", &prefilter_legTks_dEta2 );
  outTree->Branch("prefilter_legTks_dPhi2", &prefilter_legTks_dPhi2 );
  outTree->Branch("prefilter_legTks_dxy",   &prefilter_legTks_dxy );
  outTree->Branch("prefilter_legTks_dz",    &prefilter_legTks_dz );
  outTree->Branch("prefilter_legTks_d0",    &prefilter_legTks_d0 );
  outTree->Branch("prefilter_legTks_d0err", &prefilter_legTks_d0err );
  outTree->Branch("prefilter_legTks_bool",  &prefilter_legTks_bool );

  outTree->Branch("myPix_Track1_pT", &myPix_Track1_pT );
  outTree->Branch("myPix_Track1_eta", &myPix_Track1_eta );
  outTree->Branch("myPix_Track1_phi", &myPix_Track1_phi );
  outTree->Branch("myPix_Track1_d0", &myPix_Track1_d0 );
  outTree->Branch("myPix_Track2_pT", &myPix_Track2_pT );
  outTree->Branch("myPix_Track2_eta", &myPix_Track2_eta );
  outTree->Branch("myPix_Track2_phi", &myPix_Track2_phi );
  outTree->Branch("myPix_Track2_d0", &myPix_Track2_d0 );

  outTree->Branch("myPix_Track1_deltaR", &myPix_Track1_deltaR );
  outTree->Branch("myPix_Track2_deltaR", &myPix_Track2_deltaR );

  outTree->Branch("myPix_Track1_dxy", &myPix_Track1_dxy );
  outTree->Branch("myPix_Track1_dz",  &myPix_Track1_dz );
  outTree->Branch("myPix_Track1_dEtaMu1", &myPix_Track1_dEtaMu1 );
  outTree->Branch("myPix_Track1_dPhiMu1", &myPix_Track1_dPhiMu1 );
  outTree->Branch("myPix_Track1_dEtaMu2", &myPix_Track1_dEtaMu2 );
  outTree->Branch("myPix_Track1_dPhiMu2", &myPix_Track1_dPhiMu2 );
  outTree->Branch("myPix_Track2_dxy", &myPix_Track2_dxy );
  outTree->Branch("myPix_Track2_dz",  &myPix_Track2_dz );
  outTree->Branch("myPix_Track2_dEtaMu1", &myPix_Track2_dEtaMu1 );
  outTree->Branch("myPix_Track2_dPhiMu1", &myPix_Track2_dPhiMu1 );
  outTree->Branch("myPix_Track2_dEtaMu2", &myPix_Track2_dEtaMu2 );
  outTree->Branch("myPix_Track2_dPhiMu2", &myPix_Track2_dPhiMu2 );

  outTree->Branch("myPix_Track1_nHits", &myPix_Track1_nHits );
  outTree->Branch("myPix_Track1_normChi2", &myPix_Track1_normChi2 );
  outTree->Branch("myPix_Track2_nHits", &myPix_Track2_nHits );
  outTree->Branch("myPix_Track2_normChi2", &myPix_Track2_normChi2 );

  outTree->Branch("myPix_collSize", &myPix_collSize );

  outTree->Branch("reg_Track1_pT", &reg_Track1_pT );
  outTree->Branch("reg_Track1_eta", &reg_Track1_eta );
  outTree->Branch("reg_Track1_phi", &reg_Track1_phi );
  outTree->Branch("reg_Track1_d0", &reg_Track1_d0 );
  outTree->Branch("reg_Track2_pT", &reg_Track2_pT );
  outTree->Branch("reg_Track2_eta", &reg_Track2_eta );
  outTree->Branch("reg_Track2_phi", &reg_Track2_phi );
  outTree->Branch("reg_Track2_d0", &reg_Track2_d0 );

  outTree->Branch("reg_Track1_deltaR", &reg_Track1_deltaR );
  outTree->Branch("reg_Track2_deltaR", &reg_Track2_deltaR );

  outTree->Branch("reg_Track1_dxy", &reg_Track1_dxy );
  outTree->Branch("reg_Track1_dz",  &reg_Track1_dz );
  outTree->Branch("reg_Track1_dEtaMu1", &reg_Track1_dEtaMu1 );
  outTree->Branch("reg_Track1_dPhiMu1", &reg_Track1_dPhiMu1 );
  outTree->Branch("reg_Track1_dEtaMu2", &reg_Track1_dEtaMu2 );
  outTree->Branch("reg_Track1_dPhiMu2", &reg_Track1_dPhiMu2 );
  outTree->Branch("reg_Track2_dxy", &reg_Track2_dxy );
  outTree->Branch("reg_Track2_dz",  &reg_Track2_dz );
  outTree->Branch("reg_Track2_dEtaMu1", &reg_Track2_dEtaMu1 );
  outTree->Branch("reg_Track2_dPhiMu1", &reg_Track2_dPhiMu1 );
  outTree->Branch("reg_Track2_dEtaMu2", &reg_Track2_dEtaMu2 );
  outTree->Branch("reg_Track2_dPhiMu2", &reg_Track2_dPhiMu2 );

  outTree->Branch("reg_Track1_nHits", &reg_Track1_nHits );
  outTree->Branch("reg_Track1_normChi2", &reg_Track1_normChi2 );
  outTree->Branch("reg_Track2_nHits", &reg_Track2_nHits );
  outTree->Branch("reg_Track2_normChi2", &reg_Track2_normChi2 );

  outTree->Branch("reg_collSize", &reg_collSize );
  outTree->Branch("regionalFrac_collSize", &regionalFrac_collSize );

  outTree->Branch("ctf0_Track1_pT", &ctf0_Track1_pT );
  outTree->Branch("ctf0_Track1_eta", &ctf0_Track1_eta );
  outTree->Branch("ctf0_Track1_phi", &ctf0_Track1_phi );
  outTree->Branch("ctf0_Track1_d0", &ctf0_Track1_d0 );
  outTree->Branch("ctf0_Track2_pT", &ctf0_Track2_pT );
  outTree->Branch("ctf0_Track2_eta", &ctf0_Track2_eta );
  outTree->Branch("ctf0_Track2_phi", &ctf0_Track2_phi );
  outTree->Branch("ctf0_Track2_d0", &ctf0_Track2_d0 );

  outTree->Branch("ctf0_Track1_deltaR", &ctf0_Track1_deltaR );
  outTree->Branch("ctf0_Track2_deltaR", &ctf0_Track2_deltaR );

  outTree->Branch("ctf0_Track1_dxy", &ctf0_Track1_dxy );
  outTree->Branch("ctf0_Track1_dz",  &ctf0_Track1_dz );
  outTree->Branch("ctf0_Track1_dEtaMu1", &ctf0_Track1_dEtaMu1 );
  outTree->Branch("ctf0_Track1_dPhiMu1", &ctf0_Track1_dPhiMu1 );
  outTree->Branch("ctf0_Track1_dEtaMu2", &ctf0_Track1_dEtaMu2 );
  outTree->Branch("ctf0_Track1_dPhiMu2", &ctf0_Track1_dPhiMu2 );
  outTree->Branch("ctf0_Track2_dxy", &ctf0_Track2_dxy );
  outTree->Branch("ctf0_Track2_dz",  &ctf0_Track2_dz );
  outTree->Branch("ctf0_Track2_dEtaMu1", &ctf0_Track2_dEtaMu1 );
  outTree->Branch("ctf0_Track2_dPhiMu1", &ctf0_Track2_dPhiMu1 );
  outTree->Branch("ctf0_Track2_dEtaMu2", &ctf0_Track2_dEtaMu2 );
  outTree->Branch("ctf0_Track2_dPhiMu2", &ctf0_Track2_dPhiMu2 );

  outTree->Branch("ctf0_Track1_nHits", &ctf0_Track1_nHits );
  outTree->Branch("ctf0_Track1_normChi2", &ctf0_Track1_normChi2 );
  outTree->Branch("ctf0_Track2_nHits", &ctf0_Track2_nHits );
  outTree->Branch("ctf0_Track2_normChi2", &ctf0_Track2_normChi2 );

  outTree->Branch("ctf0_collSize", &ctf0_collSize );

  outTree->Branch("iter0_Track1_pT", &iter0_Track1_pT );
  outTree->Branch("iter0_Track1_eta", &iter0_Track1_eta );
  outTree->Branch("iter0_Track1_phi", &iter0_Track1_phi );
  outTree->Branch("iter0_Track1_d0", &iter0_Track1_d0 );
  outTree->Branch("iter0_Track2_pT", &iter0_Track2_pT );
  outTree->Branch("iter0_Track2_eta", &iter0_Track2_eta );
  outTree->Branch("iter0_Track2_phi", &iter0_Track2_phi );
  outTree->Branch("iter0_Track2_d0", &iter0_Track2_d0 );

  outTree->Branch("iter0_collSize", &iter0_collSize );

  outTree->Branch("ctf1_Track1_pT", &ctf1_Track1_pT );
  outTree->Branch("ctf1_Track1_eta", &ctf1_Track1_eta );
  outTree->Branch("ctf1_Track1_phi", &ctf1_Track1_phi );
  outTree->Branch("ctf1_Track1_d0", &ctf1_Track1_d0 );
  outTree->Branch("ctf1_Track2_pT", &ctf1_Track2_pT );
  outTree->Branch("ctf1_Track2_eta", &ctf1_Track2_eta );
  outTree->Branch("ctf1_Track2_phi", &ctf1_Track2_phi );
  outTree->Branch("ctf1_Track2_d0", &ctf1_Track2_d0 );

  outTree->Branch("ctf1_Track1_deltaR", &ctf1_Track1_deltaR );
  outTree->Branch("ctf1_Track2_deltaR", &ctf1_Track2_deltaR );

  outTree->Branch("ctf1_Track1_dxy", &ctf1_Track1_dxy );
  outTree->Branch("ctf1_Track1_dz",  &ctf1_Track1_dz );
  outTree->Branch("ctf1_Track1_dEtaMu1", &ctf1_Track1_dEtaMu1 );
  outTree->Branch("ctf1_Track1_dPhiMu1", &ctf1_Track1_dPhiMu1 );
  outTree->Branch("ctf1_Track1_dEtaMu2", &ctf1_Track1_dEtaMu2 );
  outTree->Branch("ctf1_Track1_dPhiMu2", &ctf1_Track1_dPhiMu2 );
  outTree->Branch("ctf1_Track2_dxy", &ctf1_Track2_dxy );
  outTree->Branch("ctf1_Track2_dz",  &ctf1_Track2_dz );
  outTree->Branch("ctf1_Track2_dEtaMu1", &ctf1_Track2_dEtaMu1 );
  outTree->Branch("ctf1_Track2_dPhiMu1", &ctf1_Track2_dPhiMu1 );
  outTree->Branch("ctf1_Track2_dEtaMu2", &ctf1_Track2_dEtaMu2 );
  outTree->Branch("ctf1_Track2_dPhiMu2", &ctf1_Track2_dPhiMu2 );

  outTree->Branch("ctf1_Track1_nHits", &ctf1_Track1_nHits );
  outTree->Branch("ctf1_Track1_normChi2", &ctf1_Track1_normChi2 );
  outTree->Branch("ctf1_Track2_nHits", &ctf1_Track2_nHits );
  outTree->Branch("ctf1_Track2_normChi2", &ctf1_Track2_normChi2 );

  outTree->Branch("ctf1_collSize", &ctf1_collSize );

  outTree->Branch("iter1_Track1_pT", &iter1_Track1_pT );
  outTree->Branch("iter1_Track1_eta", &iter1_Track1_eta );
  outTree->Branch("iter1_Track1_phi", &iter1_Track1_phi );
  outTree->Branch("iter1_Track1_d0", &iter1_Track1_d0 );
  outTree->Branch("iter1_Track2_pT", &iter1_Track2_pT );
  outTree->Branch("iter1_Track2_eta", &iter1_Track2_eta );
  outTree->Branch("iter1_Track2_phi", &iter1_Track2_phi );
  outTree->Branch("iter1_Track2_d0", &iter1_Track2_d0 );

  outTree->Branch("iter1_collSize", &iter1_collSize );

  outTree->Branch("ctf2_Track1_pT", &ctf2_Track1_pT );
  outTree->Branch("ctf2_Track1_eta", &ctf2_Track1_eta );
  outTree->Branch("ctf2_Track1_phi", &ctf2_Track1_phi );
  outTree->Branch("ctf2_Track1_d0", &ctf2_Track1_d0 );
  outTree->Branch("ctf2_Track2_pT", &ctf2_Track2_pT );
  outTree->Branch("ctf2_Track2_eta", &ctf2_Track2_eta );
  outTree->Branch("ctf2_Track2_phi", &ctf2_Track2_phi );
  outTree->Branch("ctf2_Track2_d0", &ctf2_Track2_d0 );

  outTree->Branch("ctf2_Track1_deltaR", &ctf2_Track1_deltaR );
  outTree->Branch("ctf2_Track2_deltaR", &ctf2_Track2_deltaR );

  outTree->Branch("ctf2_Track1_dxy", &ctf2_Track1_dxy );
  outTree->Branch("ctf2_Track1_dz",  &ctf2_Track1_dz );
  outTree->Branch("ctf2_Track1_dEtaMu1", &ctf2_Track1_dEtaMu1 );
  outTree->Branch("ctf2_Track1_dPhiMu1", &ctf2_Track1_dPhiMu1 );
  outTree->Branch("ctf2_Track1_dEtaMu2", &ctf2_Track1_dEtaMu2 );
  outTree->Branch("ctf2_Track1_dPhiMu2", &ctf2_Track1_dPhiMu2 );
  outTree->Branch("ctf2_Track2_dxy", &ctf2_Track2_dxy );
  outTree->Branch("ctf2_Track2_dz",  &ctf2_Track2_dz );
  outTree->Branch("ctf2_Track2_dEtaMu1", &ctf2_Track2_dEtaMu1 );
  outTree->Branch("ctf2_Track2_dPhiMu1", &ctf2_Track2_dPhiMu1 );
  outTree->Branch("ctf2_Track2_dEtaMu2", &ctf2_Track2_dEtaMu2 );
  outTree->Branch("ctf2_Track2_dPhiMu2", &ctf2_Track2_dPhiMu2 );

  outTree->Branch("ctf2_Track1_nHits", &ctf2_Track1_nHits );
  outTree->Branch("ctf2_Track1_normChi2", &ctf2_Track1_normChi2 );
  outTree->Branch("ctf2_Track2_nHits", &ctf2_Track2_nHits );
  outTree->Branch("ctf2_Track2_normChi2", &ctf2_Track2_normChi2 );

  outTree->Branch("ctf2_collSize", &ctf2_collSize );

  outTree->Branch("iter2_Track1_pT", &iter2_Track1_pT );
  outTree->Branch("iter2_Track1_eta", &iter2_Track1_eta );
  outTree->Branch("iter2_Track1_phi", &iter2_Track1_phi );
  outTree->Branch("iter2_Track1_d0", &iter2_Track1_d0 );
  outTree->Branch("iter2_Track1_d0err", &iter2_Track1_d0err );
  outTree->Branch("iter2_Track2_pT", &iter2_Track2_pT );
  outTree->Branch("iter2_Track2_eta", &iter2_Track2_eta );
  outTree->Branch("iter2_Track2_phi", &iter2_Track2_phi );
  outTree->Branch("iter2_Track2_d0", &iter2_Track2_d0 );
  outTree->Branch("iter2_Track2_d0err", &iter2_Track2_d0err );

  outTree->Branch("iter2_collSize", &iter2_collSize );

  outTree->Branch("prefilter_newTks_pT",    &prefilter_newTks_pT );
  outTree->Branch("prefilter_newTks_eta",   &prefilter_newTks_eta );
  outTree->Branch("prefilter_newTks_phi",   &prefilter_newTks_phi );
  outTree->Branch("prefilter_newTks_dEta1", &prefilter_newTks_dEta1 );
  outTree->Branch("prefilter_newTks_dPhi1", &prefilter_newTks_dPhi1 );
  outTree->Branch("prefilter_newTks_dEta2", &prefilter_newTks_dEta2 );
  outTree->Branch("prefilter_newTks_dPhi2", &prefilter_newTks_dPhi2 );
  outTree->Branch("prefilter_newTks_dxy",   &prefilter_newTks_dxy );
  outTree->Branch("prefilter_newTks_dz",    &prefilter_newTks_dz );
  outTree->Branch("prefilter_newTks_d0",    &prefilter_newTks_d0 );
  outTree->Branch("prefilter_newTks_d0err", &prefilter_newTks_d0err );
  outTree->Branch("prefilter_newTks_bool",  &prefilter_newTks_bool );

  outTree->Branch("dummy", &dummy);

  // muonFilterMap = new std::map<std::string, std::string>();

  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_3_Bs_v","hltDisplacedmumuFilterDoubleMu4Bs"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_3_Jpsi_Displaced_v","hltDisplacedmumuFilterDoubleMu43Jpsi"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4LowMassNonResonant"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_v","hltDisplacedmumuFilterDoubleMu3Tau3mu"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4PsiPrime"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Jpsi_v","hltL3fSQMu7p5L2Mu2L3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Upsilon_v","hltL3fSQMu7p5L2Mu2L3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Onia_v","hltL3fL1sMu22orMu20erorMu25L1f0L2f0L3Filtered25"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu30_TkMu0_Onia_v","hltL3fL1sMu22orMu20erorMu25L1f0L2f0L3Filtered30"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu20_TkMu0_Phi_v","hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Phi_v","hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered25"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_L1_NoOS_v","hltDisplacedmumuFilterDimuon0JpsiL1sNoOS"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v","hltDimuon0JpsiNoVtxNoOSL3Filtered"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_v","hltDisplacedmumuFilterDimuon0Jpsi"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_v","hltDimuon0JpsiL3Filtered"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v","hltDisplacedmumuFilterDimuon0JpsiL1s4R0er1p5R"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v","hltDimuon0JpsiL1s4R0er1p5RL3Filtered"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi3p5_Muon2_v","hltVertexmumuFilterJpsiMuon3p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_5_v","hltDisplacedmumuFilterDimuon0UpsilonL1s5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5NoOS_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5NoOS"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5er2p0_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5er2p0"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5er2p0M"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_NoVertexing_v","hltDimuon0UpsilonL1s4p5er2p0ML3Filtered"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_5M_v","hltDisplacedmumuFilterDimuon0UpsilonL1s5M"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_0er1p5R_v","hltDisplacedmumuFilterDimuon0LowMassL1s0er1p5R"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_0er1p5_v","hltDisplacedmumuFilterDimuon0LowMassL1s0er1p5"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_v","hltDisplacedmumuFilterDimuon0LowMass"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_4_v","hltDisplacedmumuFilterDimuon0LowMassL1s4"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_4R_v","hltDisplacedmumuFilterDimuon0LowMassL1s4R"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_TM530_v","hltDisplacedmumuFilterDimuon0LowMassL1sTM530"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_Muon_L1_TM0_v","hltVertexmumuFilterUpsilon0MuonL1sTM0"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v","hltVertexmumuFilterUpsilon0MuonNoL1Mass"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v","hltDisplacedmumuFilterDoubleMu3Tau3muNoL1Mass"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_Jpsi_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_Jpsi_NoVertexing_v","hltDoubleMu4JpsiDisplacedL3Filtered"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon10PsiPrimeBarrelnoCow"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon10_Upsilon_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon10UpsilonBarrelnoCow"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon12_Upsilon_eta1p5_v","hltDisplacedmumuFilterDimuon12Upsilons"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon14_Phi_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon14PhiBarrelnoCow"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon18_PsiPrime_v","hltDisplacedmumuFilterDimuon18PsiPrimes"));
  // muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon25_Jpsi_v","hltDisplacedmumuFilterDimuon25Jpsis"));

  // extraFilterMap = new std::map<std::string, std::string>();

  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Jpsi_v","hltSQMu7p5L2Mu2JpsiTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Upsilon_v","hltSQMu7p5L2Mu2UpsilonTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Jpsi_v","hltMu7p5Track2JpsiTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Jpsi_v","hltMu7p5Track3p5JpsiTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Jpsi_v","hltMu7p5Track7JpsiTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Upsilon_v","hltMu7p5Track2UpsilonTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Upsilon_v","hltMu7p5Track3p5UpsilonTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Upsilon_v","hltMu7p5Track7UpsilonTrackMassFiltered"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Onia_v","hltDiMuonGlbFiltered25TrkFiltered0"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu30_TkMu0_Onia_v","hltDiMuonGlbFiltered30TrkFiltered0"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu20_TkMu0_Phi_v","hltDiMuonGlbFiltered20TrkFiltered0"));
  // extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Phi_v","hltDiMuonGlbFiltered25PhiTrkFiltered0"));

  // trackFilterMap = new std::map<std::string, std::string>();

  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrk_Displaced_v","hltJpsiTkVertexFilter"));
  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v","hltLowMassNonResonantTkVertexFilter"));
  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_v","hltTau3muTkVertexFilter"));
  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v","hltPsiPrimeTkVertexFilter"));
  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v","hltTau3muNoL1MassTkVertexFilter"));
  // trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v","hltJpsiTkTkVertexFilterPhiKstar"));
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VanillaHLTAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VanillaHLTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VanillaHLTAnalyzer);
