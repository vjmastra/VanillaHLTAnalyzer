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

   passHLT_DoubleMu4_Jpsi_Displaced = 0;
   passHLT_DoubleMu4_JpsiTrk_Displaced = 0;

   muonsMatched = 0;

   hltMatchedLeadingTrack = false;
   hltMatchedSubleadingTrack = false;

   JpsiTrackMatched = 0;

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

   Mu1_deltaR = -10;
   Mu2_deltaR = -10;
   Track1_deltaR = -10;
   Track2_deltaR = -10;

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
     if (hltNames[i].find("HLT_DoubleMu4_Jpsi_Displaced_v") != std::string::npos)
       hltResults[i] == 1 ? passHLT_DoubleMu4_Jpsi_Displaced = 1 : passHLT_DoubleMu4_Jpsi_Displaced = 0;
     if (hltNames[i].find("HLT_DoubleMu4_JpsiTrk_Displaced_v") != std::string::npos)
       hltResults[i] == 1 ? passHLT_DoubleMu4_JpsiTrk_Displaced = 1 : passHLT_DoubleMu4_JpsiTrk_Displaced = 0;
   }    

   //Look for gen particles

   edm::Handle< std::vector<reco::GenParticle> >  genParticles;
   iEvent.getByToken(genParticlesToken_, genParticles);

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

     //Match with trigger object

     edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
     iEvent.getByToken(triggerObjectsToken_, triggerObjectsSummary);
     trigger::TriggerObjectCollection selectedObjects;

     bool matchMuon1 = false;
     bool matchMuon2 = false;
     bool matchTrack = false;

     std::vector<float> hltTriggerObjects_pT;
     std::vector<float> hltTriggerObjects_eta;
     std::vector<float> hltTriggerObjects_phi;
 
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

       float dRthreshold = 0.03;
       for (unsigned int i = 0; i < hltTriggerObjects_pT.size(); i++) {
         float dR1 = deltaR(Muon1_eta, Muon1_phi, hltTriggerObjects_eta[i], hltTriggerObjects_phi[i]);
         float dR2 = deltaR(Muon2_eta, Muon2_phi, hltTriggerObjects_eta[i], hltTriggerObjects_phi[i]);
         if (dR1 < dRthreshold) {
           hlt_Mu1_pT = hltTriggerObjects_pT[i];
           hlt_Mu1_eta = hltTriggerObjects_eta[i];
           hlt_Mu1_phi = hltTriggerObjects_phi[i];
           Mu1_deltaR = dR1;
           matchMuon1 = true;
         }
         if (dR2 < dRthreshold) {
           hlt_Mu2_pT = hltTriggerObjects_pT[i];
           hlt_Mu2_eta = hltTriggerObjects_eta[i];
           hlt_Mu2_phi = hltTriggerObjects_phi[i];
           Mu2_deltaR = dR2;
           matchMuon2 = true;
         }
       }      
     }//end passHLT

     //al momento non èescluso che i muoni vengano da due diverse jpsi che hanno fatto scattare il trigger
     muonsMatched = matchMuon1 && matchMuon2;

     hltTriggerObjects_pT.clear();
     hltTriggerObjects_eta.clear();
     hltTriggerObjects_phi.clear();

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

       float dRthreshold = 0.03;
       for (unsigned int i = 2; i < hltTriggerObjects_pT.size(); i=i+3) {//loop on tracks only
         float dR1 = deltaR(track1_eta, track1_phi, hltTriggerObjects_eta[i], hltTriggerObjects_phi[i]);
         float dR2 = deltaR(track2_eta, track2_phi, hltTriggerObjects_eta[i], hltTriggerObjects_phi[i]);
         if (dR1 < dRthreshold) {
           hlt_Track1_pT = hltTriggerObjects_pT[i];
           hlt_Track1_eta = hltTriggerObjects_eta[i];
           hlt_Track1_phi = hltTriggerObjects_phi[i];
           Track1_deltaR = dR1;
           hltMatchedLeadingTrack = true;
         }
         if (dR2 < dRthreshold) {
           hlt_Track2_pT = hltTriggerObjects_pT[i];
           hlt_Track2_eta = hltTriggerObjects_eta[i];
           hlt_Track2_phi = hltTriggerObjects_phi[i];
           Track2_deltaR = dR2;
           hltMatchedSubleadingTrack = true;
         }
       } 

       matchTrack = hltMatchedLeadingTrack || hltMatchedSubleadingTrack;
     }//end passHLT

     JpsiTrackMatched = muonsMatched && matchTrack;

     outTree->Fill();

   }//end Bs for

}


float
VanillaHLTAnalyzer::deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dPhi = fabs(phi1 - phi2);
  if ( dPhi > TMath::Pi() )  dPhi = 2*TMath::Pi() - dPhi;

  return sqrt( (eta1-eta2)*(eta1-eta2) + dPhi*dPhi );
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

  outTree->Branch("passHLT_DoubleMu4_Jpsi_Displaced", &passHLT_DoubleMu4_Jpsi_Displaced );
  outTree->Branch("passHLT_DoubleMu4_JpsiTrk_Displaced", &passHLT_DoubleMu4_JpsiTrk_Displaced );

  outTree->Branch("muonsMatched", &muonsMatched );

  outTree->Branch("hltMatchedLeadingTrack", &hltMatchedLeadingTrack );
  outTree->Branch("hltMatchedSubleadingTrack", &hltMatchedSubleadingTrack );

  outTree->Branch("JpsiTrackMatched", &JpsiTrackMatched );

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

  outTree->Branch("Mu1_deltaR", &Mu1_deltaR );
  outTree->Branch("Mu2_deltaR", &Mu2_deltaR );
  outTree->Branch("Track1_deltaR", &Track1_deltaR );
  outTree->Branch("Track2_deltaR", &Track2_deltaR );

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
