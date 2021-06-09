/** \class VanillaHLTAnalyzer
 */

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

#include <memory>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class VanillaHLTAnalyzer : public edm::one::EDAnalyzer<>  {
   public:
      explicit VanillaHLTAnalyzer(const edm::ParameterSet&);
      ~VanillaHLTAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef std::map<uint, reco::Track > selTracksDef ;

      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

   private:

      float deltaR(float eta1, float phi1, float eta2, float phi2);

      float deltaPtOverPt(float pt1, float pt2);

      void printProgeny( const reco::GenParticle& part);

      FreeTrajectoryState initialFreeState( const reco::Track& tk, const MagneticField* field);


      // ----------member data ---------------------------

      // Trigger process
      edm::InputTag triggerResultTag_;
      edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;

      /* edm::InputTag triggerTag_; */
      /* edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > triggerToken_; */
      /* edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_; */

      /* edm::InputTag l1candTag_; */
      /* edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; */
      /* edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_; */

      edm::InputTag genParticlesTag_;
      edm::EDGetTokenT< std::vector< reco::GenParticle > > genParticlesToken_;

      edm::InputTag triggerObjectsTag_;
      edm::EDGetTokenT<trigger::TriggerEvent>  triggerObjectsToken_;

      edm::InputTag jpsiVerticesTag_;
      edm::EDGetTokenT<vector<reco::Vertex>> jpsiVerticesToken_;

      edm::InputTag legacyPixelTracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacyPixelTracksToken_;

      edm::InputTag legacyIter0TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacyIter0TracksToken_;

      edm::InputTag legacyIter1TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacyIter1TracksToken_;

      edm::InputTag legacyIter2TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacyIter2TracksToken_;

      edm::InputTag legacy3RecoTracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacy3RecoTracksToken_;

      edm::InputTag legacy2RecoTracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> legacy2RecoTracksToken_;

      edm::InputTag myPixelTracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> myPixelTracksToken_;

      edm::InputTag regionalPixelTracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> regionalPixelTracksToken_;

      edm::InputTag ctfIter0TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> ctfIter0TracksToken_;

      edm::InputTag triggerIter0TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> triggerIter0TracksToken_;

      edm::InputTag ctfIter1TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> ctfIter1TracksToken_;

      edm::InputTag triggerIter1TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> triggerIter1TracksToken_;

      edm::InputTag ctfIter2TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> ctfIter2TracksToken_;

      edm::InputTag triggerIter2TracksTag_;
      edm::EDGetTokenT<vector<reco::Track>> triggerIter2TracksToken_;

      edm::InputTag beamSpotTag_;
      edm::EDGetTokenT <reco::BeamSpot> beamSpotToken_;

      std::string hltTag_;
      /* HLTPrescaleProvider* prescaleProvider; */

      //tree and content

      edm::Service<TFileService> outfile_;

      TTree* outTree;

      std::vector<std::string> hltNames;
      std::vector<bool> hltResults;

      float Muon1_pT;
      float Muon1_eta;
      float Muon1_phi;
      float Muon1_cha;
      float Muon2_pT;
      float Muon2_eta;
      float Muon2_phi;
      float Muon2_cha;

      float track1_pT;
      float track1_eta;
      float track1_phi;
      float track1_cha;
      float track2_pT;
      float track2_eta;
      float track2_phi;
      float track2_cha;

      float bsFlightLength2D;
      float bsFlightLength3D;

      int passHLT_DoubleMu4_Jpsi_Displaced;
      int passHLT_DoubleMu4_JpsiTrk_Displaced;
      int passHLT_DoubleMu4_JpsiTrk_Displaced_GPU;

      int muonsMatched;

      bool hltMatch_Track1;
      bool hltMatch_Track2;
      bool gpuMatch_Track1;
      bool gpuMatch_Track2;

      bool legPixMatch_Track1;
      bool legPixMatch_Track2;
      bool legI0Match_Track1;
      bool legI0Match_Track2;
      bool legI1Match_Track1;
      bool legI1Match_Track2;
      bool legI2Match_Track1;
      bool legI2Match_Track2;
      bool leg3RMatch_Track1;
      bool leg3RMatch_Track2;
      bool leg2RMatch_Track1;
      bool leg2RMatch_Track2;

      bool myPixMatch_Track1;
      bool myPixMatch_Track2;
      bool regMatch_Track1;
      bool regMatch_Track2;
      bool ctf0Match_Track1;
      bool ctf0Match_Track2;
      bool iter0Match_Track1;
      bool iter0Match_Track2;
      bool ctf1Match_Track1;
      bool ctf1Match_Track2;
      bool iter1Match_Track1;
      bool iter1Match_Track2;
      bool ctf2Match_Track1;
      bool ctf2Match_Track2;
      bool iter2Match_Track1;
      bool iter2Match_Track2;

      int hltMatch_Track;
      int gpuMatch_Track;

      std::vector<float> hltTriggerTrack_pt;
      std::vector<float> hltTriggerTrack_eta;
      std::vector<float> hltTriggerFake_pt;
      std::vector<float> hltTriggerFake_eta;
      
      std::vector<float> gpuTriggerTrack_pt;
      std::vector<float> gpuTriggerTrack_eta;
      std::vector<float> gpuTriggerFake_pt;
      std::vector<float> gpuTriggerFake_eta;

      float j_vtx_x;
      float j_vtx_y;
      float j_vtx_z;

      unsigned int j_vtx_size;

      float hlt_Mu1_pT;
      float hlt_Mu1_eta;
      float hlt_Mu1_phi;
      float hlt_Mu1_cha;
      float hlt_Mu2_pT;
      float hlt_Mu2_eta;
      float hlt_Mu2_phi;
      float hlt_Mu2_cha; 

      float hlt_Track1_pT;
      float hlt_Track1_eta;
      float hlt_Track1_phi;   
      float hlt_Track2_pT;
      float hlt_Track2_eta;
      float hlt_Track2_phi; 

      float hlt_Track1_deltaR;
      float hlt_Track2_deltaR;

      float gpu_Track1_pT;
      float gpu_Track1_eta;
      float gpu_Track1_phi;
      float gpu_Track2_pT;
      float gpu_Track2_eta;
      float gpu_Track2_phi;

      float gpu_Track1_deltaR;
      float gpu_Track2_deltaR;

      float Mu1_deltaR;
      float Mu2_deltaR;

      float legPix_Track1_pT;
      float legPix_Track1_eta;
      float legPix_Track1_phi;
      float legPix_Track1_d0;
      float legPix_Track2_pT;
      float legPix_Track2_eta;
      float legPix_Track2_phi;
      float legPix_Track2_d0;

      float legPix_Track1_deltaR;
      float legPix_Track2_deltaR;

      float legPix_Track1_dxy;
      float legPix_Track1_dz;
      float legPix_Track1_dEtaMu1;
      float legPix_Track1_dPhiMu1;
      float legPix_Track1_dEtaMu2;
      float legPix_Track1_dPhiMu2;
      float legPix_Track2_dxy;
      float legPix_Track2_dz;
      float legPix_Track2_dEtaMu1;
      float legPix_Track2_dPhiMu1;
      float legPix_Track2_dEtaMu2;
      float legPix_Track2_dPhiMu2;

      unsigned int legPix_Track1_nHits;
      float        legPix_Track1_normChi2;
      unsigned int legPix_Track2_nHits;
      float        legPix_Track2_normChi2;

      unsigned int legPix_collSize;

      float legI0_Track1_pT;
      float legI0_Track1_eta;
      float legI0_Track1_phi;
      float legI0_Track1_d0;
      float legI0_Track2_pT;
      float legI0_Track2_eta;
      float legI0_Track2_phi;
      float legI0_Track2_d0;

      unsigned int legI0_collSize;

      float legI1_Track1_pT;
      float legI1_Track1_eta;
      float legI1_Track1_phi;
      float legI1_Track1_d0;
      float legI1_Track2_pT;
      float legI1_Track2_eta;
      float legI1_Track2_phi;
      float legI1_Track2_d0;

      unsigned int legI1_collSize;

      float legI2_Track1_pT;
      float legI2_Track1_eta;
      float legI2_Track1_phi;
      float legI2_Track1_d0;
      float legI2_Track2_pT;
      float legI2_Track2_eta;
      float legI2_Track2_phi;
      float legI2_Track2_d0;

      unsigned int legI2_collSize;

      float leg3R_Track1_pT;
      float leg3R_Track1_eta;
      float leg3R_Track1_phi;
      float leg3R_Track1_d0;
      float leg3R_Track2_pT;
      float leg3R_Track2_eta;
      float leg3R_Track2_phi;
      float leg3R_Track2_d0;

      float leg2R_Track1_pT;
      float leg2R_Track1_eta;
      float leg2R_Track1_phi;
      float leg2R_Track1_d0;
      float leg2R_Track1_d0err;
      float leg2R_Track2_pT;
      float leg2R_Track2_eta;
      float leg2R_Track2_phi;
      float leg2R_Track2_d0;
      float leg2R_Track2_d0err;
 
      std::vector<float> prefilter_legTks_pT;
      std::vector<float> prefilter_legTks_eta;
      std::vector<float> prefilter_legTks_phi;
      std::vector<float> prefilter_legTks_dEta1;
      std::vector<float> prefilter_legTks_dPhi1;
      std::vector<float> prefilter_legTks_dEta2;
      std::vector<float> prefilter_legTks_dPhi2;
      std::vector<float> prefilter_legTks_dxy;
      std::vector<float> prefilter_legTks_dz;
      std::vector<float> prefilter_legTks_d0;
      std::vector<float> prefilter_legTks_d0err;
      std::vector<bool>  prefilter_legTks_bool;

      float myPix_Track1_pT;
      float myPix_Track1_eta;
      float myPix_Track1_phi;
      float myPix_Track1_d0;
      float myPix_Track2_pT;
      float myPix_Track2_eta;
      float myPix_Track2_phi;
      float myPix_Track2_d0;

      float myPix_Track1_deltaR;
      float myPix_Track2_deltaR;

      float myPix_Track1_dxy;
      float myPix_Track1_dz;
      float myPix_Track1_dEtaMu1;
      float myPix_Track1_dPhiMu1;
      float myPix_Track1_dEtaMu2;
      float myPix_Track1_dPhiMu2;
      float myPix_Track2_dxy;
      float myPix_Track2_dz;
      float myPix_Track2_dEtaMu1;
      float myPix_Track2_dPhiMu1;
      float myPix_Track2_dEtaMu2;
      float myPix_Track2_dPhiMu2;

      unsigned int myPix_Track1_nHits;
      float myPix_Track1_normChi2;
      unsigned int myPix_Track2_nHits;
      float myPix_Track2_normChi2;

      unsigned int myPix_collSize;

      float reg_Track1_pT;
      float reg_Track1_eta;
      float reg_Track1_phi;
      float reg_Track1_d0;
      float reg_Track2_pT;
      float reg_Track2_eta;
      float reg_Track2_phi;
      float reg_Track2_d0;

      float reg_Track1_deltaR;
      float reg_Track2_deltaR;

      float reg_Track1_dxy;
      float reg_Track1_dz;
      float reg_Track1_dEtaMu1;
      float reg_Track1_dPhiMu1;
      float reg_Track1_dEtaMu2;
      float reg_Track1_dPhiMu2;
      float reg_Track2_dxy;
      float reg_Track2_dz;
      float reg_Track2_dEtaMu1;
      float reg_Track2_dPhiMu1;
      float reg_Track2_dEtaMu2;
      float reg_Track2_dPhiMu2;

      unsigned int reg_Track1_nHits;
      float        reg_Track1_normChi2;
      unsigned int reg_Track2_nHits;
      float        reg_Track2_normChi2;

      unsigned int reg_collSize;
      float        regionalFrac_collSize;

      float ctf0_Track1_pT;
      float ctf0_Track1_eta;
      float ctf0_Track1_phi;
      float ctf0_Track1_d0;
      float ctf0_Track2_pT;
      float ctf0_Track2_eta;
      float ctf0_Track2_phi;
      float ctf0_Track2_d0;

      float ctf0_Track1_deltaR;
      float ctf0_Track2_deltaR;

      float ctf0_Track1_dxy;
      float ctf0_Track1_dz;
      float ctf0_Track1_dEtaMu1;
      float ctf0_Track1_dPhiMu1;
      float ctf0_Track1_dEtaMu2;
      float ctf0_Track1_dPhiMu2;
      float ctf0_Track2_dxy;
      float ctf0_Track2_dz;
      float ctf0_Track2_dEtaMu1;
      float ctf0_Track2_dPhiMu1;
      float ctf0_Track2_dEtaMu2;
      float ctf0_Track2_dPhiMu2;

      unsigned int ctf0_Track1_nHits;
      float        ctf0_Track1_normChi2;
      unsigned int ctf0_Track2_nHits;
      float        ctf0_Track2_normChi2;

      unsigned int ctf0_collSize;

      float iter0_Track1_pT;
      float iter0_Track1_eta;
      float iter0_Track1_phi;
      float iter0_Track1_d0;
      float iter0_Track2_pT;
      float iter0_Track2_eta;
      float iter0_Track2_phi;
      float iter0_Track2_d0;

      unsigned int iter0_collSize;

      float ctf1_Track1_pT;
      float ctf1_Track1_eta;
      float ctf1_Track1_phi;
      float ctf1_Track1_d0;
      float ctf1_Track2_pT;
      float ctf1_Track2_eta;
      float ctf1_Track2_phi;
      float ctf1_Track2_d0;

      float ctf1_Track1_deltaR;
      float ctf1_Track2_deltaR;

      float ctf1_Track1_dxy;
      float ctf1_Track1_dz;
      float ctf1_Track1_dEtaMu1;
      float ctf1_Track1_dPhiMu1;
      float ctf1_Track1_dEtaMu2;
      float ctf1_Track1_dPhiMu2;
      float ctf1_Track2_dxy;
      float ctf1_Track2_dz;
      float ctf1_Track2_dEtaMu1;
      float ctf1_Track2_dPhiMu1;
      float ctf1_Track2_dEtaMu2;
      float ctf1_Track2_dPhiMu2;

      unsigned int ctf1_Track1_nHits;
      float        ctf1_Track1_normChi2;
      unsigned int ctf1_Track2_nHits;
      float        ctf1_Track2_normChi2;

      unsigned int ctf1_collSize;

      float iter1_Track1_pT;
      float iter1_Track1_eta;
      float iter1_Track1_phi;
      float iter1_Track1_d0;
      float iter1_Track2_pT;
      float iter1_Track2_eta;
      float iter1_Track2_phi;
      float iter1_Track2_d0;

      unsigned int iter1_collSize;

      float ctf2_Track1_pT;
      float ctf2_Track1_eta;
      float ctf2_Track1_phi;
      float ctf2_Track1_d0;
      float ctf2_Track2_pT;
      float ctf2_Track2_eta;
      float ctf2_Track2_phi;
      float ctf2_Track2_d0;

      float ctf2_Track1_deltaR;
      float ctf2_Track2_deltaR;

      float ctf2_Track1_dxy;
      float ctf2_Track1_dz;
      float ctf2_Track1_dEtaMu1;
      float ctf2_Track1_dPhiMu1;
      float ctf2_Track1_dEtaMu2;
      float ctf2_Track1_dPhiMu2;
      float ctf2_Track2_dxy;
      float ctf2_Track2_dz;
      float ctf2_Track2_dEtaMu1;
      float ctf2_Track2_dPhiMu1;
      float ctf2_Track2_dEtaMu2;
      float ctf2_Track2_dPhiMu2;

      unsigned int ctf2_Track1_nHits;
      float        ctf2_Track1_normChi2;
      unsigned int ctf2_Track2_nHits;
      float        ctf2_Track2_normChi2;

      unsigned int ctf2_collSize;

      float iter2_Track1_pT;
      float iter2_Track1_eta;
      float iter2_Track1_phi;
      float iter2_Track1_d0;
      float iter2_Track1_d0err;
      float iter2_Track2_pT;
      float iter2_Track2_eta;
      float iter2_Track2_phi;
      float iter2_Track2_d0;
      float iter2_Track2_d0err; 

      unsigned int iter2_collSize;

      std::vector<float> prefilter_newTks_pT;
      std::vector<float> prefilter_newTks_eta;
      std::vector<float> prefilter_newTks_phi;
      std::vector<float> prefilter_newTks_dEta1;
      std::vector<float> prefilter_newTks_dPhi1;
      std::vector<float> prefilter_newTks_dEta2;
      std::vector<float> prefilter_newTks_dPhi2;
      std::vector<float> prefilter_newTks_dxy;
      std::vector<float> prefilter_newTks_dz;
      std::vector<float> prefilter_newTks_d0;
      std::vector<float> prefilter_newTks_d0err;
      std::vector<bool>  prefilter_newTks_bool;

      int dummy;

      /* std::map < std::string, std::string > * muonFilterMap; */
      /* std::map < std::string, std::string > * extraFilterMap; */
      /* std::map < std::string, std::string > * trackFilterMap; */

      /* l1t::L1TGlobalUtil *fGtUtil; */
      /* HLTPrescaleProvider hltPrescale_; */
      /* bool changedFlag  = true; */
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
