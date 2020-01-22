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

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

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

      FreeTrajectoryState initialFreeState( const reco::Track& tk, const MagneticField* field);


      // ----------member data ---------------------------

      // Trigger process
      edm::InputTag triggerResultTag_;
      edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
      edm::InputTag triggerTag_;
      edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > triggerToken_;

      edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;

      /* edm::InputTag l1candTag_; */
      /* edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; */
      edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_;

      edm::InputTag offlinePVTag_;
      edm::EDGetTokenT<reco::VertexCollection> offlinePVToken_;
      edm::InputTag beamspotTag_;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

      edm::InputTag offlineMuonsTag_;
      edm::EDGetTokenT< std::vector< pat::Muon > > offlineMuonsToken_;
      edm::InputTag offlineTksTag_;
      edm::EDGetTokenT<reco::TrackCollection> offlineTksToken_;

      /* HLTPrescaleProvider* prescaleProvider; */

      //tree and content

      edm::Service<TFileService> outfile_;

      TTree* outTree;

      std::vector<std::string> hltNames;
      std::vector<bool> hltResults;
      std::vector<bool> hltMatchDimu;
      std::vector<bool> hltMatchB;
      std::vector<int> hltPrescales;
      std::vector<std::string> l1tNames;
      std::vector<int> l1tPrescales;

      bool doBp = true;

      float Dimuon_CL;
      float Dimuon_LS;
      float Dimuon_CosAlpha;
      float Dimuon_Mass;
      float Dimuon_pT;
      float Dimuon_eta;
      float Dimuon_phi;

      float Muon1_pT;
      float Muon1_eta;
      float Muon1_phi;
      float Muon1_cha;
      float Muon2_pT;
      float Muon2_eta;
      float Muon2_phi;
      float Muon2_cha;

      float Bp_CL;
      float Bp_LS;
      float Bp_CosAlpha;
      float Bp_Mass;
      float Bp_pT;
      float Bp_eta;
      float Bp_phi;

      float track_pT;
      float track_eta;
      float track_phi;
      float track_cha;
      float track_d0;

      float nOffVtx;

      std::map < std::string, std::string > * muonFilterMap;
      std::map < std::string, std::string > * extraFilterMap;
      std::map < std::string, std::string > * trackFilterMap;

      l1t::L1TGlobalUtil *fGtUtil;
      HLTPrescaleProvider hltPrescale_;
      /* bool changedFlag  = true; */
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
