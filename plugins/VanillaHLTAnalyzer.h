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
