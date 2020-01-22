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
  triggerTag_         (iConfig.getUntrackedParameter<edm::InputTag>("slimmedPatTrigger")), 
  triggerToken_       (consumes< std::vector<pat::TriggerObjectStandAlone> >(triggerTag_)),

  triggerPrescales_ (consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),

  // l1candTag_              (iConfig.getUntrackedParameter<edm::InputTag>("L1Candidates")),
  // l1candToken_            (consumes<l1t::MuonBxCollection>(l1candTag_)),
  l1results_      (consumes<GlobalAlgBlkBxCollection>  (iConfig.getParameter<edm::InputTag>("l1results"))),

  offlinePVTag_           (iConfig.getParameter<edm::InputTag>("offlineVtx")), 
  offlinePVToken_         (consumes<reco::VertexCollection>(offlinePVTag_)), 
  beamspotTag_            (iConfig.getParameter<edm::InputTag>("beamspot")), 
  beamspotToken_          (consumes<reco::BeamSpot>(beamspotTag_)), 

  offlineMuonsTag_        (iConfig.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag")),
  offlineMuonsToken_      (consumes< std::vector<pat::Muon> >(offlineMuonsTag_)), 
  offlineTksTag_          (iConfig.getUntrackedParameter<edm::InputTag>("OfflineTkTag")),
  offlineTksToken_        (consumes<reco::TrackCollection>(offlineTksTag_)),

  hltPrescale_ (iConfig, consumesCollector(), *this)
  // hltPrescale_ (new HLTPrescaleProvider(iConfig, consumesCollector(), *this))

{
  //now do what ever initialization is needed
  // usesResource("TFileService");
  fGtUtil = new l1t::L1TGlobalUtil(iConfig, consumesCollector(), *this, iConfig.getParameter<edm::InputTag>("l1results"), iConfig.getParameter<edm::InputTag>("l1results"));

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

   // hltPrescale_->init(iEvent.getRun(),iSetup,"HLT",changedFlag);

   hltNames.clear();
   hltResults.clear();
   hltMatchDimu.clear();
   hltMatchB.clear();
   hltPrescales.clear();
   l1tNames.clear();
   l1tPrescales.clear();

   Dimuon_CL = 0;
   Dimuon_LS = 0;
   Dimuon_CosAlpha = 0;
   Dimuon_Mass = 0;
   Dimuon_pT = 0;
   Dimuon_eta = 0;
   Dimuon_phi = 0;

   Muon1_pT = 0;
   Muon1_eta = 0;
   Muon1_phi = 0;
   Muon1_cha = 0;
   Muon2_pT = 0;
   Muon2_eta = 0;
   Muon2_phi = 0;
   Muon2_cha = 0;

   Bp_CL       = 0;
   Bp_LS       = 0;
   Bp_CosAlpha = 0;
   Bp_Mass     = 0;
   Bp_pT       = 0;
   Bp_eta      = 0;
   Bp_phi      = 0;
       
   track_pT  = 0;
   track_eta = 0;
   track_phi = 0;
   track_d0  = 0;
   track_cha = 0;

   nOffVtx = 0;

   edm::Handle<edm::TriggerResults> triggerResults;
   edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerEvent;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH;

   if ( !(iEvent.getByToken(triggerResultToken_, triggerResults) &&
	  iEvent.getByToken(triggerToken_  , triggerEvent)) ) {
     edm::LogError("") << "Trigger collections not found !!!";
     return;
   }
   if ( !iEvent.getByToken(triggerPrescales_, triggerPrescalesH) ) {
     edm::LogError("") << "Prescale collections not found !!!";
     return;
   }
   pat::PackedTriggerPrescales triggerPrescales = *triggerPrescalesH;


   const edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResults);
   for (unsigned int itrig=0; itrig < triggerNames_.size(); ++itrig) {
     std::string iName = triggerNames_.triggerName(itrig);
     iName.erase( std::remove_if(iName.end()-2, iName.end(), (int(*)(int))std::isdigit), iName.end() );
     if ( iName.find ("Jpsi")     != std::string::npos ||
	  iName.find ("Upsilon")  != std::string::npos ||
	  iName.find ("PsiPrime") != std::string::npos ||
	  iName.find ("LowMass")  != std::string::npos ||
	  iName.find ("Phi")      != std::string::npos ||
	  iName.find ("Bs")       != std::string::npos ||
	  iName.find ("Tau3Mu")   != std::string::npos ||
	  iName.find ("Tau3mu")   != std::string::npos ||
	  iName.find ("Onia")     != std::string::npos ) {


       
       // auto PSdetails = hltPrescale_->prescaleValuesInDetail(iEvent,iSetup,triggerNames_.triggerName(itrig));
       // std::cout << triggerNames_.triggerName(itrig) << " " 
       // 		 << hltPrescale_->prescaleValue (iEvent,iSetup,triggerNames_.triggerName(itrig)) << " "
       // 		 << hltPrescale_->prescaleValues(iEvent,iSetup,triggerNames_.triggerName(itrig)).first << " " 
       // 		 << hltPrescale_->prescaleValues(iEvent,iSetup,triggerNames_.triggerName(itrig)).second << " "
       // 		 << hltPrescale_->prescaleValuesInDetail(iEvent,iSetup,triggerNames_.triggerName(itrig)).second << " "
       // 		 << PSdetails.second << std::endl;

       // std::cout << "L1 prescales:";
       // float EffectiveL1PS = 1.e9;
       // for (size_t iSeed=0; iSeed < PSdetails.first.size(); ++iSeed) {
       // 	 std::string l1_den = PSdetails.first.at(iSeed).first;
       // 	 int l1_denp = PSdetails.first.at(iSeed).second;
       // 	 std::cout << " " << l1_denp;

       // 	 if (l1_denp==1) {EffectiveL1PS = 1; break;}
       // 	 if (l1_denp==0) continue;

       // 	 EffectiveL1PS = 1. / ( 1./l1_denp + 1./EffectiveL1PS );
       // }
       // std::cout << " -> " << EffectiveL1PS << std::endl;

       triggerPrescales.setTriggerNames(triggerNames_);
       int hltPres = triggerPrescales.getPrescaleForName(iName,true);
       if ( hltPres < 1 ) hltPres = 1;

       hltNames    .push_back(iName);
       hltResults  .push_back(triggerResults->accept(itrig));
       hltMatchDimu.push_back(false);
       hltMatchB   .push_back(false);
       hltPrescales.push_back(hltPres);
       // hltPrescales.push_back(PSdetails.second);
       // l1tPrescales.push_back(EffectiveL1PS);
     }

   }

   //// L1 information
   edm::Handle<GlobalAlgBlkBxCollection> l1results;
   iEvent.getByToken(l1results_, l1results);
   
   if (l1results.isValid()) {  
     fGtUtil->retrieveL1(iEvent, iSetup, l1results_);
     const std::vector<std::pair<std::string, bool> > finalDecisions = fGtUtil->decisionsFinal();
     const std::vector<std::pair<std::string, int> >  prescales = fGtUtil->prescales();
     
     for (unsigned int i = 0; i < finalDecisions.size(); ++i) {
       std::string name = (finalDecisions.at(i)).first;
       if (name == "NULL" ||
	   (name.find("DoubleMu")==std::string::npos &&
	    name.find("SingleMu")==std::string::npos &&
	    name.find("TripleMu")==std::string::npos)
	   ) continue;
       bool resultFin = (finalDecisions.at(i)).second;
       if (resultFin){
	 l1tNames    .push_back(name);
	 l1tPrescales.push_back((prescales.at(i)).second);
       }
     }
   }
   

   edm::Handle< std::vector<pat::Muon> >  muons;
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByToken(offlineMuonsToken_, muons);
   iEvent.getByToken(offlineTksToken_,   tracks);

   //get offline beamspot position
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByToken(beamspotToken_,recoBeamSpotHandle);
   const reco::BeamSpot& vertexBeamSpot = *recoBeamSpotHandle;

   //get the transient track builder
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

   edm::Handle<reco::VertexCollection> vtxColl; 
   iEvent.getByToken(offlinePVToken_, vtxColl);
   const reco::VertexCollection pvColl = *(vtxColl.product()) ;
   nOffVtx = pvColl.size();

   //get the b field
   std::string mfName_ = "";
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get("", bFieldHandle);  
   const MagneticField* magField = bFieldHandle.product();
   TSCBLBuilderNoMaterial blsBuilder;

   int nDimuon = 0;
   const pat::Muon* muObj1 = 0;
   const pat::Muon* muObj2 = 0;

   for(std::vector<pat::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) { 
     if( muon::isSoftMuon( (*mu1), pvColl[0] ) && (*mu1).pt() > 0 && fabs( (*mu1).eta() ) < 2.4 ) {
       for(std::vector<pat::Muon>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2) {
	 if( mu2->pt() >= mu1->pt() ) continue; 
	 if( muon::isSoftMuon( (*mu2), pvColl[0]) && (*mu2).pt() > 0 && fabs( (*mu2).eta() ) < 2.4 ) {
	   if(!( mu1->charge() * mu2->charge() < 0 )) continue; 

	   // do dimuon vertex fit
	   std::vector<reco::TransientTrack> j_tks;
	   j_tks.push_back((*theB).build(mu1->track().get()));
	   j_tks.push_back((*theB).build(mu2->track().get()));
	   if (j_tks.size()!=2) continue;

	   KalmanVertexFitter jkvf;
	   TransientVertex jtv = jkvf.vertex(j_tks);
	   if (!jtv.isValid()) continue;

	   reco::Vertex dimuonvertex = jtv;
	   float DiMu_CL = 0;
	   if( (dimuonvertex.chi2()>=0.0) && (dimuonvertex.ndof()>0) ) DiMu_CL = TMath::Prob(dimuonvertex.chi2(), dimuonvertex.ndof() );

	   if ( DiMu_CL < 0.001 ) continue;

	   math::XYZVector dimuonperp(mu1->px() + mu2->px() ,
				  mu1->py() + mu2->py() ,
				  0.);
         
	   reco::Particle::LorentzVector p1 = reco::Particle::LorentzVector( mu1->px(), mu1->py(), mu1->pz(), sqrt( mu1->momentum().Mag2() + 0.106*0.106 ) );
	   reco::Particle::LorentzVector p2 = reco::Particle::LorentzVector( mu2->px(), mu2->py(), mu2->pz(), sqrt( mu2->momentum().Mag2() + 0.106*0.106 ) );
	   reco::Particle::LorentzVector pDimuon = p1 + p2;

	   GlobalPoint jVertex = jtv.position();
	   GlobalError jerr    = jtv.positionError();
          
	   //calculate decay length  significance w.r.t. the beamspot
	   GlobalPoint displacementFromBeamspotDimuon( -1*((vertexBeamSpot.x0() - jVertex.x()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
						       -1*((vertexBeamSpot.y0() - jVertex.y()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
						       0);
	   reco::Vertex::Point vperpj(displacementFromBeamspotDimuon.x(), displacementFromBeamspotDimuon.y(), 0.);

	   float DiMu_ls   = displacementFromBeamspotDimuon.perp() /sqrt(jerr.rerr(displacementFromBeamspotDimuon));
	   float DiMu_cos  = vperpj.Dot(dimuonperp)/(vperpj.R()*dimuonperp.R());
	   float DiMu_mass = (mu1->p4() + mu2->p4()).mass();

	   // std::cout << "DiMu CL: "         << DiMu_CL   << std::endl;
	   // std::cout << "DiMu L_sig: "      << DiMu_ls   << std::endl;
	   // std::cout << "DiMu cos(alpha): " << DiMu_cos  << std::endl;
	   // std::cout << "DiMu mass: "       << DiMu_mass << std::endl;
	   // std::cout << "DiMu pT: "         << dimuonperp.R() << std::endl;


	   nDimuon++;

	   if (DiMu_CL > Dimuon_CL) {

	     Dimuon_CL       = DiMu_CL;
	     Dimuon_LS       = DiMu_ls;
	     Dimuon_CosAlpha = DiMu_cos;
	     Dimuon_Mass     = DiMu_mass;
	     Dimuon_pT       = pDimuon.pt();
	     Dimuon_eta      = pDimuon.eta();
	     Dimuon_phi      = pDimuon.phi();

	     Muon1_pT  = mu1->pt();
	     Muon1_eta = mu1->eta();
	     Muon1_phi = mu1->phi();
	     Muon1_cha = mu1->charge();
	     Muon2_pT  = mu2->pt();
	     Muon2_eta = mu2->eta();
	     Muon2_phi = mu2->phi();
	     Muon2_cha = mu2->charge();

	     muObj1 = &(*mu1);
	     muObj2 = &(*mu2);

	     for (unsigned iTrig=0; iTrig<hltNames.size(); ++iTrig)
	       if ( hltResults.at(iTrig) ) {
		 auto mapIter = muonFilterMap->find(hltNames.at(iTrig));
		 if ( mapIter == muonFilterMap->end() ) continue;
		 bool m1_matched = false;
		 bool m2_matched = false;
		 for (unsigned i=0; i<triggerEvent->size(); ++i) {
		   pat::TriggerObjectStandAlone obj = triggerEvent->at(i);
		   obj.unpackFilterLabels(iEvent,*triggerResults);
		   if (obj.hasFilterLabel(mapIter->second)) {
		     if ( deltaR(mu1->eta(),mu1->phi(),obj.eta(),obj.phi())<0.01 ) m1_matched = true;
		     if ( deltaR(mu2->eta(),mu2->phi(),obj.eta(),obj.phi())<0.01 ) m2_matched = true;
		     if ( m1_matched && m2_matched ) {
		       hltMatchDimu[iTrig] = true;
		       break;
		     }
		   }
		 }
	       }

	   }
	 }
       }
     }
   }

   // std::cout << "# of Dimuon: " << nDimuon << std::endl;
   if ( nDimuon == 0 ) return;

   int nBp = 0;

   if ( doBp && fabs(Dimuon_Mass-3.0969) < 0.2 ) {
     // cout<<"Start dimuon vtx fit. Mu1: "<<muObj1->pt()<<" "<<muObj1->eta()<<" "<<muObj1->phi()
     // 	 <<" Mu2:  "<<muObj2->pt()<<" "<<muObj2->eta()<<" "<<muObj2->phi()<<" "<<endl;
     // preselection on tracks
     selTracksDef qualityTracks;
     for (uint tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++) {
       reco::TrackRef checkTrk(tracks,tracksIt) ;                                                
       reco::Track itrk  = tracks->at(tracksIt) ;                                                
    
       if (!checkTrk->quality(reco::TrackBase::highPurity))                continue;
       if (checkTrk->pt() < 0.8 || fabs(checkTrk->eta())  > 2.4) continue;            
    
       // FreeTrajectoryState InitialFTS2 = initialFreeState(itrk, magField);
       // TrajectoryStateClosestToBeamLine tscb2( blsBuilder(InitialFTS2, *recoBeamSpotHandle) );
       // float trk_d0sig = tscb2.transverseImpactParameter().significance();
       // if (trk_d0sig  < 2)    continue;

       // check overlap with muon collection
       bool flag = false ;                                                         
       for (std::vector<pat::Muon>::const_iterator mu =muons->begin(); mu!=muons->end(); mu++)                    
	 {                                                                
	   if (mu->track().get()!= 0 && mu->track().get() == checkTrk.get()) { flag=true; break; }                    
	 }                                                                
       if (flag)   continue;            

       qualityTracks[tracksIt] = itrk;     
     }

     // Loop on track collection - trk 1
     for (selTracksDef::const_iterator tracksIt=qualityTracks.begin(); tracksIt!=qualityTracks.end(); ++tracksIt) {
       reco::Track itrk1((*tracksIt).second) ;
       if (deltaR(muObj1->eta(),muObj1->phi(),itrk1.eta(),itrk1.phi())<.001)    continue;
       if (deltaR(muObj2->eta(),muObj2->phi(),itrk1.eta(),itrk1.phi())<.001)    continue;
            
       FreeTrajectoryState InitialFTS = initialFreeState(itrk1, magField);
       TrajectoryStateClosestToBeamLine tscb( blsBuilder(InitialFTS, *recoBeamSpotHandle) );

       // hists_["trkPt"] -> Fill( itrk1.pt() );
       // hists_["D0sig"] -> Fill( trk1_d0sig );
 
       reco::Particle::LorentzVector pB, p1, p2, p3;

       // Combined system
       double thirdTrackMass2 = 0.493677 * 0.493677;
       double MuMass2         = 0.105658 * 0.105658;
       double e1 = sqrt( muObj1->momentum().Mag2() + MuMass2         );
       double e2 = sqrt( muObj2->momentum().Mag2() + MuMass2         );
       double e3 = sqrt( itrk1.momentum().Mag2()   + thirdTrackMass2 );
            
       p1 = reco::Particle::LorentzVector(muObj1->px(), muObj1->py(), muObj1->pz(), e1);
       p2 = reco::Particle::LorentzVector(muObj2->px(), muObj2->py(), muObj2->pz(), e2);
       p3 = reco::Particle::LorentzVector(itrk1.px()  , itrk1.py()  , itrk1.pz()  , e3);
            
       pB = p1 + p2 + p3;
            
       if (pB.mass() > 5.5 || pB.mass() < 5) continue;

       // do the vertex fit
       std::vector<reco::TransientTrack> t_tks;
       t_tks.push_back((*theB).build(muObj1->track().get()));
       t_tks.push_back((*theB).build(muObj2->track().get()));
       t_tks.push_back((*theB).build(&itrk1));
       if (t_tks.size()!=3) continue;
            
       KalmanVertexFitter kvf; //cout<<"Vertex fit. Trk: "<<itrk1.pt()<<" "<<itrk1.eta()<<" "<<itrk1.phi()<<" "<<muObj1->track()->pt()<<" "<<muObj1->track()->eta()<<" "<<muObj1->track()->phi()<<" "<<muObj2->track()->pt()<<" "<<muObj2->track()->eta()<<" "<<muObj2->track()->phi()<<endl;
       for (std::vector<reco::TransientTrack>::iterator iter = t_tks.begin(); iter!=t_tks.end(); ++iter)
	 // cout<<iter->track().pt()<<" "<<iter->track().eta()<<" "<<iter->track().phi()<<" "<<iter->track().lost()<<" "<<iter->track().found()<<" "<<iter->track().numberOfValidHits()<<endl;
       TransientVertex tv  = kvf.vertex(t_tks);
       // cout<<"Done!"<<endl;
       if (!tv.isValid()) continue;
       reco::Vertex vertex = tv;
       // hists_["B0InvMass"]->Fill( pB.mass() );
       float JpsiTkCL = 0;
       if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   
	 JpsiTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );

       if ( JpsiTkCL < 0.001 ) continue;
              
       // calculate four-track transverse momentum
       math::XYZVector pperp(muObj1->px() + muObj2->px() + itrk1.px(),
			     muObj1->py() + muObj2->py() + itrk1.py(),
			     0.);
       // get vertex position and error to calculate the decay length significance
       GlobalPoint secondaryVertex = tv.position();
       GlobalError err             = tv.positionError();
       GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + 
						 (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
					     -1*((vertexBeamSpot.y0() - secondaryVertex.y()) + 
						 (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 
					     0);
       reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);

       ++nBp;

       if (JpsiTkCL > Bp_CL) {
	 Bp_CL       = JpsiTkCL;
	 Bp_LS       = displacementFromBeamspot.perp() / sqrt(err.rerr(displacementFromBeamspot));
	 Bp_CosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	 Bp_Mass     = pB.mass();
	 Bp_pT       = pB.pt();
	 Bp_eta      = pB.eta();
	 Bp_phi      = pB.phi();
       
	 track_pT  = itrk1.pt();
	 track_eta = itrk1.eta();
	 track_phi = itrk1.phi();
	 track_cha = itrk1.charge();
	 track_d0  = tscb.transverseImpactParameter().significance();
 
	 for (unsigned iTrig=0; iTrig<hltNames.size(); ++iTrig)
	   if ( hltResults.at(iTrig) && hltMatchDimu.at(iTrig) ) {
	     auto mapIter = trackFilterMap->find(hltNames.at(iTrig));
	     if ( mapIter == trackFilterMap->end() ) continue;
	     for (unsigned i=0; i<triggerEvent->size(); ++i) {
	       pat::TriggerObjectStandAlone obj = triggerEvent->at(i);
	       obj.unpackFilterLabels(iEvent,*triggerResults);
	       if ( obj.hasFilterLabel(mapIter->second) && deltaR(itrk1.eta(),itrk1.phi(),obj.eta(),obj.phi())<0.01 ) {
		 hltMatchB[iTrig] = true;
		 break;
	       }
	     }
	   }

       }
     }
   }

   // std::cout << "# of B+: " << nBp << ", CLmax=" <<Bp_CL << ", trk_d0=" << track_d0 << std::endl;
   // if ( nBp > 0 ) {
   //   std::cout << Bp_CL << std::endl;
   //   std::cout << Bp_LS << std::endl;
   //   std::cout << Bp_CosAlpha << std::endl;
   //   std::cout << Bp_Mass << std::endl;
   //   std::cout << Bp_pT << std::endl;
   //   std::cout << Bp_eta << std::endl;
   //   std::cout << Bp_phi << std::endl;

   //   std::cout << track_pT << std::endl;
   //   std::cout << track_eta << std::endl;
   //   std::cout << track_phi << std::endl;
   //   std::cout << track_cha << std::endl;
   // }
   // if ( nBp == 0 ) return;


   // edm::Handle<l1t::MuonBxCollection> l1cands;

   // if (iEvent.getByToken(l1candToken_, l1cands))
   //   for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
   //     if (ibx != 0) continue;
   //     for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++){

   // 	 l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );
   // 	 // std::cout << muon->pt() << std::endl;

   //     }
   //   }
   // else
   //   edm::LogWarning("") << "Online L1 muon collection not found !!!";

   outTree->Fill();

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}

float
VanillaHLTAnalyzer::deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dPhi = fabs(phi1 - phi2);
  if ( dPhi > TMath::Pi() )  dPhi = 2*TMath::Pi() - dPhi;

  return sqrt( (eta1-eta2)*(eta1-eta2) + dPhi*dPhi );
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
  bool changed(true);
  if (hltPrescale_.init(run,iSetup,"HLT",changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run
      std::cout << "Initializing HLTConfigProvider"  << std::endl;
    }
  } 
  else {
    std::cout << " HLT config extraction failure with process name HLT" << std::endl;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
VanillaHLTAnalyzer::beginJob()
{
  outTree = outfile_-> make<TTree>("ntupleTree","ntupleTree");

  outTree->Branch("hltNames"    , &hltNames     );
  outTree->Branch("hltResults"  , &hltResults   );
  outTree->Branch("hltMatchDimu", &hltMatchDimu );
  outTree->Branch("hltMatchB"   , &hltMatchB    );
  outTree->Branch("hltPrescales", &hltPrescales );
  outTree->Branch("l1tNames"    , &l1tNames     );
  outTree->Branch("l1tPrescales", &l1tPrescales );

  outTree->Branch("Dimuon_CL"      , &Dimuon_CL       );
  outTree->Branch("Dimuon_LS"      , &Dimuon_LS       );
  outTree->Branch("Dimuon_CosAlpha", &Dimuon_CosAlpha );
  outTree->Branch("Dimuon_Mass"    , &Dimuon_Mass     );
  outTree->Branch("Dimuon_pT"      , &Dimuon_pT       );
  outTree->Branch("Dimuon_eta"     , &Dimuon_eta      );
  outTree->Branch("Dimuon_phi"     , &Dimuon_phi      );

  outTree->Branch("Muon1_pT" , &Muon1_pT  );
  outTree->Branch("Muon1_eta", &Muon1_eta );
  outTree->Branch("Muon1_phi", &Muon1_phi );
  outTree->Branch("Muon1_cha", &Muon1_cha );
  outTree->Branch("Muon2_pT" , &Muon2_pT  );
  outTree->Branch("Muon2_eta", &Muon2_eta );
  outTree->Branch("Muon2_phi", &Muon2_phi );
  outTree->Branch("Muon2_cha", &Muon2_cha );

  outTree->Branch("Bp_CL"      , &Bp_CL       );
  outTree->Branch("Bp_LS"      , &Bp_LS       );
  outTree->Branch("Bp_CosAlpha", &Bp_CosAlpha );
  outTree->Branch("Bp_Mass"    , &Bp_Mass     );
  outTree->Branch("Bp_pT"      , &Bp_pT       );
  outTree->Branch("Bp_eta"     , &Bp_eta      );
  outTree->Branch("Bp_phi"     , &Bp_phi      );

  outTree->Branch("track_pT" , &track_pT  );
  outTree->Branch("track_eta", &track_eta );
  outTree->Branch("track_phi", &track_phi );
  outTree->Branch("track_cha", &track_cha );
  outTree->Branch("track_d0" , &track_d0  );

  outTree->Branch("nOffVtx", &nOffVtx );

  muonFilterMap = new std::map<std::string, std::string>();

  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_3_Bs_v","hltDisplacedmumuFilterDoubleMu4Bs"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_3_Jpsi_Displaced_v","hltDisplacedmumuFilterDoubleMu43Jpsi"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4LowMassNonResonant"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_v","hltDisplacedmumuFilterDoubleMu3Tau3mu"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4PsiPrime"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Jpsi_v","hltL3fSQMu7p5L2Mu2L3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Upsilon_v","hltL3fSQMu7p5L2Mu2L3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Jpsi_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Upsilon_v","hltL3fLMu7p5TrackL3Filtered7p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Onia_v","hltL3fL1sMu22orMu20erorMu25L1f0L2f0L3Filtered25"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu30_TkMu0_Onia_v","hltL3fL1sMu22orMu20erorMu25L1f0L2f0L3Filtered30"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu20_TkMu0_Phi_v","hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Phi_v","hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered25"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_L1_NoOS_v","hltDisplacedmumuFilterDimuon0JpsiL1sNoOS"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v","hltDimuon0JpsiNoVtxNoOSL3Filtered"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_v","hltDisplacedmumuFilterDimuon0Jpsi"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_v","hltDimuon0JpsiL3Filtered"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v","hltDisplacedmumuFilterDimuon0JpsiL1s4R0er1p5R"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v","hltDimuon0JpsiL1s4R0er1p5RL3Filtered"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Jpsi3p5_Muon2_v","hltVertexmumuFilterJpsiMuon3p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_5_v","hltDisplacedmumuFilterDimuon0UpsilonL1s5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5NoOS_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5NoOS"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5er2p0_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5er2p0"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v","hltDisplacedmumuFilterDimuon0UpsilonL1s4p5er2p0M"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_NoVertexing_v","hltDimuon0UpsilonL1s4p5er2p0ML3Filtered"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_L1_5M_v","hltDisplacedmumuFilterDimuon0UpsilonL1s5M"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_0er1p5R_v","hltDisplacedmumuFilterDimuon0LowMassL1s0er1p5R"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_0er1p5_v","hltDisplacedmumuFilterDimuon0LowMassL1s0er1p5"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_v","hltDisplacedmumuFilterDimuon0LowMass"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_4_v","hltDisplacedmumuFilterDimuon0LowMassL1s4"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_4R_v","hltDisplacedmumuFilterDimuon0LowMassL1s4R"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_LowMass_L1_TM530_v","hltDisplacedmumuFilterDimuon0LowMassL1sTM530"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_Muon_L1_TM0_v","hltVertexmumuFilterUpsilon0MuonL1sTM0"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v","hltVertexmumuFilterUpsilon0MuonNoL1Mass"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v","hltDisplacedmumuFilterDoubleMu3Tau3muNoL1Mass"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_Jpsi_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_Jpsi_NoVertexing_v","hltDoubleMu4JpsiDisplacedL3Filtered"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v","hltDisplacedmumuFilterDoubleMu4Jpsi"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon10PsiPrimeBarrelnoCow"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon10_Upsilon_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon10UpsilonBarrelnoCow"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon12_Upsilon_eta1p5_v","hltDisplacedmumuFilterDimuon12Upsilons"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon14_Phi_Barrel_Seagulls_v","hltDisplacedmumuFilterDimuon14PhiBarrelnoCow"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon18_PsiPrime_v","hltDisplacedmumuFilterDimuon18PsiPrimes"));
  muonFilterMap->insert(std::pair<std::string, std::string>("HLT_Dimuon25_Jpsi_v","hltDisplacedmumuFilterDimuon25Jpsis"));

  extraFilterMap = new std::map<std::string, std::string>();

  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Jpsi_v","hltSQMu7p5L2Mu2JpsiTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_L2Mu2_Upsilon_v","hltSQMu7p5L2Mu2UpsilonTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Jpsi_v","hltMu7p5Track2JpsiTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Jpsi_v","hltMu7p5Track3p5JpsiTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Jpsi_v","hltMu7p5Track7JpsiTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track2_Upsilon_v","hltMu7p5Track2UpsilonTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track3p5_Upsilon_v","hltMu7p5Track3p5UpsilonTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu7p5_Track7_Upsilon_v","hltMu7p5Track7UpsilonTrackMassFiltered"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Onia_v","hltDiMuonGlbFiltered25TrkFiltered0"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu30_TkMu0_Onia_v","hltDiMuonGlbFiltered30TrkFiltered0"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu20_TkMu0_Phi_v","hltDiMuonGlbFiltered20TrkFiltered0"));
  extraFilterMap->insert(std::pair<std::string, std::string>("HLT_Mu25_TkMu0_Phi_v","hltDiMuonGlbFiltered25PhiTrkFiltered0"));

  trackFilterMap = new std::map<std::string, std::string>();

  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrk_Displaced_v","hltJpsiTkVertexFilter"));
  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v","hltLowMassNonResonantTkVertexFilter"));
  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_v","hltTau3muTkVertexFilter"));
  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v","hltPsiPrimeTkVertexFilter"));
  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v","hltTau3muNoL1MassTkVertexFilter"));
  trackFilterMap->insert(std::pair<std::string, std::string>("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v","hltJpsiTkTkVertexFilterPhiKstar"));
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
