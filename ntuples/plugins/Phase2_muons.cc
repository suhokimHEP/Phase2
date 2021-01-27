////#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Phase2/ntuples/interface/Phase2.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
using namespace std;

// (local) variables associated with tree branches
Int_t            nMu_                             ; 
vector<float>    muPt_                            ; 
vector<float>    muEn_                            ; 
vector<float>    muEta_                           ; 
vector<float>    muPhi_                           ; 
vector<int>      muCharge_                        ; 
vector<int>      muType_                          ; 
vector<UShort_t> muIDbit_                         ; 
vector<bool>     muPassLooseID_                   ; 
vector<bool>     muPassHipID_                     ; 
vector<bool>     muPassTighID_                    ; 
// vector<bool>     muPassSoftID_                   ; 
// vector<bool>     muPassHighPtID_                 ; 
// vector<float>    muPVIndex_                      ; 
vector<int>     muNumberOfMissingInnerHits_       ; 
vector<int>     muNumberOfMissingMiddleHits_      ; 
vector<int>     muNumberOfMissingOuterHits_       ; 
//vector<int>     muNumberOfValidHits_              ; 
//vector<float>   muNormalizedChi2_                 ; 
vector<int>     muNumberOfMatchedStations_        ; 
vector<int>     muNumberOfValidPixelHits_         ; 
vector<int>     muTrackerLayersWithMeasurement_   ; 
vector<int>     muIsGlobalMuon_                   ; 
vector<int>     muIsPFMuon_                       ; 
//vector<int>     muIsTightMuonWRTVtx_              ; 
vector<float>   muPFdBetaIsolation_               ; 
vector<float>   muTrackdR_               ; 
float   MinmuTrackdR_               ; 

Int_t nTrack_;

vector<TVector3> AllTrackPositions; // x,y,z

vector<float>    AllTrackPt_;
vector<float>    AllTrackEta_;
vector<float>    AllTrackPhi_;
vector<float>    AllTrackdEdx_;





void Phase2::branchesMuons(TTree* tree) {
 tree->Branch("nMu",                            &nMu_                            ) ; 
 tree->Branch("muPt",                           &muPt_                           ) ; 
 tree->Branch("muEn",                           &muEn_                           ) ; 
 tree->Branch("muEta",                          &muEta_                          ) ; 
 tree->Branch("muPhi",                          &muPhi_                          ) ; 
 tree->Branch("muCharge",                       &muCharge_                       ) ; 
 tree->Branch("muType",                         &muType_                         ) ; 
 tree->Branch("muIDbit",                        &muIDbit_                        ) ; 
 tree->Branch("muPassLooseID",                  &muPassLooseID_                  ) ; 
 tree->Branch("muPassHipID",                    &muPassHipID_                    ) ; 
 tree->Branch("muPassTighID",                   &muPassTighID_                   ) ; 
 //tree->Branch("muPVIndex",                      &muPVIndex_                      ) ; 
 tree->Branch("muNumberOfMissingInnerHits",     &muNumberOfMissingInnerHits_     ) ; 
 tree->Branch("muNumberOfMissingMiddleHits",    &muNumberOfMissingMiddleHits_    ) ; 
 tree->Branch("muNumberOfMissingOuterHits",     &muNumberOfMissingOuterHits_     ) ; 
 //tree->Branch("muNumberOfValidHits",            &muNumberOfValidHits_            ) ; 
 //tree->Branch("muNormalizedChi2",               &muNormalizedChi2_               ) ; 
 tree->Branch("muNumberOfMatchedStations",      &muNumberOfMatchedStations_      ) ; 
 tree->Branch("muNumberOfValidPixelHits",       &muNumberOfValidPixelHits_       ) ; 
 tree->Branch("muTrackerLayersWithMeasurement", &muTrackerLayersWithMeasurement_ ) ; 
 tree->Branch("muIsGlobalMuon",                 &muIsGlobalMuon_                 ) ; 
 tree->Branch("muIsPFMuon",                     &muIsPFMuon_                     ) ; 
 //tree->Branch("muIsTightMuonWRTVtx",            &muIsTightMuonWRTVtx_            ) ; 
 tree->Branch("muPFdBetaIsolation",             &muPFdBetaIsolation_             ) ; 
 tree->Branch("muTrackdR",                     &muTrackdR_                     ) ; 
 tree->Branch("MinmuTrackdR",                     &MinmuTrackdR_                     ) ; 
}

// initialize branches
void Phase2::branchesTracks(TTree* tree) {
  tree->Branch("nTrack"                   , &nTrack_);
  tree->Branch("AllTrackPt"               , &AllTrackPt_);
  tree->Branch("AllTrackEta"               , &AllTrackEta_);
  tree->Branch("AllTrackPhi"               , &AllTrackPhi_);
  tree->Branch("AllTrackdEdx"               , &AllTrackdEdx_);

}







void Phase2::fillMuons(const edm::Event& e, reco::Vertex vtx) {

 // cleanup from previous execution
 nMu_ = 0;
 muPt_                          .clear() ; 
 muEn_                          .clear() ; 
 muEta_                         .clear() ; 
 muPhi_                         .clear() ; 
 muCharge_                      .clear() ; 
 muType_                        .clear() ; 
 muIDbit_                       .clear() ; 
 muPassLooseID_                 .clear() ; 
 muPassHipID_                .clear() ; 
 muPassTighID_                  .clear() ; 
 //muPVIndex_                     .clear() ; 
 muNumberOfMissingInnerHits_    .clear() ; 
 muNumberOfMissingMiddleHits_   .clear() ; 
 muNumberOfMissingOuterHits_    .clear() ; 
 //muNumberOfValidHits_           .clear() ; 
 //muNormalizedChi2_              .clear() ; 
 muNumberOfMatchedStations_     .clear() ; 
 muNumberOfValidPixelHits_      .clear() ; 
 muTrackerLayersWithMeasurement_.clear() ; 
 muIsGlobalMuon_                .clear() ; 
 muIsPFMuon_                    .clear() ; 
 //muIsTightMuonWRTVtx_           .clear() ; 
 muPFdBetaIsolation_            .clear() ; 
 muTrackdR_                    .clear() ; 
 MinmuTrackdR_=0.;
 edm::Handle<edm::View<pat::Muon> > muonHandle;
 e.getByToken(muonCollection_, muonHandle);

 if (!muonHandle.isValid()) {
  edm::LogWarning("Phase2") << "no pat::Muons in event";
  return;
 }

 for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {

  nMu_++;
  Float_t pt = iMu->pt();
  Float_t eta = iMu->eta();
  Float_t phi = iMu->phi();

  //if (pt < 5) continue;
  //if (fabs(eta) > 3.0) continue;
  if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;

  const reco::Muon &recoMu = dynamic_cast<const reco::Muon &>(*iMu);

  if( !recoMu.innerTrack().isNull() ){
   muNumberOfMissingInnerHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_INNER_HITS) ) ; 
   muNumberOfMissingMiddleHits_   .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::TRACK_HITS        ) ) ; 
   muNumberOfMissingOuterHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_OUTER_HITS) ) ; 
   muNumberOfValidPixelHits_      .push_back( recoMu.innerTrack ()->hitPattern ().numberOfValidPixelHits       ()   ) ; 
   muTrackerLayersWithMeasurement_.push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithMeasurement ()   ) ; 
  }
  //muNumberOfValidHits_           .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::)              ) ; 
  //muNormalizedChi2_              .push_back( recoMu.globalTrack().normalizedChi2 ()   ) ; 

  bool goodGlob = iMu->isGlobalMuon() && 
                  iMu->globalTrack()->normalizedChi2() < 3 && 
                  iMu->combinedQuality().chi2LocalPosition < 12 && 
                  iMu->combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
                  iMu->innerTrack()->validFraction() > 0.49 && 
                  muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 

  muPt_    .push_back(pt);
  muEn_    .push_back(iMu->energy());
  muEta_   .push_back(iMu->eta());
  muPhi_   .push_back(iMu->phi());
  muCharge_.push_back(iMu->charge());
  muType_  .push_back(iMu->type());
  // muD0_    .push_back(iMu->muonBestTrack()->dxy(pv));
  // muDz_    .push_back(iMu->muonBestTrack()->dz(pv));
  // muSIP_   .push_back(fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D));

  UShort_t tmpmuIDbit = 0;

  if (iMu->isLooseMuon())     setbit(tmpmuIDbit, 1);
  if (iMu->isMediumMuon())    setbit(tmpmuIDbit, 2);
  if (iMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
  if (iMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
  if (iMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);
  muIDbit_.push_back(tmpmuIDbit);

  muPassHipID_ .push_back( isMedium )   ;

  muIsGlobalMuon_                .push_back(iMu->  isGlobalMuon () ) ; 
  muIsPFMuon_                    .push_back(iMu->  isPFMuon     () ) ; 
  //muIsTightMuonWRTVtx_           .push_back(iMu->  isTightMuonWRTVtx ()  ) ; 

  Float_t muPFChIso      = iMu->pfIsolationR04().sumChargedHadronPt ;
  Float_t muPFPhoIso     = iMu->pfIsolationR04().sumPhotonEt        ;
  Float_t muPFNeuIso     = iMu->pfIsolationR04().sumNeutralHadronEt ;
  Float_t muPFPUIso      = iMu->pfIsolationR04().sumPUPt            ;
  Float_t pfdBetaIso     = ( muPFChIso + max(0.0,muPFNeuIso + muPFPhoIso - 0.5*muPFPUIso ) ) / pt ;

  muPFdBetaIsolation_     .push_back( pfdBetaIso     ) ; 

//  const pat::PackedCandidate &ppfMu = dynamic_cast<const pat::PackedCandidate &>(*iMu);
//
//  if (!ppfMu.vertexRef().isNull() && ppfMu.vertexRef().isAvailable()){
//   muPVIndex_ .push_back(ppfMu.vertexRef().index());
//  }
 muTrackdR_.clear();
 muTrackdR_ = CalTrackdR(eta,phi);
 //MinmuTrackdR_ = *min_element(muTrackdR_.begin(), muTrackdR_.end());
 }
}


void Phase2::fillTracks(const edm::Event& e, const edm::EventSetup& es) {

 // bool dodebug = false;
 // cleanup from previous execution

 nTrack_=0;
 AllTrackPt_.clear();
 AllTrackEta_.clear();
 AllTrackPhi_.clear();

 //e.getByToken(rhoLabel_, rhoHandle);
 //float rho = *(rhoHandle.product());

// e.getByToken(vtxLabel_, vtxHandle);
// if (!vtxHandle.isValid()) edm::LogWarning("Phase2") << "Primary vertices info not unavailable";

 // ----------------------------------------------------------
 // AOD Section ----------------------------------------------

 bool verbose_AOD = false;

 // AOD Jet Handles
 //e.getByToken( AODak4CaloJetsLabel_ ,  AODak4CaloJetsHandle  );
 //e.getByToken( AODak4CaloCorrectedJetsLabel_,  AODak4CaloJetsCorrectedHandle  );
 //e.getByToken( AODak4JetCorrectorLabel_, JetCorrector );

 //e.getByToken( AODak4PFJetsLabel_   ,  AODak4PFJetsHandle    );
 //e.getByToken( AODak4PFJetsCHSLabel_,  AODak4PFJetsCHSHandle );
 //e.getByToken( selectedPatJetsLabel_,  selectedPatJetsHandle );
 //e.getByToken( AODVertexLabel_      ,  AODVertexHandle );
edm::Handle<std::vector<pat::IsolatedTrack>   >  TrackHandle;
 e.getByToken( TrackLabel_       ,  TrackHandle );
 //// Magnetic field
 //es.get<IdealMagneticFieldRecord>().get(magneticField);
 //magneticField_ = &*magneticField;
/*

 //------------------------------
 //JEC uncertainties for CaloJets
 //------------------------------
 es.get<JetCorrectionsRecord>().get("AK4Calo", JetCorrParColl_Handle);
 JetCorrectorParameters const & JetCorrPar = (*JetCorrParColl_Handle)["Uncertainty"];
 JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorrPar);

 // beamspot
 e.getByToken(beamspotLabel_, beamspotHandle_);

 // clear values
 int nMissingInner = 0;
 int nMissingOuter = 0;

 // propagator
 std::string thePropagatorName_ = "PropagatorWithMaterial";
 es.get<TrackingComponentsRecord>().get(thePropagatorName_,thePropagator_);
 StateOnTrackerBound stateOnTracker(thePropagator_.product());

 // Vertex
 vector<int> whichVertex_;
 whichVertex_.clear();
 whichVertex_ = vector<int>(AODTrackHandle->size(),-1);

 vtxfitter_ = new ConfigurableVertexReconstructor(lldj_pset_.getParameter<edm::ParameterSet>("vertexFitterConfig"));
*/
 // clear master track vectors before starting track loop
 //AODallTrackPositions       .clear();
 //AODallTrackPt              .clear();
 //AODtempTrackVector              .clear(); 
 //AODallTrackEta             .clear();
 //AODallTrackPhi             .clear();
 //AODallTrackIFSPt           .clear();
 //AODallTracknMissingInner   .clear();
 //AODallTracknMissingOuter   .clear();
 //AODallTrackAngle           .clear();
 //AODwhichVertexByTrack      .clear();
 //AODallTrackdxy             .clear();
 //AODallTrackdxyerr          .clear();

 //for(int j = 0; j < (int)TrackHandle->size(); j++){
 for (std::vector<pat::IsolatedTrack>::const_iterator itrack = TrackHandle->begin(); itrack != TrackHandle->end(); ++itrack) {
  nTrack_++;
  Float_t pt = itrack->pt();
  //std::cout<<"pt:"<<pt<<std::endl;
  Float_t eta = itrack->eta();
  //std::cout<<"eta:"<<eta<<std::endl;
  Float_t phi = itrack->phi();
  //std::cout<<"phi:"<<phi<<std::endl;
  Float_t pfNeutralSum = itrack->pfNeutralSum();
  //std::cout<<"pfNeutralSum:"<<pfNeutralSum<<std::endl;
  Float_t dEdxPixel = itrack->dEdxPixel();
  //std::cout<<"dEdxPixel:"<<dEdxPixel<<std::endl;
  AllTrackPt_.push_back(pt);
  AllTrackEta_.push_back(eta);
  AllTrackPhi_.push_back(phi);
  AllTrackdEdx_.push_back(dEdxPixel);

  //// find where track (no B field) would hit outer tracker
  //FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState(itrack,magneticField_);
  //TrajectoryStateOnSurface outer = stateOnTracker(fts);
  //if(!outer.isValid()) continue;
  //GlobalPoint outerPos = outer.globalPosition();
  //TVector3 trackPos(outerPos.x(),outerPos.y(),outerPos.z());

  //// push back track position to save in master vector
  //AODallTrackPositions.push_back(trackPos);




/*
  // get track j using the AOD track handle
  reco::TrackBaseRef tref(TrackHandle,j);
  // make transient track (unfolding effects of B field ?)
  //reco::TransientTrack tt(AODTrackHandle->at(j),magneticField_);

  //if(!tt.isValid())continue;

  // track pt first
  float trackpt  = tref->pt();

  // make selections on track
  // for alphaMax, track angle we use ttIFSpt, not tref->pt()
  //   ---------!!!!--------------
  if (trackpt < minTrackPt_) continue;  // minimum pT for track
  //if (!tref->quality(reco::TrackBase::highPurity)) continue; // track must be highPurity

  //// find where track (no B field) would hit outer tracker
  //FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState(TrackHandle->at(j),magneticField_);
  //TrajectoryStateOnSurface outer = stateOnTracker(fts);
  //if(!outer.isValid()) continue;
  //GlobalPoint outerPos = outer.globalPosition();
  //TVector3 trackPos(outerPos.x(),outerPos.y(),outerPos.z());

  //// push back track position to save in master vector
  //AODallTrackPositions.push_back(trackPos);

  // track basics (trackpt above)
  float tracketa = tref->eta();
  float trackphi = tref->phi();
  TLorentzVector tempTrackVector(0,0,0,0);
  tempTrackVector.SetPtEtaPhiM(trackpt,tracketa,trackphi,.139);
  AODtempTrackVector.push_back(tempTrackVector);
  //float ttIFSpt  = tt.initialFreeState().momentum().transverse();
  AODallTrackPt .push_back(trackpt );
  AODallTrackEta.push_back(tracketa);
  AODallTrackPhi.push_back(trackphi);}
  //AODallTrackIFSPt.push_back(ttIFSpt);

  //// Number of missing tracker hits
  //nMissingInner = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
  //nMissingOuter = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
  //AODallTracknMissingInner.push_back(nMissingInner) ;
  //AODallTracknMissingOuter.push_back(nMissingOuter) ;
*/ /*
  /// For track angle
  // get track trajectory info
  static GetTrackTrajInfo getTrackTrajInfo;
  vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(es, (*tref));
  if ( trajInfo.size() > 0 && trajInfo[0].valid) {
   // get inner tracker hit from trajectory state
   const TrajectoryStateOnSurface& tsosInnerHit = trajInfo[0].detTSOS;

   //  here's the track angle
   // find beamspot x,y coordinates
   const reco::BeamSpot& pat_beamspot = (*beamspotHandle_);
   TVector2 bmspot(pat_beamspot.x0(),pat_beamspot.y0());
   // find track trajectory state on surface inner hit
   GlobalPoint  innerPos = tsosInnerHit.globalPosition();
   float tempx = innerPos.x();
   float tempy = innerPos.y();
   InPosx_.push_back(tempx);
   InPosy_.push_back(tempy);
   GlobalVector innerMom = tsosInnerHit.globalMomentum();

   // calculate the difference between inner hit and beamspot
   TVector2 sv(innerPos.x(),innerPos.y());
   TVector2 diff = (sv-bmspot);
   //cout<<"bs x: "<<bmspot.X()<<" y: "<<bmspot.Y()<<endl;
   //cout<<" sv x: "<<sv.X()<<" y: "<<sv.Y()<<endl;
   //cout<<" diff phi: "<<diff.Phi()<<endl;
   TVector2 momentum(innerMom.x(),innerMom.y());
   //cout<<" p x: "<<momentum.X()<<" y: "<<momentum.Y()<<endl;
   //cout<<" p phi: "<<momentum.Phi()<<endl;
   //cout<<" dPhi: "<<diff.DeltaPhi(momentum)<<endl;
   float ta = fabs( diff.DeltaPhi(momentum) ) ;

   AODallTrackAngle.push_back(ta);
  }
  else{ AODallTrackAngle.push_back(0); }

  // beamspot info, track impact parameter
  float dxy = fabs(tref->dxy(*beamspotHandle_));
  float dxyerr = tref->dxyError();
  if(verbose_AOD) printf(" dxy dxyerr: %0.4f %0.4f\n", dxy, dxyerr);
  AODallTrackdxy   .push_back(dxy   ) ;
  AODallTrackdxyerr.push_back(dxyerr) ;

 }//end track loop

 //Debug printing
 if(verbose_AOD){
  for(int j = 0; j < (int)AODTrackHandle->size(); j++){
   reco::TrackBaseRef tref(AODTrackHandle,j);
   printf("AOD track pt eta phi: %f %f %f\n",tref->pt(),tref->eta(),tref->phi());
  }

  printf("  AODallTrackPositions      %i \n",  (int)AODallTrackPositions    .size() );
  printf("  AODallTrackPt             %i \n",  (int)AODallTrackPt           .size() );
  printf("  AODallTrackEta            %i \n",  (int)AODallTrackEta          .size() );
  printf("  AODallTrackPhi            %i \n",  (int)AODallTrackPhi          .size() );
  printf("  AODallTrackIFSPt          %i \n",  (int)AODallTrackIFSPt        .size() );
  printf("  AODallTracknMissingInner  %i \n",  (int)AODallTracknMissingInner.size() );
  printf("  AODallTracknMissingOuter  %i \n",  (int)AODallTracknMissingOuter.size() );
  printf("  AODallTrackAngle          %i \n",  (int)AODallTrackAngle        .size() );
  printf("  AODwhichVertexByTrack     %i \n",  (int)AODwhichVertexByTrack   .size() );
  printf("  AODallTrackdxy            %i \n",  (int)AODallTrackdxy          .size() );
  printf("  AODallTrackdxyerr         %i \n",  (int)AODallTrackdxyerr       .size() );
 }

*/
}
}

vector<float> Phase2::CalTrackdR( float jeteta, float jetphi )
{
 vector<float> idvector;
   for( int i=0; i<(int)AllTrackEta_.size(); i++){
     float tracketa = AllTrackEta_.at(i); 
       float trackphi = AllTrackPhi_.at(i); 
         float drt = deltaR( jeteta, jetphi, tracketa, trackphi );
           idvector.push_back(drt); 
            }
   return idvector;
}

