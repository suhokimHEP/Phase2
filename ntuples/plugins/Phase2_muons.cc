#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Phase2/ntuples/interface/Phase2.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
using namespace std;
using namespace reco;
using namespace pat;

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
vector<float> EtaPhi;





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
 tree->Branch("muEtaPhi",                     &muEtaPhi_                     ) ; 
}

// initialize branches
void Phase2::branchesTracks(TTree* tree) {
  tree->Branch("nTrack"                   , &nTrack_);
  tree->Branch("AllTrackPt"               , &AllTrackPt_);
  tree->Branch("AllTrackEta"               , &AllTrackEta_);
  tree->Branch("AllTrackPhi"               , &AllTrackPhi_);
  tree->Branch("AllTrackdEdx"               , &AllTrackdEdx_);

}





void Phase2::fillMuons(const edm::Event& e, reco::Vertex vtx, string Mode) {

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
 muEtaPhi_.clear();
 edm::Handle<edm::View<reco::Muon> > muonHandle;
 edm::Handle<edm::View<pat::Muon> > slimmuonHandle;

  if(Mode.find("reco")!= std::string::npos){
 e.getByToken(muonCollection_, muonHandle);
 for (edm::View<reco::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) 


{
  EtaPhi.clear();
  Float_t pt = iMu->pt();
  Float_t eta = iMu->eta();
  Float_t phi = iMu->phi();

  if (pt < 5) continue;
  //if (fabs(eta) > 3.0) continue;
  //if (fabs(eta) < 1.4) continue;
  if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;

  const reco::Muon &recoMu = dynamic_cast<const reco::Muon &>(*iMu);
  nMu_++;

  if( !recoMu.innerTrack().isNull() ){
   muNumberOfMissingInnerHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_INNER_HITS) ) ; 
   muNumberOfMissingMiddleHits_   .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::TRACK_HITS        ) ) ; 
   muNumberOfMissingOuterHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_OUTER_HITS) ) ; 
   muNumberOfValidPixelHits_      .push_back( recoMu.innerTrack ()->hitPattern ().numberOfValidPixelHits       ()   ) ; 
   muTrackerLayersWithMeasurement_.push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithMeasurement ()   ) ; 
  }


  muPt_    .push_back(pt);
  muEn_    .push_back(iMu->energy());
  muEta_   .push_back(iMu->eta());
  muPhi_   .push_back(iMu->phi());
  muCharge_.push_back(iMu->charge());
  muType_  .push_back(iMu->type());
  EtaPhi.push_back(iMu->eta());
  EtaPhi.push_back(iMu->phi());
  muEtaPhi_.push_back(EtaPhi);
	//	std::cout<<nMu_<<"Eta:"<<EtaPhi.at(0)<<",Phi:"<<EtaPhi.at(1)<<std::endl;
  muIsGlobalMuon_                .push_back(iMu->  isGlobalMuon () ) ; 
  muIsPFMuon_                    .push_back(iMu->  isPFMuon     () ) ; 

  Float_t muPFChIso      = iMu->pfIsolationR04().sumChargedHadronPt ;
  Float_t muPFPhoIso     = iMu->pfIsolationR04().sumPhotonEt        ;
  Float_t muPFNeuIso     = iMu->pfIsolationR04().sumNeutralHadronEt ;
  Float_t muPFPUIso      = iMu->pfIsolationR04().sumPUPt            ;
  Float_t pfdBetaIso     = ( muPFChIso + max(0.0,muPFNeuIso + muPFPhoIso - 0.5*muPFPUIso ) ) / pt ;

  muPFdBetaIsolation_     .push_back( pfdBetaIso     ) ; 

 muTrackdR_.clear();
 if(AllTrackEta_.size()>0) muTrackdR_ = CalTrackdR(eta,phi);
 //MinmuTrackdR_ = *min_element(muTrackdR_.begin(), muTrackdR_.end());
} 
}



else {
 e.getByToken(slimmuonCollection_, slimmuonHandle);
 for (edm::View<pat::Muon>::const_iterator iMu = slimmuonHandle->begin(); iMu != slimmuonHandle->end(); ++iMu) 


{
  EtaPhi.clear();
  Float_t pt = iMu->pt();
  Float_t eta = iMu->eta();
  Float_t phi = iMu->phi();

  if (pt < 5) continue;
  //if (fabs(eta) > 3.0) continue;
  //if (fabs(eta) < 1.4) continue;
  if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;

  const reco::Muon &recoMu = dynamic_cast<const reco::Muon &>(*iMu);
  nMu_++;

  if( !recoMu.innerTrack().isNull() ){
   muNumberOfMissingInnerHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_INNER_HITS) ) ; 
   muNumberOfMissingMiddleHits_   .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::TRACK_HITS        ) ) ; 
   muNumberOfMissingOuterHits_    .push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithoutMeasurement (reco::HitPattern::MISSING_OUTER_HITS) ) ; 
   muNumberOfValidPixelHits_      .push_back( recoMu.innerTrack ()->hitPattern ().numberOfValidPixelHits       ()   ) ; 
   muTrackerLayersWithMeasurement_.push_back( recoMu.innerTrack ()->hitPattern ().trackerLayersWithMeasurement ()   ) ; 
  }

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
  EtaPhi.push_back(iMu->eta());
  EtaPhi.push_back(iMu->phi());
  muEtaPhi_.push_back(EtaPhi);
  UShort_t tmpmuIDbit = 0;
	//	std::cout<<nMu_<<"Eta:"<<EtaPhi.at(0)<<",Phi:"<<EtaPhi.at(1)<<std::endl;
  if(Mode.find("MiniAOD")!= std::string::npos){

  if (iMu->isLooseMuon())     setbit(tmpmuIDbit, 1);
  if (iMu->isMediumMuon())    setbit(tmpmuIDbit, 2);
  if (iMu->isTightMuon(vtx))  setbit(tmpmuIDbit, 2);
  if (iMu->isSoftMuon(vtx))   setbit(tmpmuIDbit, 3);
  if (iMu->isHighPtMuon(vtx)) setbit(tmpmuIDbit, 4);
  muIDbit_.push_back(tmpmuIDbit);

  muPassHipID_ .push_back( isMedium )   ;
}
  muIsGlobalMuon_                .push_back(iMu->  isGlobalMuon () ) ; 
  muIsPFMuon_                    .push_back(iMu->  isPFMuon     () ) ; 

  Float_t muPFChIso      = iMu->pfIsolationR04().sumChargedHadronPt ;
  Float_t muPFPhoIso     = iMu->pfIsolationR04().sumPhotonEt        ;
  Float_t muPFNeuIso     = iMu->pfIsolationR04().sumNeutralHadronEt ;
  Float_t muPFPUIso      = iMu->pfIsolationR04().sumPUPt            ;
  Float_t pfdBetaIso     = ( muPFChIso + max(0.0,muPFNeuIso + muPFPhoIso - 0.5*muPFPUIso ) ) / pt ;

  muPFdBetaIsolation_     .push_back( pfdBetaIso     ) ; 

 muTrackdR_.clear();
 if(AllTrackEta_.size()>0) muTrackdR_ = CalTrackdR(eta,phi);
 //MinmuTrackdR_ = *min_element(muTrackdR_.begin(), muTrackdR_.end());
} 
}
 if (!muonHandle.isValid()&&!slimmuonHandle.isValid() ) {
  edm::LogWarning("Phase2") << "no Muons in event";
  return;
 }

}


void Phase2::fillTracks(const edm::Event& e, const edm::EventSetup& es, string Mode) {

 nTrack_=0;
 AllTrackPt_.clear();
 AllTrackEta_.clear();
 AllTrackPhi_.clear();
edm::Handle<std::vector<reco::Track>   >  TrackHandle;
edm::Handle<std::vector<pat::IsolatedTrack>   >  isoTrackHandle;


  if(Mode.find("reco")!= std::string::npos){
 e.getByToken( TrackLabel_       ,  TrackHandle );
 for (std::vector<reco::Track>::const_iterator itrack = TrackHandle->begin(); itrack != TrackHandle->end(); ++itrack) {
  nTrack_++;
  Float_t pt = itrack->pt();
  Float_t eta = itrack->eta();
  Float_t phi = itrack->phi();
  //Float_t pfNeutralSum = itrack->pfNeutralSum();
  AllTrackPt_.push_back(pt);
  AllTrackEta_.push_back(eta);
  AllTrackPhi_.push_back(phi);

}
}


else{ e.getByToken( isoTrackLabel_       ,  isoTrackHandle );
 for (std::vector<pat::IsolatedTrack>::const_iterator itrack = isoTrackHandle->begin(); itrack != isoTrackHandle->end(); ++itrack) {
  nTrack_++;
  Float_t pt = itrack->pt();
  Float_t eta = itrack->eta();
  Float_t phi = itrack->phi();
  //Float_t pfNeutralSum = itrack->pfNeutralSum();
  Float_t dEdxPixel = itrack->dEdxPixel();
  AllTrackPt_.push_back(pt);
  AllTrackEta_.push_back(eta);
  AllTrackPhi_.push_back(phi);
  AllTrackdEdx_.push_back(dEdxPixel);

}
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

