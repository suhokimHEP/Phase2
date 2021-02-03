#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Phase2/ntuples/interface/Phase2.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"


#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
edm::Handle<std::vector<PCaloHit>   >  RechHitHandle;
hgcal::RecHitTools rhtools_;


// AOD ---------------------------------------------
vector<float>  RechEnergy_;
vector<float>  RechMIP_;
vector<int>  RechTrackID_;
vector<float>  RechTime_;
vector<float>  RechDepth_;
vector<int>  RechID_;

// initialize branches
void Phase2::branchesRecHit(TTree* tree) {
  tree->Branch("RechEnergy"                   , &RechEnergy_);
  tree->Branch("RechMIP"                   , &RechMIP_);
  tree->Branch("RechTrackID"                   , &RechTrackID_);
  tree->Branch("RechTime"                   , &RechTime_);
  tree->Branch("RechDepth"                   , &RechDepth_);
  tree->Branch("RechID"                   , &RechID_);

}

//fills slimmedJets .clear() to empty vector of old data
void Phase2::fillRecHit(const edm::Event& e, const edm::EventSetup& es) {
 RechEnergy_.clear();
 RechMIP_.clear();
 RechTrackID_.clear();
 RechTime_.clear();
 RechDepth_.clear();
 RechID_.clear();
 e.getByToken( RechHitLabel_       ,  RechHitHandle );
 //rhtools_.getEventSetup(es);
 int maxlayer_ = rhtools_.lastLayerBH();
 // cleanup from previous execution


 for (std::vector<PCaloHit>::const_iterator rech = RechHitHandle->begin(); rech != RechHitHandle->end(); ++rech) {
  Float_t tempe = rech->energy();
  Float_t tempmip = tempe*10000.;
  Float_t tempt = rech->time();
  Int_t tempTrackID = rech->geantTrackId();
  Float_t tempd = rech->depth();
  uint32_t detId = rech->id();
  RechEnergy_.push_back(tempe);
  RechMIP_.push_back(tempmip);
  RechTime_.push_back(tempt);
  RechTrackID_.push_back(tempTrackID);
  RechDepth_.push_back(tempd);
  RechID_.push_back(detId);
}
}
