#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Phase2/ntuples/interface/Phase2.h"


#include "SimDataFormats/CaloHit/interface/PCaloHit.h"



using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
edm::Handle<std::vector<PCaloHit>   >  RechHitHandle;

// AOD ---------------------------------------------
vector<float>  RechEnergy_;
vector<int>  RechID_;
vector<float>  RechTime_;

// initialize branches
void Phase2::branchesRecHit(TTree* tree) {
  tree->Branch("RechEnergy"                   , &RechEnergy_);
  tree->Branch("RechID"                   , &RechID_);
  tree->Branch("RechTime"                   , &RechTime_);

}

//fills slimmedJets .clear() to empty vector of old data
void Phase2::fillRecHit(const edm::Event& e, const edm::EventSetup& es) {
 RechEnergy_.clear();
 RechID_.clear();
 RechTime_.clear();
 e.getByToken( RechHitLabel_       ,  RechHitHandle );

 // bool dodebug = false;
 // cleanup from previous execution


 //for(int j = 0; j < (int)TrackHandle->size(); j++){
 for (std::vector<PCaloHit>::const_iterator rech = RechHitHandle->begin(); rech != RechHitHandle->end(); ++rech) {
  Float_t tempe = rech->energy();
  //std::cout<<"eta:"<<eta<<std::endl;
  Float_t tempt = rech->time();
  //std::cout<<"pfNeutralSum:"<<pfNeutralSum<<std::endl;
  Int_t tempID = rech->geantTrackId();
  //std::cout<<"dEdxPixel:"<<dEdxPixel<<std::endl;
  RechEnergy_.push_back(tempe);
  RechTime_.push_back(tempt);
  RechID_.push_back(tempID);
}
}
