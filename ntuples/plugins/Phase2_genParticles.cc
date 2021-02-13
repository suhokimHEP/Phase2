#include "Phase2/ntuples/interface/Phase2.h"
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;


//Variables for branches
vector<int>   llpId;
vector<int>   llpStatus;
vector<float> llpPt;
vector<float> llpEta;
vector<float> llpPhi;
vector<float> llpMass;
vector<int>   llpDaughterId;
vector<int>   llpDaughterStatus;
vector<float> llpDaughterPt;
vector<float> llpDaughterEta;
vector<float> llpDaughterPhi;
vector<float> llpDaughterMass;

vector<float> Wpt;
vector<float> Wmass;
vector<vector<int>>   W_daughterID;
vector<vector<int>>   W_daughterPt;
vector<vector<int>>   W_daughterEta;
vector<vector<int>>   W_daughterPhi;
vector<float> Genmupt;
vector<float> Genmueta;
vector<float> Genmuenergy;

vector<float> llpvX;
vector<float> llpvY;
vector<float> llpvZ;
vector<float> llpDaughtervX;
vector<float> llpDaughtervY;
vector<float> llpDaughtervZ;
void Phase2::branchesGenPart(TTree* tree) {

  //tree->Branch("llpId",             &llpId);
  //tree->Branch("llpStatus",         &llpStatus);
  //tree->Branch("llpPt",             &llpPt);
  //tree->Branch("llpEta",            &llpEta);
  //tree->Branch("llpPhi",            &llpPhi);
  //tree->Branch("llpMass",           &llpMass);
  //tree->Branch("llpDaughterId",     &llpDaughterId);
  //tree->Branch("llpDaughterStatus", &llpDaughterStatus);
  //tree->Branch("llpDaughterPt",     &llpDaughterPt);
  //tree->Branch("llpDaughterEta",    &llpDaughterEta);
  //tree->Branch("llpDaughterPhi",    &llpDaughterPhi);
  //tree->Branch("llpDaughterMass",   &llpDaughterMass);
  tree->Branch("Wpt",    &Wpt);
  tree->Branch("Wmass",  &Wmass);
  tree->Branch("W_daughterID",   &W_daughterID);
  tree->Branch("W_daughterPt",   &W_daughterPt);
  tree->Branch("W_daughterEta",  &W_daughterEta);
  tree->Branch("W_daughterPhi",  &W_daughterPhi);
  tree->Branch("Genmupt",    &Genmupt);
  tree->Branch("Genmueta",    &Genmueta);
  tree->Branch("Genmuenergy",    &Genmuenergy);

  //tree->Branch("llpvX",             &llpvX);
  //tree->Branch("llpvY",             &llpvY);
  //tree->Branch("llpvZ",             &llpvZ);
  //tree->Branch("llpDaughtervX",             &llpDaughtervX);
  //tree->Branch("llpDaughtervY",             &llpDaughtervY);
  //tree->Branch("llpDaughtervZ",             &llpDaughtervZ);

}

void Phase2::fillGenPart(const edm::Event& e) {

  //Initialize -- set numbers to e.g. 0 and clear vectors 
  llpId.clear();
  llpStatus.clear();
  llpPt.clear();
  llpEta.clear();
  llpPhi.clear();
  llpMass.clear();
  llpDaughterId.clear();
  llpDaughterStatus.clear();
  llpDaughterPt.clear();
  llpDaughterEta.clear();
  llpDaughterPhi.clear();
  llpDaughterMass.clear();
  Wpt.clear();
  Wmass.clear();
  W_daughterID.clear();
  W_daughterPt.clear();
  W_daughterEta.clear();
  W_daughterPhi.clear();
  Genmupt.clear();
  Genmueta.clear();
  Genmuenergy.clear();

  llpvX.clear();
  llpvY.clear();
  llpvZ.clear();
  llpDaughtervX.clear();
  llpDaughtervY.clear();
  llpDaughtervZ.clear();

  //Gen particles handle
  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  //Loop over gen particles
  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    if( abs(ip->pdgId()) == 13 && ip->isLastCopy()){
	Genmupt.push_back(ip->pt());
	Genmueta.push_back(ip->eta());
	Genmuenergy.push_back(ip->energy());




 	} 
    
    //Save W particles
    vector<int> W_daughterID_;
    vector<int> W_daughterPt_;
    vector<int> W_daughterEta_;
    vector<int> W_daughterPhi_;
    if( abs(ip->pdgId()) == 24 && ip->isLastCopy()){
     Wpt.push_back  ( ip->pt() );
     Wmass.push_back( ip->mass() );
     for(size_t i=0; i<ip->numberOfDaughters(); ++i){
       const reco::Candidate* W_d = ip->daughter(i);
       W_daughterID_.push_back(W_d->pdgId());
       W_daughterPt_.push_back(W_d->pt());
       W_daughterEta_.push_back(W_d->eta());
       W_daughterPhi_.push_back(W_d->phi());
     }
     W_daughterID.push_back(W_daughterID_);
     W_daughterPt.push_back(W_daughterPt_);
     W_daughterEta.push_back(W_daughterEta_);
     W_daughterPhi.push_back(W_daughterPhi_);
    }

  }//end gen loop

}
