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

    //Save long lived BSM particles
    /*if( abs(ip->pdgId()) == 9000006 ){
      llpId.push_back(      ip->pdgId() );
      llpStatus.push_back(  ip->status() );
      llpPt.push_back(      ip->pt()    );
      llpEta.push_back(     ip->eta()   );
      llpPhi.push_back(     ip->phi()   );
      llpMass.push_back(    ip->mass()  );
      llpvX.push_back(    ip->vx()  );
      llpvY.push_back(    ip->vy()  );
      llpvZ.push_back(    ip->vz()  );
      for(size_t j=0; j<ip->numberOfDaughters(); ++j){
	const reco::Candidate* d = ip->daughter(j);
	  llpDaughtervX.push_back(d->vx());
	  llpDaughtervY.push_back(d->vy());
	  llpDaughtervZ.push_back(d->vz());
	} 

    }*/
    /*else if ( particleHistory.hasRealParent() ) {
      reco::GenParticleRef momRef = particleHistory.parent();
      if ( momRef.isNonnull() && momRef.isAvailable() ) {
	if( abs(momRef->pdgId()) == 9000006 ){
	  llpDaughterId.push_back(     ip->pdgId() );
	  llpDaughterStatus.push_back( ip->status() );
	  llpDaughterPt.push_back(     ip->pt()    );
	  llpDaughterEta.push_back(    ip->eta()   );
	  llpDaughterPhi.push_back(    ip->phi()   );
	  llpDaughterMass.push_back(   ip->mass()  );
	}
      }
    }*/

  }//end gen loop

}
