#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Common/interface/TriggerNames.h"
#include "Phase2/ntuples/interface/Phase2.h"
#include <iomanip>
#include <bitset>
using namespace std;

// (local) variables associated with tree branches
Int_t       run_;
Long64_t    event_;
Int_t       lumis_;
Bool_t      isData_;
Int_t       nVtx_;
Int_t       nGoodVtx_;
Int_t       nTrksPV_;
Bool_t      isPVGood_;
float       vtx_;
float       vty_;
float       vtz_;
float       rho_;
float       rhoCentral_;

Int_t       nTruePU_;

void Phase2::branchesGlobalEvent(TTree* tree) {

  tree->Branch("run",     	       &run_);
  tree->Branch("event",    	       &event_);
  tree->Branch("lumis",   	       &lumis_);
  tree->Branch("isData",  	       &isData_);
  tree->Branch("nVtx",                 &nVtx_);
  tree->Branch("nGoodVtx",             &nGoodVtx_);
  tree->Branch("nTrksPV",              &nTrksPV_);
  tree->Branch("isPVGood",             &isPVGood_);
  tree->Branch("vtx",                  &vtx_);
  tree->Branch("vty",                  &vty_);
  tree->Branch("vtz",                  &vtz_);
  //tree->Branch("rho",                  &rho_);
  //tree->Branch("rhoCentral",           &rhoCentral_);


  tree->Branch("nTruePU",  &nTruePU_);

}

void Phase2::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  //phoPrescale_.clear();

  //edm::Handle<double> rhoHandle;
  //e.getByToken(rhoLabel_, rhoHandle);

  //edm::Handle<double> rhoCentralHandle;
  //e.getByToken(rhoCentralLabel_, rhoCentralHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  //rho_    = *(rhoHandle.product());
  //if (rhoCentralHandle.isValid()) rhoCentral_ = *(rhoCentralHandle.product());
  //else rhoCentral_ = -99.;
  //std::cout<<"here1"<<std::endl;
  nTruePU_ = -1 ;
  if (!e.isRealData()) {
   edm::Handle<vector<PileupSummaryInfo> > puInfoHandle;
   e.getByToken(puCollection_, puInfoHandle);
   if ( puInfoHandle->size() > 0 ){
    nTruePU_ = puInfoHandle->at(0).getTrueNumInteractions() ;
   }
  }
  //std::cout<<"here2"<<std::endl;

  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  nVtx_     = -1;
  nGoodVtx_ = -1;
  if (vtxHandle.isValid())
  {
   nVtx_     = 0;
   nGoodVtx_ = 0;
  //std::cout<<"here3"<<std::endl;

   for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v)
   {
    if (nVtx_ == 0)
    {
     nTrksPV_ = v->nTracks();
     vtx_     = v->x();
     vty_     = v->y();
     vtz_     = v->z();

     isPVGood_ = false;
     if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) isPVGood_ = true;
    }

    if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
    nVtx_++;
   }
  }
  else {edm::LogWarning("Phase2") << "Primary vertices info not unavailable";}

}
