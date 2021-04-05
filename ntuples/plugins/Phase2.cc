#include "Phase2/ntuples/interface/Phase2.h"
//#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
//#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

//Phase2::Phase2(const edm::ParameterSet& ps) :
//  hltPrescale_(ps,consumesCollector(),*this) 
Phase2::Phase2(const edm::ParameterSet& ps)
{  
  lldj_pset_ = ps;
  doMiniAOD_               = ps.getParameter<string>("doMiniAOD");
  HGCMode_               = ps.getParameter<string>("HGCMode");
  

  trgResultsLabel_         = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("triggerResults"));
  trgResultsProcess_       =                                          ps.getParameter<InputTag>("triggerResults").process();
  //// met
  patTrgResultsLabel_      = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("patTriggerResults"));
  isoTrackLabel_           = consumes<std::vector<pat::IsolatedTrack> >       (ps.getParameter<InputTag>("isoTrackSrc"));
  TrackLabel_           = consumes<std::vector<reco::Track> >       (ps.getParameter<InputTag>("TrackSrc"));
  EESimHitLabel_           = consumes<std::vector<PCaloHit> >       (ps.getParameter<InputTag>("EESimHitSrc"));
  FHSimHitLabel_           = consumes<std::vector<PCaloHit> >       (ps.getParameter<InputTag>("FHSimHitSrc"));
  BHSimHitLabel_           = consumes<std::vector<PCaloHit> >       (ps.getParameter<InputTag>("BHSimHitSrc"));
  recHitSourceEE_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("EERecHits"));
  recHitSourceFH_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("FHRecHits"));
  recHitSourceBH_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("BHRecHits"));
  // muons
  slimmuonCollection_          = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("slimmuonSrc"));
  muonCollection_          = consumes<View<reco::Muon> >              (ps.getParameter<InputTag>("muonSrc"));
  // trigger
  triggerBits_                    = consumes <edm::TriggerResults>                     (ps.getParameter<edm::InputTag>("bits"));
  triggerObjects_                 = consumes <edm::View<pat::TriggerObjectStandAlone>> (ps.getParameter<edm::InputTag>("objects"));
  triggerPrescales_               = consumes <pat::PackedTriggerPrescales>             (ps.getParameter<edm::InputTag>("prescales"));

  puCollection_            = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
  slimpuCollection_            = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("slimpileupCollection"));
  vtxLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));
  slimvtxLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("slimVtxLabel"));


  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));
  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed events",   1,  0,   2);

  branchesGlobalEvent(tree_);
  branchesGenPart(tree_);
  branchesMuons(tree_);
  branchesTracks(tree_);
  branchesRecHit(tree_);
}

Phase2::~Phase2() {
}


void Phase2::beginRun(edm::Run const& run, edm::EventSetup const& eventsetup) {

  //if(hltConfig_.init(run,eventsetup,"HLT",changed)){
  //}
  //hltPrescale_.init(run,eventsetup,"HLT",changed);

}

void Phase2::analyze(const edm::Event& e, const edm::EventSetup& es) {
  if(doMiniAOD_.find("reco")!= std::string::npos){
  muEtaPhi_.clear();
   fillGlobalEvent(e, es,doMiniAOD_);
  fillGenPart(e);
  fillTracks(e, es,doMiniAOD_);
  // muons use vtx for isolation
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  reco::Vertex vtx;
  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0); 
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    bool isFake = -(v->chi2() == 0 && v->ndof() == 0); 
 
    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v; 
      break;
    }   
  }
   fillMuons(e,es, vtx,doMiniAOD_); //muons use vtx for isolation
   fillHGCalHit(e,es,doMiniAOD_,HGCMode_);
}
  else if(doMiniAOD_.find("RAW")!= std::string::npos){
   fillGlobalEvent(e, es,doMiniAOD_);
   fillGenPart(e);
   fillHGCalHit(e,es,doMiniAOD_,HGCMode_);
}
 else {
   fillGlobalEvent(e, es,doMiniAOD_);
  fillGenPart(e);

  fillTracks(e, es,doMiniAOD_);
  // muons use vtx for isolation
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  reco::Vertex vtx;
  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0); 
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    bool isFake = -(v->chi2() == 0 && v->ndof() == 0); 
 
    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v; 
      break;
    }   
  }
   fillMuons(e, es,vtx,doMiniAOD_); //muons use vtx for isolation

}
 hEvents_->Fill(1.);
 tree_->Fill();
}

DEFINE_FWK_MODULE(Phase2);
