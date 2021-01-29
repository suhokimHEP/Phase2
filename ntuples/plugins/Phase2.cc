#include "Phase2/ntuples/interface/Phase2.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

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
  doMiniAOD_               = ps.getParameter<bool>("doMiniAOD");


  trgResultsLabel_         = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("triggerResults"));
  trgResultsProcess_       =                                          ps.getParameter<InputTag>("triggerResults").process();
  //// met
  patTrgResultsLabel_      = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("patTriggerResults"));
  //BadChCandFilterToken_    = consumes<bool>                          (ps.getParameter<InputTag>("BadChargedCandidateFilter"));
  //BadPFMuonFilterToken_    = consumes<bool>                          (ps.getParameter<edm::InputTag>("BadPFMuonFilter"));
  //pfMETlabel_              = consumes<View<pat::MET> >               (ps.getParameter<InputTag>("pfMETLabel"));
  //AODCaloMETlabel_         = consumes<edm::View<reco::CaloMET> >     (ps.getParameter<InputTag>("AODCaloMETlabel"));
  //AODpfChMETlabel_         = consumes<edm::View<reco::PFMET> >       (ps.getParameter<InputTag>("AODpfChMETlabel"));
  //AODpfMETlabel_           = consumes<edm::View<reco::PFMET> >       (ps.getParameter<InputTag>("AODpfMETlabel"));
  //TrackLabel_           = consumes<edm::View<pat::IsolatedTrack> >       (ps.getParameter<InputTag>("TrackSrc"));
  TrackLabel_           = consumes<std::vector<pat::IsolatedTrack> >       (ps.getParameter<InputTag>("TrackSrc"));
  RechHitLabel_           = consumes<std::vector<PCaloHit> >       (ps.getParameter<InputTag>("RechHitSrc"));

  // muons
  muonCollection_          = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("muonSrc"));
  // trigger
  triggerBits_                    = consumes <edm::TriggerResults>                     (ps.getParameter<edm::InputTag>("bits"));
  triggerObjects_                 = consumes <edm::View<pat::TriggerObjectStandAlone>> (ps.getParameter<edm::InputTag>("objects"));
  triggerPrescales_               = consumes <pat::PackedTriggerPrescales>             (ps.getParameter<edm::InputTag>("prescales"));







  puCollection_            = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
  vtxLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));


  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));













  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  hTTSF_   = fs->make<TH1F>("hTTSF",      "TTbar scalefactors",   200,  0,   2);
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

  bool changed(true);
  //if(hltConfig_.init(run,eventsetup,"HLT",changed)){
  //}
  //hltPrescale_.init(run,eventsetup,"HLT",changed);

}

void Phase2::analyze(const edm::Event& e, const edm::EventSetup& es) {
  if(doMiniAOD_){fillRecHit(e,es);}
  else {fillGlobalEvent(e, es);
  // muons use vtx for isolation
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  reco::Vertex vtx;
  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0); 
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = -(v->chi2() == 0 && v->ndof() == 0); 
 
    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v; 
      break;
    }   
  }
  fillGenPart(e);
  fillTracks(e, es);
  fillMuons(e, vtx); //muons use vtx for isolation
}
 hEvents_->Fill(1.);
 tree_->Fill();
}

DEFINE_FWK_MODULE(Phase2);
