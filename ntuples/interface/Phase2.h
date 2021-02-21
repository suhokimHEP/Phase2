#ifndef Phase2_h
#define Phase2_h

#include "TTree.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

//#include "RecoTracker/DebugTools/interface/GetTrackTrajInfo.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "MagneticField/Engine/interface/MagneticField.h" 

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"



using namespace std;

void setbit(UShort_t& x, UShort_t bit);

class Phase2 : public edm::EDAnalyzer {
 public:

  explicit Phase2(const edm::ParameterSet&);
  ~Phase2();
  
  //   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  
  edm::ParameterSet lldj_pset_;

  //   virtual void beginJob() {};
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);//for trigger
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //   virtual void endJob() {};
  
  void branchesGlobalEvent (TTree*);
  void branchesMuons       (TTree*);
  void branchesTrigger     (TTree*);
  void branchesGenPart     (TTree*);
  void branchesTracks        (TTree*);
  void branchesRecHit       (TTree*);

  void fillGlobalEvent (const edm::Event&, const edm::EventSetup&, string Mode);
  void fillMuons       (const edm::Event&, const reco::Vertex, string Mode);
  void fillTrigger     (const edm::Event&, const edm::EventSetup&);
  void fillGenPart     (const edm::Event&);
  void fillTracks        (const edm::Event&, const edm::EventSetup&, string Mode);
  void fillHGCalHit   (const edm::Event&, const edm::EventSetup&, string Mode, string HGCMode);

  bool doAOD_     ; 
  string doMiniAOD_ ; 
  string HGCMode_ ; 

  // collections
  // met
    edm::EDGetTokenT<edm::TriggerResults>            patTrgResultsLabel_;

  // global event
  edm::EDGetTokenT<double>                         rhoCentralLabel_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     slimpuCollection_;
  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  edm::EDGetTokenT<reco::VertexCollection>         slimvtxLabel_;
  edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
  string                                           trgResultsProcess_;

  // beamspot
  edm::EDGetTokenT<reco::BeamSpot>                 beamspotLabel_;

  edm::EDGetTokenT<std::vector<reco::Track>  >       TrackLabel_;
  edm::EDGetTokenT<std::vector<pat::IsolatedTrack>  >       isoTrackLabel_;
  edm::EDGetTokenT<std::vector<PCaloHit>  >       EESimHitLabel_;
  edm::EDGetTokenT<std::vector<PCaloHit>  >       FHSimHitLabel_;
  edm::EDGetTokenT<std::vector<PCaloHit>  >       BHSimHitLabel_;
  edm::EDGetToken recHitSourceEE_;
  edm::EDGetToken recHitSourceFH_;
  edm::EDGetToken recHitSourceBH_;



  const MagneticField*                             magneticField_;
  edm::ESHandle<Propagator>                        thePropagator_;
  //edm::ESHandle<TransientTrackBuilder>             theBuilder_;

  // muons
  edm::EDGetTokenT<edm::View<reco::Muon> >          muonCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> >          slimmuonCollection_;

  //gen
  edm::EDGetTokenT<vector<reco::GenParticle> >     genParticlesCollection_;



  // trigger
  edm::EDGetTokenT<edm::TriggerResults>                     triggerBits_;
  edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>             triggerPrescales_;


  edm::InputTag AODTriggerLabel_;
  edm::EDGetTokenT<edm::TriggerResults> AODTriggerToken_;
  edm::InputTag AODTriggerEventLabel_;
  edm::EDGetTokenT<trigger::TriggerEvent> AODTriggerEventToken_;
  edm::Handle<edm::TriggerResults> AODTriggerHandle_;
  edm::Handle<trigger::TriggerEvent> AODTriggerEventHandle_;
  //HLTConfigProvider hltConfig_;
  //HLTPrescaleProvider hltPrescale_;


  TTree   *tree_;
  TH1F    *hEvents_;
  // shared between miniAOD jets and AOD jets modules
  edm::Handle<double>                 rhoHandle;
  edm::Handle<reco::VertexCollection> vtxHandle;
  vector<vector<float>> muEtaPhi_;
  // jet functions
  vector<int> getJetTrackIndexs( float jeteta, float jetphi);
  vector<float> CalTrackdR( float jeteta, float jetphi);
  void fill_simhit_tree_(const edm::Event& e, const edm::EventSetup&, const std::vector<PCaloHit>& hits);
  void fill_rechit_tree_(const edm::Event& e, const edm::EventSetup&, const HGCRecHitCollection& hits);













};

#endif
