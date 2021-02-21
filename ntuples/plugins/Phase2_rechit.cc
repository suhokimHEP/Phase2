//#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Phase2/ntuples/interface/Phase2.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

#include "DataFormats/Math/interface/deltaR.h"
using namespace std;
using namespace hgcal;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// AOD ---------------------------------------------
vector<float>  SimhitEnergy_;
vector<float>  SimhitEMEnergy_;
vector<float>  SimhitHadEnergy_;
vector<float>  SimhitMIP_;
vector<int>  SimhitTrackID_;
vector<int>  SimhitID_;
vector<float>  SimhitTime_;
vector<int>  SimhitDepth_;
vector<float>  SimhitGPx_;
vector<float>  SimhitGPy_;
vector<float>  SimhitGPr_;
vector<float>  SimhitGPz_;

vector<int>  SimhitLayer_;

vector<float>  HGCRechitEnergy_;
vector<DetId>  HGCRechitID_;
vector<GlobalPoint>  HGCRechitGP_;
vector<float>  HGCRechitEta_;
vector<float>  HGCRechitPhi_;
vector<int>  HGCRechitLayer_;
vector<int>  HGCRechitSithick_;

vector<float>  muHGCRHEn_;
vector<float>  muHGCRHdR_;
vector<float>  muHGCRHmuEta_;
// initialize branches
void Phase2::branchesRecHit(TTree* tree) {
  tree->Branch("SimhitEnergy"                   , &SimhitEnergy_);
  tree->Branch("SimhitEMEnergy"                   , &SimhitEMEnergy_);
  tree->Branch("SimhitHadEnergy"                   , &SimhitHadEnergy_);
  tree->Branch("SimhitMIP"                   , &SimhitMIP_);
  tree->Branch("SimhitTrackID"                   , &SimhitTrackID_);
  tree->Branch("SimhitID"                   , &SimhitID_);
  tree->Branch("SimhitTime"                   , &SimhitTime_);
  tree->Branch("SimhitDepth"                   , &SimhitDepth_);
  tree->Branch("SimhitGPx"                   , &SimhitGPx_);
  tree->Branch("SimhitGPy"                   , &SimhitGPy_);
  tree->Branch("SimhitGPr"                   , &SimhitGPr_);
  tree->Branch("SimhitGPz"                   , &SimhitGPz_);
  tree->Branch("SimhitLayer"                   , &SimhitLayer_);
  tree->Branch("HGCRechitEnergy"                   , &HGCRechitEnergy_);
  tree->Branch("HGCRechitID"                   , &HGCRechitID_);
  tree->Branch("HGCRechitGP"                   , &HGCRechitGP_);
  tree->Branch("HGCRechitEta"                   , &HGCRechitEta_);
  tree->Branch("HGCRechitPhi"                   , &HGCRechitPhi_);
  tree->Branch("HGCRechitLayer"                   , &HGCRechitLayer_);
  tree->Branch("HGCRechitSithick"                   , &HGCRechitSithick_);
  tree->Branch("muHGCRHEn"                   , &muHGCRHEn_);
  tree->Branch("muHGCRHdR"                   , &muHGCRHdR_);
  tree->Branch("muHGCRHmuEta"                   , &muHGCRHmuEta_);
}

void Phase2::fillHGCalHit(const edm::Event& e, const edm::EventSetup& es, string Mode, string HGCMode) {
TString Whichendcap;
Whichendcap = HGCMode;

edm::Handle<std::vector<PCaloHit>   >  EESimHitHandle;
edm::Handle<std::vector<PCaloHit>   >  FHSimHitHandle;
edm::Handle<std::vector<PCaloHit>   >  BHSimHitHandle;
edm::Handle<HGCRecHitCollection> handleTheRecHitsEE;
edm::Handle<HGCRecHitCollection> handleTheRecHitsFH;
edm::Handle<HGCRecHitCollection> handleTheRecHitsBH;
 e.getByToken( EESimHitLabel_       ,  EESimHitHandle );
 e.getByToken( FHSimHitLabel_       ,  FHSimHitHandle );
 e.getByToken( BHSimHitLabel_       ,  BHSimHitHandle );

  if(doMiniAOD_.find("reco")!= std::string::npos){
e.getByToken(recHitSourceEE_, handleTheRecHitsEE);
e.getByToken(recHitSourceFH_, handleTheRecHitsFH);
e.getByToken(recHitSourceBH_, handleTheRecHitsBH);
}

if(Whichendcap.Contains("EE"))
{
   fill_simhit_tree_(e,es, *EESimHitHandle);
}
else if(Whichendcap.Contains("FH"))

{
   fill_simhit_tree_(e,es, *FHSimHitHandle);

}
else{
   fill_simhit_tree_(e,es, *BHSimHitHandle);

}
  if(doMiniAOD_.find("reco")!= std::string::npos){
if(Whichendcap.Contains("EE"))
{
   fill_rechit_tree_(e,es, *handleTheRecHitsEE);
}
else if(Whichendcap.Contains("FH"))

{
   fill_rechit_tree_(e,es, *handleTheRecHitsFH);

}
else{
   fill_rechit_tree_(e,es, *handleTheRecHitsBH);

}
}
}

void Phase2::fill_simhit_tree_(const edm::Event& event, const edm::EventSetup& es ,const std::vector<PCaloHit>& hits) {
 SimhitEnergy_.clear();
 SimhitEMEnergy_.clear();
 SimhitHadEnergy_.clear();
 SimhitMIP_.clear();
 SimhitTrackID_.clear();
 SimhitID_.clear();
 SimhitTime_.clear();
 SimhitDepth_.clear();
 SimhitGPx_.clear();
 SimhitGPy_.clear();
 SimhitGPr_.clear();
 SimhitGPz_.clear();
 SimhitLayer_.clear();
    edm::ESHandle<CaloGeometry> sgeom;
    es.get<CaloGeometryRecord>().get(sgeom);

  for (unsigned int i = 0; i < hits.size(); ++i) {

    const PCaloHit& hit = hits[i];
  Float_t tempe = hit.energy();
  Float_t tempEM = hit.energyEM();
  Float_t tempHad = hit.energyHad();
  Float_t tempmip = tempe*10000.;
  Float_t tempt = hit.time();
  Int_t tempTrackID = hit.geantTrackId();
   uint32_t detId = hit.id();
   DetId detID = DetId(detId);
  uint16_t tempd = hit.depth();
   //std::cout<<tempd<<std::endl;
  SimhitEnergy_.push_back(tempe);
  SimhitEMEnergy_.push_back(tempEM);
  SimhitHadEnergy_.push_back(tempHad);
  SimhitMIP_.push_back(tempmip);
  SimhitTime_.push_back(tempt);
  SimhitTrackID_.push_back(tempTrackID);
  SimhitID_.push_back(detId);
  SimhitDepth_.push_back(tempd);
    GlobalPoint Position = sgeom->getGeometry(detId)->getPosition();
    float px = Position.x();
    float py = Position.y();
    float pz = Position.z();
    float radius = sqrt(px*px+py*py);
  SimhitGPx_.push_back(px);
  SimhitGPy_.push_back(py);
  SimhitGPr_.push_back(radius);
  SimhitGPz_.push_back(pz);
    unsigned int layer = std::numeric_limits<unsigned int>::max();
    if (detID.det() == DetId::HGCalEE || detID.det() == DetId::HGCalHSi) {
    layer = HGCSiliconDetId(detID).layer();
         } 
    else if (detID.det() == DetId::HGCalHSc) {
    layer = HGCScintillatorDetId(detID).layer();
  	}

   // HGCScintillatorDetId tempnum = HGCScintillatorDetId(detId).layer(); 
   // HGCSiliconDetId tempnum = HGCSiliconDetId(detId).layer(); 
   // int templay = int(tempnum);
  SimhitLayer_.push_back(layer);
}

}

void Phase2::fill_rechit_tree_(const edm::Event& event, const edm::EventSetup& es ,const HGCRecHitCollection& hits) {

 HGCRechitEnergy_.clear();
 HGCRechitID_.clear();
 HGCRechitGP_.clear();
 HGCRechitEta_.clear();
 HGCRechitPhi_.clear();
 HGCRechitLayer_.clear();
 HGCRechitSithick_.clear();
 muHGCRHEn_.clear();
 muHGCRHdR_.clear();
 muHGCRHmuEta_.clear();
    edm::ESHandle<CaloGeometry> geom;
    es.get<CaloGeometryRecord>().get(geom);
    //es.get<IdealGeometryRecord>().get(geom);
    RecHitTools rhtools_;
    rhtools_.setGeometry(*geom);
  for (unsigned int i = 0; i < hits.size(); ++i) {

    const HGCRecHit& hit = hits[i];
    float Energy = hit.energy();
    DetId detId = hit.detid();
    //GlobalPoint Position = geom->getGeometry(detId)->getPosition();
    GlobalPoint Position = rhtools_.getPosition(detId);
    float eta = Position.eta();
    float phi = Position.phi();


    HGCRechitEnergy_.push_back(Energy);
    HGCRechitID_.push_back(detId);
    HGCRechitGP_.push_back(Position);   
    HGCRechitEta_.push_back(eta);
    HGCRechitPhi_.push_back(phi);
    //HGCRechitLayer_.push_back(layer);   

    int ithickness = rhtools_.getSiThickIndex(detId);
    HGCRechitSithick_.push_back(ithickness);   
    int ilayer = rhtools_.getLayerWithOffset(detId);
    HGCRechitLayer_.push_back(ilayer);   

  for (unsigned int j = 0; j < muEtaPhi_.size(); ++j) {
    float muEta = muEtaPhi_.at(j).at(0);
    float muPhi = muEtaPhi_.at(j).at(1);
    float drt = deltaR( eta, phi, muEta, muPhi );
//    if(drt<.1){
muHGCRHEn_.push_back(Energy);	
muHGCRHdR_.push_back(drt);	
muHGCRHmuEta_.push_back(muEta);	
//	}
	}






    }
  
}

