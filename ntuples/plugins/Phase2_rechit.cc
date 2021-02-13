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
vector<float>  SimhitMIP_;
vector<int>  SimhitTrackID_;
vector<float>  SimhitTime_;
vector<float>  SimhitDepth_;
vector<int>  SimhitID_;
vector<float>  HGCRechitEnergy_;
vector<DetId>  HGCRechitID_;
vector<GlobalPoint>  HGCRechitGP_;
vector<float>  HGCRechitEta_;
vector<float>  HGCRechitPhi_;
vector<int>  HGCRechitLayer_;

vector<float>  muHGCRHEn_;
vector<float>  muHGCRHdR_;
vector<float>  muHGCRHmuEta_;
// initialize branches
void Phase2::branchesRecHit(TTree* tree) {
  tree->Branch("SimhitEnergy"                   , &SimhitEnergy_);
  tree->Branch("SimhitMIP"                   , &SimhitMIP_);
  tree->Branch("SimhitTrackID"                   , &SimhitTrackID_);
  tree->Branch("SimhitTime"                   , &SimhitTime_);
  tree->Branch("SimhitDepth"                   , &SimhitDepth_);
  tree->Branch("SimhitID"                   , &SimhitID_);
  tree->Branch("HGCRechitEnergy"                   , &HGCRechitEnergy_);
  tree->Branch("HGCRechitID"                   , &HGCRechitID_);
  tree->Branch("HGCRechitGP"                   , &HGCRechitGP_);
  tree->Branch("HGCRechitEta"                   , &HGCRechitEta_);
  tree->Branch("HGCRechitPhi"                   , &HGCRechitPhi_);
  tree->Branch("HGCRechitLayer"                   , &HGCRechitLayer_);
  tree->Branch("muHGCRHEn"                   , &muHGCRHEn_);
  tree->Branch("muHGCRHdR"                   , &muHGCRHdR_);
  tree->Branch("muHGCRHmuEta"                   , &muHGCRHmuEta_);
}

void Phase2::fillRecHit(const edm::Event& e, const edm::EventSetup& es, string Mode) {
TString Whichendcap;
Whichendcap = "BH";

 SimhitEnergy_.clear();
 SimhitMIP_.clear();
 SimhitTrackID_.clear();
 SimhitTime_.clear();
 SimhitDepth_.clear();
 SimhitID_.clear();
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
    edm::ESHandle<CaloGeometry> geom;
    es.get<CaloGeometryRecord>().get(geom);
     

    //rhtools_.setGeometry(*geom);
    //edm::ESHandle<HGCalGeometry> geom;
    //es.get<CaloGeometryRecord>().get(geom);
    //rhtools_.setGeometry(*geom);
 //rhtools_.getEventSetup(es);
// int maxlayer_ = -1;
//maxlayer_= rhtools_.lastLayerBH();
  //std::cout<<maxlayer_<<std::endl;
 // cleanup from previous execution

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
 SimhitMIP_.clear();
 SimhitTrackID_.clear();
 SimhitTime_.clear();
 SimhitDepth_.clear();
 SimhitID_.clear();

  for (unsigned int i = 0; i < hits.size(); ++i) {

    const PCaloHit& hit = hits[i];
  Float_t tempe = hit.energy();
  Float_t tempmip = tempe*10000.;
  Float_t tempt = hit.time();
  Int_t tempTrackID = hit.geantTrackId();
  Float_t tempd = hit.depth();
  uint32_t detId = hit.id();
  SimhitEnergy_.push_back(tempe);
  SimhitMIP_.push_back(tempmip);
  SimhitTime_.push_back(tempt);
  SimhitTrackID_.push_back(tempTrackID);
  SimhitDepth_.push_back(tempd);
  SimhitID_.push_back(detId);
}

}

void Phase2::fill_rechit_tree_(const edm::Event& event, const edm::EventSetup& es ,const HGCRecHitCollection& hits) {

 HGCRechitEnergy_.clear();
 HGCRechitID_.clear();
 HGCRechitGP_.clear();
 HGCRechitEta_.clear();
 HGCRechitPhi_.clear();
 HGCRechitLayer_.clear();
 muHGCRHEn_.clear();
 muHGCRHdR_.clear();
 muHGCRHmuEta_.clear();
    edm::ESHandle<CaloGeometry> geom;
    es.get<CaloGeometryRecord>().get(geom);
    RecHitTools rhtools_;

    //es.get<IdealGeometryRecord>().get(geom);
    //rhtools_.setGeometry(*geom);
  //std::cout<<"FH rhsize:"<<hits.size()<<std::endl;
  for (unsigned int i = 0; i < hits.size(); ++i) {

    const HGCRecHit& hit = hits[i];
    float Energy = hit.energy();
    DetId detId = hit.detid();
    GlobalPoint Position = geom->getGeometry(detId)->getPosition();
    float eta = Position.eta();
    float phi = Position.phi();
    unsigned int layer = std::numeric_limits<unsigned int>::max();
   if (detId.det() == DetId::HGCalEE || detId.det() == DetId::HGCalHSi) {
    layer = HGCSiliconDetId(detId).layer();
  } else if (detId.det() == DetId::HGCalHSc) {
    layer = HGCScintillatorDetId(detId).layer();
  } 
   // HGCScintillatorDetId tempnum = HGCScintillatorDetId(detId).layer(); 
   //   HGCSiliconDetId tempnum = HGCSiliconDetId(detId).layer(); 
   // int templay = int(tempnum);

    HGCRechitEnergy_.push_back(Energy);
    HGCRechitID_.push_back(detId);
    HGCRechitGP_.push_back(Position);   
    HGCRechitEta_.push_back(eta);
    HGCRechitPhi_.push_back(phi);
    HGCRechitLayer_.push_back(layer);   

  for (unsigned int j = 0; j < muEtaPhi_.size(); ++j) {
    float muEta = muEtaPhi_.at(j).at(0);
    float muPhi = muEtaPhi_.at(j).at(1);
    float drt = deltaR( eta, phi, muEta, muPhi );
    //if(drt<99.){
muHGCRHEn_.push_back(Energy);	
muHGCRHdR_.push_back(drt);	
muHGCRHmuEta_.push_back(muEta);	

//	}
	}

	//int ScintId = tempScintId;
    //std::cout<<detId<<std::endl;    
    int ithickness = rhtools_.getSiThickIndex(detId);
    //std::cout<<ithickness<<std::endl;
    int ilayer = rhtools_.getLayerWithOffset(detId)-1;
    //
    //if (ithickness == -1)
    //    ithickness = 3;
    //double storedThreshold = thresholds_[ilayer][ithickness];
    //if (hit.energy() < storedThreshold) continue;
    
    //rechitsInfo.event = event;
    //rechitsInfo.energy = hit.energy();
    //rechitsInfo.time = hit.time();
    //GlobalPoint global = geom->getPosition(detId);
    //GlobalPoint global = rhtools_.getPosition(detId);
    //rechitsInfo.x = global.x();
    //rechitsInfo.y = global.y();
    //rechitsInfo.z = global.z();
    //rechitsInfo.eta = global.eta();
    //rechitsInfo.phi = global.phi();
    //rechitsInfo.layer = ilayer+1;
    //rechitsInfo.additional1 = ithickness;
    //RecHitTree->Fill();
    //if (regionEtaMin_ < global.eta() &&
    //    global.eta() < regionEtaMax_ &&
    //    regionPhiMin_ < global.phi() &&
    //    global.phi() < regionPhiMax_) {
    //  rechitsSignalRegionInfo.event = event;
    //  rechitsSignalRegionInfo.energy = hit.energy();
    //  rechitsSignalRegionInfo.time = hit.time();
    //  //GlobalPoint global = geom->getPosition(detId);
    //  //GlobalPoint global = rhtools_.getPosition(detId);
    //  rechitsSignalRegionInfo.x = global.x();
    //  rechitsSignalRegionInfo.y = global.y();
    //  rechitsSignalRegionInfo.z = global.z();
    //  rechitsSignalRegionInfo.eta = global.eta();
    //  rechitsSignalRegionInfo.phi = global.phi();
    //  rechitsSignalRegionInfo.layer = ilayer+1;
    //  rechitsSignalRegionInfo.additional1 = ithickness;
    //  RecHitSignalRegionTree->Fill();
    }
  
}

