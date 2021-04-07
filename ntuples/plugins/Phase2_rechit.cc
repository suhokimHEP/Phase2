#include "Phase2/ntuples/interface/Phase2.h"
//#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HGCalGeometryRecord.h"

//#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/Math/interface/deltaR.h"
using namespace std;
using namespace hgcal;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// AOD ---------------------------------------------
vector<float>  SimhitEnergy_;
vector<float>  SimhitLogEnergy_;
vector<float>  SimhitSeparateEnergy_;
vector<float>  SimhitSeparateEMEnergy_;
vector<float>  SimhitSeparateHadEnergy_;
vector<int>  SimhitTrackID_;
vector<uint32_t>  SimhitID_;
vector<float>  SimhitTime_;
vector<int>  SimhitDepth_;
vector<float>  SimhitGPx_;
vector<float>  SimhitGPy_;
vector<float>  SimhitGPr_;
vector<float>  SimhitGPz_;
vector<int>  SimhitLayer_;


vector<pair<uint32_t,float>> IdAccEn_;


vector<float>  HGCRechitEnergy_;
vector<DetId>  HGCRechitID_;
vector<GlobalPoint>  HGCRechitGP_;
vector<float>  HGCRechitEta_;
vector<float>  HGCRechitPhi_;
vector<tuple<float,float,float,int,float,float>> HGCRechitEnEtaPhiLayer_;
vector<int>  HGCRechitLayer_;
vector<int>  HGCRechitSithick_;

vector<float>  muHGCRHEn_;
vector<int>  muHGCRHLayer_;
vector<float>  muHGCRHdR_;
vector<float>  DenmuHGCRHEn_;
vector<int>  DenmuHGCRHLayer_;
vector<float>  DenmuHGCRHEta_;
vector<float>  DenmuHGCRHX_;
vector<float>  DenmuHGCRHY_;
vector<float>  DenmuHGCRHR_;
vector<float>  p2muHGCRHEn_;
vector<int>  p2muHGCRHLayer_;
vector<float>  p2muHGCRHEta_;
vector<float>  p2muHGCRHX_;
vector<float>  p2muHGCRHY_;
vector<float>  p2muHGCRHR_;
vector<float>  p2muHGCRHPhi_;

vector<float>  muPropEta_;

vector<DetId>  muPropId_;
vector<DetId>  RechitId_;
vector<pair<uint32_t,float>> recIdEn_;
vector<pair<uint32_t,float>> recIdEta_;

vector<float> MatchedSimE_;
vector<float> MatchedRecE_;
vector<int> MatchedELayer_;

// initialize branches
void Phase2::branchesRecHit(TTree* tree) {
  tree->Branch("SimhitEnergy"                   , &SimhitEnergy_);
  tree->Branch("SimhitLogEnergy"                   , &SimhitLogEnergy_);
  tree->Branch("SimhitSeparateEnergy"                   , &SimhitSeparateEnergy_);
  tree->Branch("SimhitSeparateEMEnergy"                   , &SimhitSeparateEMEnergy_);
  tree->Branch("SimhitSeparateHadEnergy"                   , &SimhitSeparateHadEnergy_);
  tree->Branch("SimhitTrackID"                   , &SimhitTrackID_);
  tree->Branch("SimhitID"                   , &SimhitID_);
  tree->Branch("SimhitTime"                   , &SimhitTime_);
  tree->Branch("SimhitDepth"                   , &SimhitDepth_);
  tree->Branch("SimhitGPx"                   , &SimhitGPx_);
  tree->Branch("SimhitGPy"                   , &SimhitGPy_);
  tree->Branch("SimhitGPr"                   , &SimhitGPr_);
  tree->Branch("SimhitGPz"                   , &SimhitGPz_);
  tree->Branch("SimhitLayer"                   , &SimhitLayer_);
  //tree->Branch("IdAccEn"                   , &IdAccEn_);
  tree->Branch("HGCRechitEnergy"                   , &HGCRechitEnergy_);
  tree->Branch("HGCRechitID"                   , &HGCRechitID_);
  tree->Branch("HGCRechitGP"                   , &HGCRechitGP_);
  tree->Branch("HGCRechitEta"                   , &HGCRechitEta_);
  tree->Branch("HGCRechitPhi"                   , &HGCRechitPhi_);
  tree->Branch("HGCRechitLayer"                   , &HGCRechitLayer_);
  tree->Branch("HGCRechitSithick"                   , &HGCRechitSithick_);
  tree->Branch("muHGCRHEn"                   , &muHGCRHEn_);
  tree->Branch("muHGCRHLayer"                   , &muHGCRHLayer_);
  tree->Branch("muHGCRHdR"                   , &muHGCRHdR_);
  tree->Branch("DenmuHGCRHEn"                   , &DenmuHGCRHEn_);
  tree->Branch("DenmuHGCRHLayer"                   , &DenmuHGCRHLayer_);
  tree->Branch("DenmuHGCRHEta"                   , &DenmuHGCRHEta_);
  tree->Branch("DenmuHGCRHX"                   , &DenmuHGCRHX_);
  tree->Branch("DenmuHGCRHY"                   , &DenmuHGCRHY_);
  tree->Branch("DenmuHGCRHR"                   , &DenmuHGCRHR_);
  tree->Branch("p2muHGCRHEn"                   , &p2muHGCRHEn_);
  tree->Branch("p2muHGCRHLayer"                   , &p2muHGCRHLayer_);
  tree->Branch("p2muHGCRHEta"                   , &p2muHGCRHEta_);
  tree->Branch("p2muHGCRHX"                   , &p2muHGCRHX_);
  tree->Branch("p2muHGCRHY"                   , &p2muHGCRHY_);
  tree->Branch("p2muHGCRHR"                   , &p2muHGCRHR_);
  tree->Branch("p2muHGCRHPhi"                   , &p2muHGCRHPhi_);
  tree->Branch("muPropEta"                   , &muPropEta_);
  tree->Branch("MatchedSimE"                   , &MatchedSimE_);
  tree->Branch("MatchedRecE"                   , &MatchedRecE_);
  tree->Branch("MatchedELayer"                   , &MatchedELayer_);
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

SimRecMatch(IdAccEn_,recIdEn_,MatchedSimE_,MatchedRecE_,SimhitLayer_,MatchedELayer_);
}
}

void Phase2::fill_simhit_tree_(const edm::Event& event, const edm::EventSetup& es ,const std::vector<PCaloHit>& hits) {
 SimhitEnergy_.clear();
 SimhitLogEnergy_.clear();
 SimhitSeparateEnergy_.clear();
 SimhitSeparateEMEnergy_.clear();
 SimhitSeparateHadEnergy_.clear();
 SimhitTrackID_.clear();
 SimhitID_.clear();
 SimhitTime_.clear();
 SimhitDepth_.clear();
 SimhitGPx_.clear();
 SimhitGPy_.clear();
 SimhitGPr_.clear();
 SimhitGPz_.clear();
 SimhitLayer_.clear();
 IdAccEn_.clear();
edm::ESHandle<CaloGeometry> sgeom;
es.get<CaloGeometryRecord>().get(sgeom);
std::pair <uint32_t,float> IdEn;


  for (unsigned int i = 0; i < hits.size(); ++i) {

    const PCaloHit& hit = hits[i];
  Float_t tempe = hit.energy();
  Float_t tempEM = hit.energyEM();
  Float_t tempHad = hit.energyHad();
  Float_t tempt = hit.time();
  Int_t tempTrackID = hit.geantTrackId();
   uint32_t detId = hit.id();
  uint16_t tempd = hit.depth();
  SimhitSeparateEnergy_.push_back(tempe);
  SimhitSeparateEMEnergy_.push_back(tempEM);
  SimhitSeparateHadEnergy_.push_back(tempHad);
  SimhitTime_.push_back(tempt);
  SimhitTrackID_.push_back(tempTrackID);
  SimhitID_.push_back(detId);
  SimhitDepth_.push_back(tempd);
   IdEn={}; 
   IdEn = std::make_pair(detId,tempe);
    //IdAccEn_.insert(IdEn);
    IdAccEn_.push_back(IdEn);
}


for (unsigned int i = 1; i<IdAccEn_.size(); ++i)
{
  uint32_t prevIndex = IdAccEn_.at(i-1).first;
  uint32_t currIndex = IdAccEn_.at(i).first;
  float prevVal = IdAccEn_.at(i-1).second;
  float currVal = IdAccEn_.at(i).second;
   if(prevIndex==currIndex)
	{IdAccEn_.at(i-1).second = prevVal+currVal;
	IdAccEn_.erase(IdAccEn_.begin()+i);
	i = i-1;}
}

  for (unsigned int i = 0; i < IdAccEn_.size(); ++i) {

   uint32_t AccdetId = IdAccEn_.at(i).first;
   DetId AccdetID = DetId(AccdetId);
   float AccEn = IdAccEn_.at(i).second;
   SimhitEnergy_.push_back(AccEn);
   SimhitLogEnergy_.push_back(log10(AccEn));
   
    GlobalPoint Position = sgeom->getGeometry(AccdetId)->getPosition();
    float px = Position.x();
    float py = Position.y();
    float pz = Position.z();
    float radius = sqrt(px*px+py*py);
  SimhitGPx_.push_back(px);
  SimhitGPy_.push_back(py);
  SimhitGPr_.push_back(radius);
  SimhitGPz_.push_back(pz);
    unsigned int layer = std::numeric_limits<unsigned int>::max();
    if (AccdetID.det() == DetId::HGCalEE || AccdetID.det() == DetId::HGCalHSi) {
    layer = HGCSiliconDetId(AccdetID).layer();
         } 
    else if (AccdetID.det() == DetId::HGCalHSc) {
    layer = HGCScintillatorDetId(AccdetID).layer();
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
 HGCRechitEnEtaPhiLayer_.clear();
 HGCRechitLayer_.clear();
 HGCRechitSithick_.clear();
 muHGCRHEn_.clear();
 muHGCRHLayer_.clear();
 muHGCRHdR_.clear();
 DenmuHGCRHEn_.clear();
 DenmuHGCRHLayer_.clear();
 DenmuHGCRHEta_.clear();
 DenmuHGCRHX_.clear();
 DenmuHGCRHY_.clear();
 DenmuHGCRHR_.clear();
 p2muHGCRHEn_.clear();
 p2muHGCRHLayer_.clear();
 p2muHGCRHEta_.clear();
 p2muHGCRHX_.clear();
 p2muHGCRHY_.clear();
 p2muHGCRHR_.clear();
 p2muHGCRHPhi_.clear();
 muPropEta_.clear();
 muPropId_.clear();
 RechitId_.clear();
 recIdEn_.clear();
 recIdEta_.clear();
 std::pair <uint32_t,float> Idvar;

    edm::ESHandle<CaloGeometry> geom;
    es.get<CaloGeometryRecord>().get(geom);
    edm::ESHandle<HGCalGeometry> BHgeom;
    es.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",BHgeom);
    const HGCalGeometry* geom0BH = BHgeom.product();
    DetId mudetID;
    RecHitTools rhtools_;
    rhtools_.setGeometry(*geom);
  for (unsigned int i = 0; i < hits.size(); ++i) {

    const HGCRecHit& hit = hits[i];
    float Energy = hit.energy();
    DetId detId = hit.detid();
    //std::cout<<"rechitID:"<<detId.rawId()<<",rechitEnergy:"<<Energy<<std::endl;
    RechitId_.push_back(detId);
    //GlobalPoint Position = geom->getGeometry(detId)->getPosition();
    //const Surface sur = geom->surface();
    GlobalPoint Position = rhtools_.getPosition(detId);
    //std::cout<<"rechitPosition:"<<Position<<std::endl;
    float eta = Position.eta();
    float phi = Position.phi();
    float GPx = Position.x();
    float GPy = Position.y();
    HGCRechitEnergy_.push_back(Energy);
    HGCRechitID_.push_back(detId);
    HGCRechitGP_.push_back(Position);   
    HGCRechitEta_.push_back(eta);
    HGCRechitPhi_.push_back(phi);

    int ithickness = rhtools_.getSiThickIndex(detId);
    HGCRechitSithick_.push_back(ithickness);   
    int ilayer = rhtools_.getLayerWithOffset(detId);
    HGCRechitLayer_.push_back(ilayer);   
    HGCRechitEnEtaPhiLayer_.push_back(make_tuple(Energy,eta,phi,ilayer,GPx,GPy));

   Idvar={}; 
   Idvar = std::make_pair(detId.rawId(),Energy);
    recIdEn_.push_back(Idvar);
   Idvar={}; 
   Idvar = std::make_pair(detId.rawId(),eta);
    recIdEta_.push_back(Idvar);
	 }
  

   // vector<GlobalPoint> muhittemp;
    GlobalPoint temp;
    float mueta;
    float muphi;
    float mux;
    float muy;
    float mur;
    float rechiten ;
    float rechiteta ;
    float rechitphi ;
    float rechitx ;
    float rechity ;
    float rechitr ;
    int rechitlayer ;
    float muhitdR;
    bool withineta;
  for (unsigned int l = 0; l < MuPropTrkHit_.size(); ++l) {
    temp = MuPropTrkHit_[l];
    mudetID=geom0BH->getClosestCell(temp); 
    muPropId_.push_back(mudetID);
    mueta = -999.;
    muphi = -999.;
    mueta = temp.eta();
    muphi = temp.phi();
    mux = temp.x();
    muy = temp.y();
    mur = sqrt(mux*mux+muy*muy);
    int layer = l % 14;	
    layer += 37;
    bool onlyone = true;
    withineta = TestEta(mueta,layer);
    if(withineta){
    //DenmuHGCRHEn_.push_back(rechiten);
    // std::cout<<layer<<":---------Layer"<<std::endl;
    // std::cout<<mueta<<":---------mueta"<<std::endl;
    DenmuHGCRHLayer_.push_back(layer);
    DenmuHGCRHEta_.push_back(mueta);
  for (unsigned int j = 0; j < HGCRechitEnEtaPhiLayer_.size(); ++j) {
    rechitlayer  = get<3>(HGCRechitEnEtaPhiLayer_[j]);
     //std::cout<<rechitlayer<<":rechitlayer"<<std::endl;
	if(layer==rechitlayer){
    rechiten = -999.;
    rechiteta = -999.;
    rechitphi = -999.;
    muhitdR = -999.;
    rechiten   = get<0>(HGCRechitEnEtaPhiLayer_[j]);
    rechiteta  = get<1>(HGCRechitEnEtaPhiLayer_[j]);
    rechitphi  = get<2>(HGCRechitEnEtaPhiLayer_[j]);
    rechitx  = get<4>(HGCRechitEnEtaPhiLayer_[j]);
    rechity  = get<5>(HGCRechitEnEtaPhiLayer_[j]);
    rechitr = sqrt(rechitx*rechitx+rechity*rechity);
    muhitdR = deltaR(mueta,muphi,rechiteta,rechitphi); 	
    muHGCRHEn_.push_back(rechiten);
    muHGCRHLayer_.push_back(rechitlayer);
    muHGCRHdR_.push_back(muhitdR);
    if((muhitdR<0.02) && onlyone)
	{
    p2muHGCRHEta_.push_back(mueta);
    // std::cout<<mueta<<":---------Matchedmueta"<<std::endl;
    p2muHGCRHEn_.push_back(rechiten);
    p2muHGCRHLayer_.push_back(rechitlayer);
    p2muHGCRHPhi_.push_back(rechitphi);
    p2muHGCRHX_.push_back(mux);
    p2muHGCRHY_.push_back(muy);
    p2muHGCRHR_.push_back(mur);
    //p2muHGCRHX_.push_back(rechitx);
    //p2muHGCRHY_.push_back(rechity);
    //p2muHGCRHR_.push_back(rechitr);
    mux = rechitx;
    muy = rechity;
    mur = rechitr;
	onlyone=false;
	}
	}
	}
    DenmuHGCRHX_.push_back(mux);
    DenmuHGCRHY_.push_back(muy);
    DenmuHGCRHR_.push_back(mur);
	}
	}
}


bool Phase2::TestEta( float rechiteta, int rechitlayer)
{
bool withineta = false;
switch(rechitlayer){
case 37:{ if(rechiteta>1.45 && rechiteta<1.7) withineta = true;} break;
case 38:{ if(rechiteta>1.45 && rechiteta<1.73) withineta = true;} break;
case 39:{ if(rechiteta>1.43 && rechiteta<1.75) withineta = true;} break;
case 40:{ if(rechiteta>1.40 && rechiteta<1.75) withineta = true;} break;
case 41:{ if(rechiteta>1.40 && rechiteta<1.87) withineta = true;} break;
case 42:{ if(rechiteta>1.37 && rechiteta<1.9) withineta = true;} break;
case 43:{ if(rechiteta>1.35 && rechiteta<2.05) withineta = true;} break;
case 44:{ if(rechiteta>1.35 && rechiteta<2.07) withineta = true;} break;
case 45:{ if(rechiteta>1.35 && rechiteta<2.08) withineta = true;} break;
case 46:{ if(rechiteta>1.38 && rechiteta<2.1) withineta = true;} break;
case 47:{ if(rechiteta>1.4 && rechiteta<2.25) withineta = true;} break;
case 48:{ if(rechiteta>1.4 && rechiteta<2.27) withineta = true;} break;
case 49:{ if(rechiteta>1.43 && rechiteta<2.3) withineta = true;} break;
case 50:{ if(rechiteta>1.48 && rechiteta<2.3) withineta = true;} break;

}
return withineta;
}





void Phase2::SimRecMatch( vector<pair<uint32_t,float>> simpair, vector<pair<uint32_t,float>> recpair, vector<float>& simE, vector<float>& recE, vector<int> SimLayer, vector<int>& Layer)
{
 simE.clear();
 recE.clear();
 Layer.clear();


for (unsigned int i = 0; i<simpair.size(); ++i)
{ 
uint32_t x = simpair.at(i).first;
for (unsigned int j = 0; j<recpair.size(); ++j)
{
uint32_t y = recpair.at(j).first;
if(x==y) 
{

simE.push_back(simpair.at(i).second);
recE.push_back(recpair.at(j).second);
Layer.push_back(SimLayer.at(i));

}
}


}

}
