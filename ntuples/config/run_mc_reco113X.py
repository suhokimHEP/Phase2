import FWCore.ParameterSet.Config as cms

##########################################################################################
# Setup

# this is the process run by cmsRun
process = cms.Process('Phase2')
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
# log output
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )  ## number of events -1 does all
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# input files
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(

'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-RECO/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/17a96b62-df6b-44aa-b3db-f27b89d401fe.root',
#'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/0061d0f8-9134-46dc-a78b-07d776bb3d64.root',

 ),
)
# output name
#process.TFileService = cms.Service('TFileService', fileName = cms.string('Relval_Reco.root'));
process.TFileService = cms.Service('TFileService', fileName = cms.string('FHZMMRelval_Reco.root'));

# global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '111X_mcRun4_realistic_Queue'
process.GlobalTag.globaltag = '113X_mcRun4_realistic_v3'


## cms geometry
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
##process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
####process.load('Configuration.StandardSequences.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
##process.load("Configuration.Geometry.GeometryECALHCAL_cff")



########################  BadPFMuonFilter
#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag('slimmedMuons', '', 'PAT')
#process.BadPFMuonFilter.PFCandidates = cms.InputTag('packedPFCandidates', '', 'PAT')
#
########################  BadChargedCandidateFilter             
#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag('slimmedMuons', '', 'PAT')
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag('packedPFCandidates', '', 'PAT')
#
#process.lldjMETFiltersSequence = cms.Sequence(
#     process.BadPFMuonFilter *
#     process.BadChargedCandidateFilter 
#)

###########################################################################################
## For AOD Track variables
## 
#process.MaterialPropagator = cms.ESProducer('PropagatorWithMaterialESProducer',
#    ComponentName = cms.string('PropagatorWithMaterial'),
#    Mass = cms.double(0.105),
#    MaxDPhi = cms.double(1.6),
#    PropagationDirection = cms.string('alongMomentum'),
#    SimpleMagneticField = cms.string(''),
#    ptMin = cms.double(-1.0),
#    useRungeKutta = cms.bool(False)
#)
#
#process.TransientTrackBuilderESProducer = cms.ESProducer('TransientTrackBuilderESProducer',
#    ComponentName = cms.string('TransientTrackBuilder')
#)

#NTuplizer
process.Phase2 = cms.EDAnalyzer('Phase2',
 doMiniAOD                 = cms.string('reco'),
 stageL1Trigger = cms.uint32(1),

 triggerResults            = cms.InputTag('TriggerResults', '', 'HLT'),
 patTriggerResults         = cms.InputTag('TriggerResults', '', 'PAT'),
 bits = cms.InputTag("TriggerResults","","HLT"),
 prescales = cms.InputTag("patTrigger"),
 objects = cms.InputTag("selectedPatTrigger"),

 slimmuonSrc                   = cms.InputTag('slimmedMuons'),
 muonSrc                   = cms.InputTag('muons','','RECO'),

 TrackSrc               = cms.InputTag('generalTracks', '', 'RECO'),
 isoTrackSrc               = cms.InputTag('isolatedTracks', '', 'RECO'),

 pileupCollection          = cms.InputTag('addPileupInfo'),
 slimpileupCollection          = cms.InputTag('slimmedAddPileupInfo'),

 slimVtxLabel                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
 VtxLabel                  = cms.InputTag('offlinePrimaryVertices'),

 genParticleSrc    = cms.InputTag("genParticles"),

 EESimHitSrc = cms.InputTag('g4SimHits','HGCHitsEE','SIM'),
 FHSimHitSrc = cms.InputTag('g4SimHits','HGCHitsHEfront','SIM'),
 BHSimHitSrc = cms.InputTag('g4SimHits','HGCHitsHEback','SIM'),
 EERecHits = cms.InputTag('HGCalRecHit','HGCEERecHits','RECO'),
 FHRecHits = cms.InputTag('HGCalRecHit','HGCHEFRecHits','RECO'),
 BHRecHits = cms.InputTag('HGCalRecHit','HGCHEBRecHits','RECO'),
 #EERecHits = cms.InputTag('HGCalRecHit:HGCEERecHits'),
 #FHRecHits = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
 #BHRecHits = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),

)

#builds Ntuple
process.p = cms.Path(
    process.Phase2
    )
