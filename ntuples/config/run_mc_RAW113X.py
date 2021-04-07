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

'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/020c2f63-473b-4193-8661-42208b47b980.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/013cc9de-72bc-440a-9da2-b916446d9a76.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/009c8a3f-0eff-4846-b9f0-dd65939853ce.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/00899286-0fa4-49af-a11e-80961e6aa731.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/0552c92c-f35c-45a3-bcb9-3013e40b7dd2.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/05ba6480-df1c-4a31-acc2-f6bb905fe3d2.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/05c7a6bc-668d-4abf-930b-9bc729c3e8b2.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValZMM_14/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/08676327-5bbc-4629-9824-fad8d50261c3.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/240000/09624B53-125A-624D-AC33-91368C03150B.root',
#'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/02e69a6e-598c-4570-8aef-7978f084e438.root',
#'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/066e7900-4691-4d85-b9ad-4c65f3b0a96f.root',
#'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/10000/0ab1668c-7378-46ce-993e-3aea1486fafa.root',

 ),
)
# output name
process.TFileService = cms.Service('TFileService', fileName = cms.string('BHZMM_RAW.root'));

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




process.options.numberOfThreads=cms.untracked.uint32(8)
#NTuplizer
process.Phase2 = cms.EDAnalyzer('Phase2',
 doMiniAOD                 = cms.string('RAW'),
 HGCMode                 = cms.string('BH'),
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

)

#builds Ntuple
process.p = cms.Path(
    process.Phase2
    )
