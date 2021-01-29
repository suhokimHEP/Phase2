import FWCore.ParameterSet.Config as cms

##########################################################################################
# Setup

# this is the process run by cmsRun
process = cms.Process('LLDJ')
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
# log output
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )  ## number of events -1 does all
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# input files
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(

#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/240000/09624B53-125A-624D-AC33-91368C03150B.root',
'root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_3_0_pre1/RelValTTbar_14TeV/GEN-SIM-RECO/PU_113X_mcRun4_realistic_v1_2026D49PU200-v1/10000/00c1d9f3-2d0a-4174-9627-cc6391c58dbd.root',

 ),
)
# output name
#process.TFileService = cms.Service('TFileService', fileName = cms.string('lldjntuple_200mc_miniAOD.root'));
#process.TFileService = cms.Service('TFileService', fileName = cms.string('lldjntuple_mc_RAW.root'));
process.TFileService = cms.Service('TFileService', fileName = cms.string('Relval_RAW.root'));

## cms geometry
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

# global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '111X_mcRun4_realistic_Queue'
process.GlobalTag.globaltag = '113X_mcRun4_realistic_v3'
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
 doMiniAOD                 = cms.bool(True),
 stageL1Trigger = cms.uint32(1),
 #BadChargedCandidateFilter = cms.InputTag('BadChargedCandidateFilter'),
 #BadPFMuonFilter           = cms.InputTag('BadPFMuonFilter'),
 #pfMETLabel                = cms.InputTag('slimmedMETs'),
 triggerResults            = cms.InputTag('TriggerResults', '', 'HLT'),

 muonSrc                   = cms.InputTag('slimmedMuons'),
 patTriggerResults         = cms.InputTag('TriggerResults', '', 'PAT'),
 bits = cms.InputTag("TriggerResults","","HLT"),
 prescales = cms.InputTag("patTrigger"),
 objects = cms.InputTag("selectedPatTrigger"),
 TrackSrc               = cms.InputTag('isolatedTracks', '', 'RECO'),



 pileupCollection          = cms.InputTag('slimmedAddPileupInfo'),
 VtxLabel                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
 genParticleSrc    = cms.InputTag("genParticles"),


 RechHitSrc = cms.InputTag('g4SimHits','HGCHitsEE','SIM'),
)

#builds Ntuple
process.p = cms.Path(
    process.Phase2
    )
