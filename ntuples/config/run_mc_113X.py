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

'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/100000/00D62BEB-4E52-E243-980F-D473E59EACDE.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/100000/067F7D95-C15A-B241-A2EA-B0A9AF1C512F.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/NoPU_111X_mcRun4_realistic_T15_v1-v1/100000/53807AB8-9399-8942-A573-499BA283F18C.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_111X_mcRun4_realistic_T15_v1-v1/260000/052F3259-5977-3248-88DF-1C1D50753CE6.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_111X_mcRun4_realistic_T15_v1-v1/260000/016F48F9-F967-8E46-B1DF-5B382572DB60.root',
#'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_111X_mcRun4_realistic_T15_v1-v1/260000/06256FC2-7EA3-C245-A04B-DE607D985F5D.root',

 ),
)
# output name
#process.TFileService = cms.Service('TFileService', fileName = cms.string('lldjntuple_200mc_miniAOD.root'));
process.TFileService = cms.Service('TFileService', fileName = cms.string('lldjntuple_mc_miniAOD.root'));

# cms geometry
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

# global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '111X_mcRun4_realistic_Queue'
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
 doMiniAOD                 = cms.bool(False),
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
