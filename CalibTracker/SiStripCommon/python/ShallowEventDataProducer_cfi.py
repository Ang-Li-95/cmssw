import FWCore.ParameterSet.Config as cms

shallowEventRun = cms.EDProducer(
   "ShallowEventDataProducer",
   isRECO = cms.bool(True),
   trigRecord = cms.InputTag('gtDigis'),
   lumiScalers = cms.InputTag("scalersRawToDigi"),
   metadata = cms.InputTag('onlineMetaDataDigis'),
   PUInfo = cms.InputTag('addPileupInfo'),
   )
