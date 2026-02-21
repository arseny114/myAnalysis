import os, sys, glob
from Gaudi.Configuration import *

########### k4DataSvc ####################/cefs/higgs/zhangkl/Production/2501/E240_nnHbb_0117/Combined
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="/cefs/higgs/liugeliang/CEPC/202501/Production/Hinvi/E240_e1e1Hinvi/Combined/rec_E240_e1e1Hinvi_00001.root")

##########################################


########## CEPCSWData ################# 
cepcswdatatop ="/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/data/latest"
#######################################


########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = [ "CyberPFOPID" ]
##########################################


from Configurables import myAnalysis
myAnalysis = myAnalysis("myAnalysis")
myAnalysis.numberJets = 2
myAnalysis.isolationDeltaR = 0.4
myAnalysis.outputRootFile = "test_output.root"
myAnalysis.centerOfMassEnergy = 240.0

# Новые параметры для ee_genkt
myAnalysis.jetR = 0.5  # Большой R для exclusive-like (эквивалент ee_kt)
myAnalysis.jetP = 1.0  # kt-like
myAnalysis.jetPtMin = 5.0  # Для inclusive
myAnalysis.useInclusive = False  # False для exclusive, True для inclusive
myAnalysis.pfoEnergyMin = 1.0 # Частицы с E < 1.0 ГэВ не будут участвовать в кластеризации и изоляции
myAnalysis.minPtForIsolation = 2.0   # Минимальный Pt для расчёта изоляции (ГэВ)
myAnalysis.isolationThreshold = 0.1  # Порог изоляции
myAnalysis.minConstPerJet = 6        # Минимальное количество частиц в джете

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis ],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[podioevent],
)
