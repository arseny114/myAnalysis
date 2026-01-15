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
inp.collections = [ "CyberPFOPID", "MCParticle" , "RecTofCollection", "DndxTracks"]
##########################################


from Configurables import myAnalysis
myAnalysis = myAnalysis("myAnalysis")
myAnalysis.jetClusteringAlgoName = "ee_kt_algorithm"
myAnalysis.numberJets = 2
myAnalysis.jetsR = 0.6
myAnalysis.outputRootFile = "test_output.root"

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis ],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[podioevent],
)
