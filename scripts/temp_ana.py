# Файл шаблон для генерации py скриптов

import os, sys, glob
from Gaudi.Configuration import *

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="{rec_path}")

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
myAnalysis.jetClusteringAlgoName = "ee_kt_algorithm"
myAnalysis.numberJets = 2
myAnalysis.isolationDeltaR = 0.4
myAnalysis.outputRootFile = "{ana_path}"  # Выходной файл
myAnalysis.centerOfMassEnergy = 240.0

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis ],
    EvtSel="NONE",
    EvtMax=-1,  # Обрабатываем все события в файле
    ExtSvc=[podioevent],
)
