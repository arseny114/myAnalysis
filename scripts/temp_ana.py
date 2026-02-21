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
myAnalysis.numberJets = 2  # Не используется в inclusive моде поиска джетов
myAnalysis.isolationDeltaR = 0.4
myAnalysis.outputRootFile = "{ana_path}"  # Выходной файл
myAnalysis.centerOfMassEnergy = 240.0

# Новые параметры для ee_genkt
myAnalysis.jetR = 0.5  # R для exclusive-like (эквивалент ee_kt)
myAnalysis.jetP = 1.0  # kt-like
myAnalysis.jetPtMin = 5.0  # Для inclusive
myAnalysis.useInclusive = True  # True для inclusive, False для exclusive
myAnalysis.pfoEnergyMin = 1.0 # Частицы с E < 1.0 ГэВ не будут участвовать в кластеризации и изоляции
myAnalysis.minPtForIsolation = 2.0   # Минимальный Pt для расчёта изоляции (ГэВ)
myAnalysis.isolationThreshold = 0.1  # Порог изоляции
myAnalysis.minConstPerJet = 6        # Минимальное количество частиц в джете

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis ],
    EvtSel="NONE",
    EvtMax=-1,  # Обрабатываем все события в файле
    ExtSvc=[podioevent],
)
