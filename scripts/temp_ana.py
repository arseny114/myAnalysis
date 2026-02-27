# Файл шаблон для генерации py скриптов
import os, sys, glob
from Gaudi.Configuration import *

########### k4DataSvc ####################
from Configurables import k4DataSvc
podioevent = k4DataSvc("EventDataSvc", input="{rec_path}")
##########################################

########## CEPCSWData #################
cepcswdatatop = "/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/data/latest"
#######################################

########## Podio Input ###################
from Configurables import PodioInput
inp = PodioInput("InputReader")
inp.collections = ["CyberPFOPID"]
##########################################

from Configurables import myAnalysis
myAnalysis = myAnalysis("myAnalysis")

# =============================================================================
# Общие настройки
# =============================================================================
myAnalysis.outputRootFile = "{ana_path}"           # Выходной файл
myAnalysis.centerOfMassEnergy = 240.0              # Полная энергия в системе центра масс (ГэВ)
myAnalysis.pfoEnergyMin = 1.0                      # Минимальная энергия PFO для кластеризации джетов и расчета изоляции (ГэВ)

# =============================================================================
# Настройки изоляции частиц
# =============================================================================
myAnalysis.isolationDeltaR = 0.4                   # Радиус конуса для расчёта изоляции (ΔR)
myAnalysis.minPtForIsolation = 1.0                 # Минимальный Pt частицы для расчёта изоляции (ГэВ)
myAnalysis.isolationThreshold = 2.0                # Порог изоляции (если ниже, то частица считается изолированной)

# =============================================================================
# Настройки кластеризации джетов
# =============================================================================
myAnalysis.numberJets = 2                          # Желаемое количество джетов (для exclusive режима)
myAnalysis.useInclusive = True                     # True для inclusive, False для exclusive
myAnalysis.jetR = 0.5                              # Радиус джета R для ee_genkt_algorithm
myAnalysis.jetP = 1.0                              # Параметр p для ee_genkt_algorithm (1 для kt-like)
myAnalysis.jetPtMin = 5.0                          # Минимальный pT джета для inclusive режима (ГэВ)
myAnalysis.minConstPerJet = 6                      # Минимальное количество частиц в джете

# =============================================================================
# Настройки отбора событий
# =============================================================================
myAnalysis.applyJetSelection = False               # Включить отбор по джетам
myAnalysis.applyIsolationSelection = True          # Включить отбор по изоляции

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis],
    EvtSel="NONE",
    EvtMax=-1,  # Обрабатываем все события в файле
    ExtSvc=[podioevent],
)
