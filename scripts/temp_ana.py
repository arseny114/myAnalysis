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

# =============================================================================
# Настройки изоляции лептонов (ILC-style, по энергии в конусе)
# =============================================================================
# Косинус полу-угла конуса для расчёта изоляции
myAnalysis.cosConeAngle = 0.985                    
# Включить прямоугольные критерии изоляции
myAnalysis.useRectangularIsolation = True          
# Минимальная энергия трека-кандидата (ГэВ)
myAnalysis.isoMinTrackEnergy = 5.0                 
# Максимальная энергия трека для требования изоляции
myAnalysis.isoMaxTrackEnergy = 1e20                
# Минимальная допустимая энергия в конусе изоляции
myAnalysis.isoMinConeEnergy = 0.0                  
# Максимальная суммарная энергия других частиц в конусе вокруг лептона (ГэВ)
myAnalysis.isoMaxConeEnergy = 2.0                  
# Включить полиномиальные критерии изоляции
myAnalysis.usePolynomialIsolation = False          
# Коэффициенты полинома: E_cone² < A·E_track² + B·E_track + C
myAnalysis.isoPolynomialA = 0.0                    
myAnalysis.isoPolynomialB = 20.0                   
myAnalysis.isoPolynomialC = -300.0                 

# =============================================================================
# Настройки кластеризации джетов
# =============================================================================
myAnalysis.numberJets = 2                          # Желаемое количество джетов (для exclusive режима)
myAnalysis.jetR = 0.5                              # Радиус джета R для ee_genkt_algorithm (для inclusive режима)
myAnalysis.jetP = 1.0                              # Параметр p для ee_genkt_algorithm (1 для kt-like)
myAnalysis.jetPtMin = 5.0                          # Минимальный pT джета для inclusive режима (ГэВ)

# =============================================================================
# Параметры геометрии детектора (для определения ECAL/HCAL)
# =============================================================================
myAnalysis.ECALRMax = 2130.0                       # Максимальный радиус барреля ECAL [мм]
myAnalysis.ECALZMax = 3230.0                       # Максимальная |Z| торца ECAL [мм]

# =============================================================================
# Настройки сбора статистики хитов
# =============================================================================
myAnalysis.collectHitStats = False                 # Собирать статистику хитов (True/False)

########################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg=[inp, myAnalysis],
    EvtSel="NONE",
    EvtMax=-1,  # Обрабатываем все события в файле
    ExtSvc=[podioevent],
)
