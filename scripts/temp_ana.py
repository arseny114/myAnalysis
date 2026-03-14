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
myAnalysis.pfoEnergyMin = 1.0                      # Минимальная энергия PFO для кластеризации джетов (ГэВ)

# =============================================================================
# Настройки изоляции лептонов (ILC-style, по энергии в конусе)
# =============================================================================

# Косинус полу-угла конуса для расчёта изоляции.
# cos(θ) = 0.98 соответствует углу θ ≈ 11.5° вокруг направления лептона.
# Частицы считаются "внутри конуса", если cos(Δθ) >= cosConeAngle,
# где Δθ вычисляется через скалярное произведение импульсов.
myAnalysis.cosConeAngle = 0.98                     # ~11.5° полу-угол конуса

# Включить прямоугольные (прямоугольные в пространстве [E_track, E_cone]) критерии изоляции.
# Если True, применяются пороги isoMin/MaxTrackEnergy и isoMin/MaxConeEnergy.
myAnalysis.useRectangularIsolation = True

# Минимальная энергия трека-кандидата в изолированные лептоны (ГэВ).
# Лептоны с энергией ниже этого порога не проходят отбор как "изолированные".
myAnalysis.isoMinTrackEnergy = 15.0                # ГэВ

# Максимальная энергия трека для требования изоляции.
# Обычно оставляют очень большим (1e20), чтобы не отрезать высокоэнергетические лептоны.
myAnalysis.isoMaxTrackEnergy = 1e20

# Минимальная допустимая энергия в конусе изоляции.
# Обычно 0.0, т.к. нас интересует верхний порог ("тишина" вокруг лептона).
myAnalysis.isoMinConeEnergy = 0.0

# Максимальная суммарная энергия других частиц в конусе вокруг лептона (ГэВ).
# Если энергия в конусе превышает это значение, лептон считается НЕизолированным.
# Типичные значения: 1–5 ГэВ; 2.0 ГэВ это хороший баланс между эффективностью и чистотой.
myAnalysis.isoMaxConeEnergy = 2.0                  # ГэВ. Порог "тишины" вокруг лептона

# Включить полиномиальные критерии изоляции вида:
#   E_cone² < A·E_track² + B·E_track + C
# Позволяет задавать энергозависимую границу изоляции (гибче прямоугольных порогов).
# Можно использовать совместно с rectangular (логическое И) или вместо них.
myAnalysis.usePolynomialIsolation = False

# Коэффициенты полиномиального критерия изоляции:
#   (E_cone)² < A·(E_track)² + B·(E_track) + C
#
# isoPolynomialA = 0.0: нет квадратичной зависимости от энергии трека
# isoPolynomialB = 20.0: линейный член, задаёт наклон границы
# isoPolynomialC = -300.0: сдвиг границы; отрицательное значение ужесточает отбор при малых E_track
#
# Пример: при E_track = 20 ГэВ порог на E_cone = sqrt(0*400 + 20*20 - 300) = sqrt(100) = 10 ГэВ
# При таких параметрах полиномиальный критерий мягче прямоугольного (2 ГэВ) для высокоэнергетических лептонов.
myAnalysis.isoPolynomialA = 0.0
myAnalysis.isoPolynomialB = 20.0
myAnalysis.isoPolynomialC = -300.0

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
