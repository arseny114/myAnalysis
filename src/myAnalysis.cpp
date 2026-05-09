// src/myAnalysis.cpp

/**
 * @file myAnalysis.cpp
 * @brief Реализация алгоритма myAnalysis для обработки событий в CEPCSW
 *
 * Алгоритм выполняет следующие основные задачи:
 * - чтение реконструированных частиц (PFO)
 * - расчёт изоляции частиц
 * - кластеризация на джеты с помощью fastjet
 */

#include "myAnalysis.h"

// Стандартные библиотеки
#include <chrono>
#include <cmath>

// Компоненты Gaudi
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// Интерфейсы детектора и единицы измерения
#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DetInterface/IGeomSvc.h"

// FastJet — кластеризация джетов
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

// ROOT — четырёхвекторы и сохранение в файл
#include "TLorentzVector.h"

// Макрос Gaudi для регистрации компонента
DECLARE_COMPONENT(myAnalysis)

/**
 * @brief Конструктор алгоритма
 *
 * Здесь регистрируются все настраиваемые параметры, которые можно задать
 * из python-конфигурационного файла.
 */
myAnalysis::myAnalysis(const std::string &name, ISvcLocator *pSvcLocator)
    : Algorithm(name, pSvcLocator) {
    // =========================================================================
    // Общие настройки
    // =========================================================================
    declareProperty("outputRootFile", myOutputFileName, "Имя выходного ROOT-файла");
    declareProperty("centerOfMassEnergy", myCenterOfMassEnergy,
                    "Полная энергия в системе центра масс (ГэВ)");

    // =========================================================================
    // Настройки изоляции лептонов (ILC-style, по энергии в конусе)
    // =========================================================================
    declareProperty("cosConeAngle", myCosConeAngle,
                    "Косинус полу-угла конуса для расчёта изоляции (cosθ >= cosConeAngle)");
    declareProperty("useRectangularIsolation", myUseRectangularIsolation,
                    "Использовать прямоугольные критерии изоляции");
    declareProperty("isoMinTrackEnergy", myIsoMinTrackEnergy,
                    "Минимальная энергия трека для требования изоляции (ГэВ)");
    declareProperty("isoMaxTrackEnergy", myIsoMaxTrackEnergy,
                    "Максимальная энергия трека для требования изоляции (ГэВ)");
    declareProperty("isoMinConeEnergy", myIsoMinConeEnergy,
                    "Минимальная энергия в конусе для требования изоляции (ГэВ)");
    declareProperty("isoMaxConeEnergy", myIsoMaxConeEnergy,
                    "Максимальная энергия в конусе для требования изоляции (ГэВ)");
    declareProperty("usePolynomialIsolation", myUsePolynomialIsolation,
                    "Использовать полиномиальные критерии изоляции");
    declareProperty("isoPolynomialA", myIsoPolynomialA,
                    "Коэффициент A полиномиального критерия: E_cone² < A·E² + B·E + C");
    declareProperty("isoPolynomialB", myIsoPolynomialB,
                    "Коэффициент B полиномиального критерия: E_cone² < A·E² + B·E + C");
    declareProperty("isoPolynomialC", myIsoPolynomialC,
                    "Коэффициент C полиномиального критерия: E_cone² < A·E² + B·E + C");

    // =========================================================================
    // Настройки кластеризации джетов
    // =========================================================================
    declareProperty("numberJets", myNumberJets, "Желаемое количество джетов");
    declareProperty("jetR", myJetR, "Радиус джета R для ee_genkt_algorithm");
    declareProperty("jetP", myJetP, "Параметр p для ee_genkt_algorithm (1 для kt-like)");
    declareProperty("jetPtMin", myJetPtMin, "Минимальный pT джета для inclusive режима (GeV)");

    // =========================================================================
    // Параметры геометрии детектора (для определения ECAL/HCAL)
    // =========================================================================
    declareProperty("ECALRMax", ecalRMax, "Максимальный радиус барреля ECAL [мм]");
    declareProperty("ECALZMax", ecalZMax, "Максимальная |Z| торца ECAL [мм]");

    // =========================================================================
    // Настройки сбора статистики хитов
    // =========================================================================
    declareProperty("collectHitStats", myCollectHitStats,
                    "Собирать статистику хитов в ECAL/HCAL (true) или пропускать (false)");
}

/**
 * @brief Инициализация алгоритма (вызывается один раз в начале)
 *
 * Создаёт структуру данных события, открывает выходной файл,
 * создаёт дерево и все необходимые ветви.
 */
StatusCode myAnalysis::initialize() {
    myEventData = EventData();
    myEventData.eventNumber = 0;

    // Открываем файл в режиме пересоздания
    myOutputFile = TFile::Open(myOutputFileName.value().c_str(), "RECREATE");
    if (!myOutputFile || myOutputFile->IsZombie()) {
        error() << "Не удалось создать выходной файл: " << myOutputFileName << endmsg;
        return StatusCode::FAILURE;
    }

    myOutputTree = new TTree("outputTree", "Results");

    // PFO kinematics + index
    myOutputTree->Branch("pfoIdx", &myEventData.pfoIdx);
    myOutputTree->Branch("pfoE", &myEventData.pfoE);
    myOutputTree->Branch("pfoPx", &myEventData.pfoPx);
    myOutputTree->Branch("pfoPy", &myEventData.pfoPy);
    myOutputTree->Branch("pfoPz", &myEventData.pfoPz);
    myOutputTree->Branch("pfoCharge", &myEventData.pfoCharge);

    // Hit stats
    myOutputTree->Branch("pfoNHitsEcal", &myEventData.pfoNHitsEcal);
    myOutputTree->Branch("pfoNHitsHcal", &myEventData.pfoNHitsHcal);
    myOutputTree->Branch("pfoNClustersEcal", &myEventData.pfoNClustersEcal);
    myOutputTree->Branch("pfoNClustersHcal", &myEventData.pfoNClustersHcal);

    // Isolation & particle types
    myOutputTree->Branch("leptonConeEnergy", &myEventData.leptonConeEnergy);
    myOutputTree->Branch("isIsolatedLeptonFlag", &myEventData.isIsolatedLeptonFlag);
    myOutputTree->Branch("particleType", &myEventData.particleType);
    myOutputTree->Branch("isLepton", &myEventData.isLepton);
    myOutputTree->Branch("isChargedHadron", &myEventData.isChargedHadron);

    // Inclusive jets
    myOutputTree->Branch("nJetsInclusive", &myEventData.nJetsInclusive);
    myOutputTree->Branch("inclusiveJetConstituentsPfoIdx",
                         &myEventData.inclusiveJetConstituentsPfoIdx);
    myOutputTree->Branch("inclusiveJetPx", &myEventData.inclusiveJetPx);
    myOutputTree->Branch("inclusiveJetPy", &myEventData.inclusiveJetPy);
    myOutputTree->Branch("inclusiveJetPz", &myEventData.inclusiveJetPz);
    myOutputTree->Branch("inclusiveJetE", &myEventData.inclusiveJetE);
    myOutputTree->Branch("inclusiveJetThrust", &myEventData.inclusiveJetThrust);
    myOutputTree->Branch("inclusiveJetSize", &myEventData.inclusiveJetSize);

    // Exclusive jets
    myOutputTree->Branch("nJetsExclusive", &myEventData.nJetsExclusive);
    myOutputTree->Branch("exclusiveJetConstituentsPfoIdx",
                         &myEventData.exclusiveJetConstituentsPfoIdx);
    myOutputTree->Branch("exclusiveJetPx", &myEventData.exclusiveJetPx);
    myOutputTree->Branch("exclusiveJetPy", &myEventData.exclusiveJetPy);
    myOutputTree->Branch("exclusiveJetPz", &myEventData.exclusiveJetPz);
    myOutputTree->Branch("exclusiveJetE", &myEventData.exclusiveJetE);
    myOutputTree->Branch("exclusiveJetThrust", &myEventData.exclusiveJetThrust);
    myOutputTree->Branch("exclusiveJetSize", &myEventData.exclusiveJetSize);

    // Event counter
    myOutputTree->Branch("eventNumber", &myEventData.eventNumber);

    return StatusCode::SUCCESS;
}

/**
 * @brief Основная функция обработки каждого события
 */
StatusCode myAnalysis::execute() {
    TLorentzVector totalPFO4Momentum(0., 0., 0., 0.);
    int pfoIndex = 0;
    std::vector<fastjet::PseudoJet> inputParticles;

    // Получаем коллекцию PFO
    if (!getPfoCollection()) {
        return StatusCode::SUCCESS;
    }

    // Сбрасываем данные события
    myEventData.reset();

    // ── 1. Проходим по всем PFO в событии ───────────────────────────────
    for (const auto &pfo : *myPfoCollPtr) {
        // Пропускаем невалидные
        if (isInvalidPFO(pfo))
            continue;

        // Увеличиваем индекс при обработке нового PFO
        pfoIndex++;
        myEventData.pfoIdx.push_back(pfoIndex);

        // Сохраняем кинематику
        myEventData.pfoE.push_back(pfo.getEnergy());
        myEventData.pfoPx.push_back(pfo.getMomentum()[0]);
        myEventData.pfoPy.push_back(pfo.getMomentum()[1]);
        myEventData.pfoPz.push_back(pfo.getMomentum()[2]);
        myEventData.pfoCharge.push_back(pfo.getCharge());

        // Записываем тип частицы
        myEventData.particleType.push_back(pfo.getType());
        myEventData.isLepton.push_back(pfoIsLepton(pfo) ? 1 : 0);
        myEventData.isChargedHadron.push_back(pfoIsChargedHadron(pfo) ? 1 : 0);

        // Подсчёт хитов в кластерах для этого PFO (только если включено)
        if (myCollectHitStats.value()) {
            // Подсчёт хитов в кластерах для этого PFO
            auto hitStats = countClusterHits(pfo);

            // Сохраняем в векторы (по одному значению на PFO)
            myEventData.pfoNHitsEcal.push_back(hitStats.nHitsEcal);
            myEventData.pfoNHitsHcal.push_back(hitStats.nHitsHcal);
            myEventData.pfoNClustersEcal.push_back(hitStats.nClustersEcal);
            myEventData.pfoNClustersHcal.push_back(hitStats.nClustersHcal);
        } else {
            // Заполняем заглушками, чтобы размеры векторов совпадали
            myEventData.pfoNHitsEcal.push_back(-1);
            myEventData.pfoNHitsHcal.push_back(-1);
            myEventData.pfoNClustersEcal.push_back(-1);
            myEventData.pfoNClustersHcal.push_back(-1);
        }

        // Суммируем полный четырёхимпульс события
        totalPFO4Momentum += TLorentzVector(pfo.getMomentum()[0], pfo.getMomentum()[1],
                                            pfo.getMomentum()[2], pfo.getEnergy());

        // Готовим частицу для кластеризации
        // emplace_back() это эффективный способ добавления в конец вектора
        // (он также создает объект PseudoJet с переданными параметрами импульса)
        inputParticles.emplace_back(pfo.getMomentum()[0], pfo.getMomentum()[1],
                                    pfo.getMomentum()[2], pfo.getEnergy());
        // Сохраняем для восстановления состава джета
        inputParticles.back().set_user_index(pfoIndex);

        // Считаем изоляцию только для лептонов
        if (pfoIsLepton(pfo)) {
            double coneEnergy = getConeEnergy(pfo, myPfoCollPtr);
            myEventData.leptonConeEnergy.push_back(coneEnergy);

            bool isolated = isIsolatedLepton(pfo, myPfoCollPtr);
            myEventData.isIsolatedLeptonFlag.push_back(isolated ? 1 : 0);
        } else {
            // Для не-лептонов или мягких частиц записываем заглушки
            myEventData.leptonConeEnergy.push_back(-1.0);
            myEventData.isIsolatedLeptonFlag.push_back(0);
        }
    }

    // ── 2. Кластеризация джетов ───────────────────────────────────
    fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, myJetR.value(), myJetP.value());
    fastjet::ClusterSequence cs(inputParticles, jet_def);

    auto incJets = fastjet::sorted_by_pt(cs.inclusive_jets(myJetPtMin.value()));
    auto excJets = fastjet::sorted_by_pt(cs.exclusive_jets(myNumberJets.value()));

    // Заполняем данные для каждого режима через обычную функцию-член
    fillJetData(incJets, true);  // inclusive-режим
    fillJetData(excJets, false); // exclusive-режим

    // Заполняем дерево и увеличиваем счётчик событий
    myEventData.eventNumber++;
    myOutputTree->Fill();
    return StatusCode::SUCCESS;
}

/**
 * @brief Завершение работы алгоритма (вызывается один раз в конце)
 */
StatusCode myAnalysis::finalize() {
    if (myOutputFile) {
        myOutputFile->Write();
        myOutputFile->Close();
    }
    return StatusCode::SUCCESS;
}

/**
 * @brief Проверка, является ли PFO невалидным
 */
bool myAnalysis::isInvalidPFO(const edm4hep::ReconstructedParticle &pfo) const {
    const auto &mom = pfo.getMomentum();
    return (std::isnan(mom[0]) || std::isnan(mom[1]) || std::isnan(mom[2]) || pfo.getEnergy() <= 0);
}

/**
 * @brief Проверка, является ли PFO лептоном.
 */
bool myAnalysis::pfoIsLepton(const edm4hep::ReconstructedParticle &pfo) const {
    int absPfoType = std::abs(pfo.getType());
    return (absPfoType == 11 || absPfoType == 13);
}

/**
 * @brief Проверка, является ли PFO заряженным адроном.
 */
bool myAnalysis::pfoIsChargedHadron(const edm4hep::ReconstructedParticle &pfo) const {
    int absPfoType = std::abs(pfo.getType());
    return (absPfoType == 2212 || absPfoType == 321 || absPfoType == 211);
}

/**
 * @brief Вычисляет энергию в конусе вокруг данной PFO-частицы
 *
 * Конус определяется через косинус угла между импульсами:
 * cosTheta = (P₁·P₂) / (|P₁|·|P₂|) >= cosConeAngle
 *
 * @param pfo Ссылка на объект реконструированной частицы
 * @param allPfos Коллекция всех PFO в событии
 * @return Суммарная энергия частиц в конусе (без учёта самой частицы)
 */
double myAnalysis::getConeEnergy(const edm4hep::ReconstructedParticle &pfo,
                                 const edm4hep::ReconstructedParticleCollection *allPfos) const {
    double coneEnergy = 0.0;
    TVector3 P_main(pfo.getMomentum()[0], pfo.getMomentum()[1], pfo.getMomentum()[2]);
    double pMain = P_main.Mag();

    if (pMain < 1e-9)
        return 0.0;

    for (const auto &other : *allPfos) {
        if (other.getObjectID() == pfo.getObjectID())
            continue;
        if (isInvalidPFO(other))
            continue;

        TVector3 P_other(other.getMomentum()[0], other.getMomentum()[1], other.getMomentum()[2]);
        double pOther = P_other.Mag();
        if (pOther < 1e-9)
            continue;

        double cosTheta = P_main.Dot(P_other) / (pMain * pOther);

        if (cosTheta >= myCosConeAngle.value()) {
            coneEnergy += other.getEnergy();
        }
    }
    return coneEnergy;
}

/**
 * @brief Проверка прямоугольных критериев изоляции
 */
bool myAnalysis::isIsolatedRectangular(double trackEnergy, double coneEnergy) const {
    if (trackEnergy < myIsoMinTrackEnergy.value())
        return false;
    if (trackEnergy > myIsoMaxTrackEnergy.value())
        return false;
    if (coneEnergy < myIsoMinConeEnergy.value())
        return false;
    if (coneEnergy > myIsoMaxConeEnergy.value())
        return false;
    return true;
}

/**
 * @brief Проверка полиномиальных критериев изоляции
 *
 * Критерий: E_cone² < A·E_track² + B·E_track + C
 */
bool myAnalysis::isIsolatedPolynomial(double trackEnergy, double coneEnergy) const {
    double coneE2 = coneEnergy * coneEnergy;
    double threshold = myIsoPolynomialA.value() * trackEnergy * trackEnergy +
                       myIsoPolynomialB.value() * trackEnergy + myIsoPolynomialC.value();
    return coneE2 < threshold;
}

/**
 * @brief Проверка, является ли лептон изолированным по критериям ILC
 */
bool myAnalysis::isIsolatedLepton(const edm4hep::ReconstructedParticle &pfo,
                                  const edm4hep::ReconstructedParticleCollection *allPfos) const {
    if (!pfoIsLepton(pfo))
        return false;

    // Если оба критерия выключены, то изоляция не применяется
    if (!myUseRectangularIsolation.value() && !myUsePolynomialIsolation.value())
        return false;

    double trackEnergy = pfo.getEnergy();
    double coneEnergy = getConeEnergy(pfo, allPfos);

    bool passedRectangular =
        !myUseRectangularIsolation.value() || isIsolatedRectangular(trackEnergy, coneEnergy);
    bool passedPolynomial =
        !myUsePolynomialIsolation.value() || isIsolatedPolynomial(trackEnergy, coneEnergy);

    return passedRectangular && passedPolynomial;
}

/**
 * @brief Заполняет данные джетов для одного из режимов кластеризации
 */
void myAnalysis::fillJetData(const std::vector<fastjet::PseudoJet> &jets, bool isInclusive) {
    // Выбираем, в какие векторы писать данные, в зависимости от режима
    std::vector<double> *jetPx =
        isInclusive ? &myEventData.inclusiveJetPx : &myEventData.exclusiveJetPx;
    std::vector<double> *jetPy =
        isInclusive ? &myEventData.inclusiveJetPy : &myEventData.exclusiveJetPy;
    std::vector<double> *jetPz =
        isInclusive ? &myEventData.inclusiveJetPz : &myEventData.exclusiveJetPz;
    std::vector<double> *jetE =
        isInclusive ? &myEventData.inclusiveJetE : &myEventData.exclusiveJetE;
    std::vector<double> *jetTh =
        isInclusive ? &myEventData.inclusiveJetThrust : &myEventData.exclusiveJetThrust;
    std::vector<double> *jetSz =
        isInclusive ? &myEventData.inclusiveJetSize : &myEventData.exclusiveJetSize;
    std::vector<std::vector<int>> *jetIdx = isInclusive
                                                ? &myEventData.inclusiveJetConstituentsPfoIdx
                                                : &myEventData.exclusiveJetConstituentsPfoIdx;
    int *nJets = isInclusive ? &myEventData.nJetsInclusive : &myEventData.nJetsExclusive;

    // Очищаем целевые векторы перед заполнением
    jetPx->clear();
    jetPy->clear();
    jetPz->clear();
    jetE->clear();
    jetTh->clear();
    jetSz->clear();
    jetIdx->clear();
    *nJets = static_cast<int>(jets.size());

    // Проходим по всем джетам
    for (const auto &jet : jets) {
        int sz = static_cast<int>(jet.constituents().size());
        jetSz->push_back(sz);

        // Сохраняем 4-импульс джета
        jetPx->push_back(jet.px());
        jetPy->push_back(jet.py());
        jetPz->push_back(jet.pz());
        jetE->push_back(jet.E());

        // Расчёт thrust и индексов конституентов
        double thrust = 0.0, pSum = 0.0;
        TLorentzVector j4(jet.px(), jet.py(), jet.pz(), jet.E());
        std::vector<int> cIdx;

        for (const auto &c : jet.constituents()) {
            cIdx.push_back(c.user_index());
            TLorentzVector c4(c.px(), c.py(), c.pz(), c.E());
            thrust += std::abs(c4.Vect().Dot(j4.Vect())) / j4.P();
            pSum += c4.P();
        }
        if (pSum > 1e-9)
            thrust /= pSum;

        jetTh->push_back(thrust);
        jetIdx->push_back(std::move(cIdx));
    }
}

/**
 * @brief Получение коллекции PFO из хранилища событий
 */
bool myAnalysis::getPfoCollection() {
    try {
        myPfoCollPtr = pfoCollHandler.get();
    } catch (const GaudiException &e) {
        info() << "Коллекция PFO недоступна в событии " << myEventData.eventNumber << endmsg;
        return false;
    }

    if (myPfoCollPtr->empty()) {
        info() << "Коллекция PFO пуста в событии " << myEventData.eventNumber << endmsg;
        return false;
    }

    return true;
}

/**
 * @brief Определяет, находится ли кластер в ECAL по геометрическим порогам
 */
bool myAnalysis::isClusterInEcal(const edm4hep::Cluster &cluster) const {
    // Получаем позицию кластера (Vector3f)
    // Это энергетически взвешенный центр тяжести всех хитов кластера в мм
    auto pos = cluster.getPosition();
    // Радиальное расстояние от оси пучка
    // hypot(x, y) это r = sqrt(x^2 + y^2)
    double r = std::hypot(pos.x, pos.y);
    // Абсолютное значение по оси Z
    double z = std::fabs(pos.z);

    // ECAL: внутри цилиндра радиусом ecalRMax и высотой ±ecalZMax
    return (r < ecalRMax.value() && z < ecalZMax.value());
}

/**
 * @brief Подсчитывает хиты в кластерах данного PFO с разделением на ECAL/HCAL
 */
myAnalysis::ClusterHitStats
myAnalysis::countClusterHits(const edm4hep::ReconstructedParticle &pfo) const {

    // Создаем экземпляр структуры для подсчета хитов
    ClusterHitStats stats;

    // Перебираем все кластеры, ассоциированные с этим PFO
    for (const auto &cluster : pfo.getClusters()) {
        // Количество хитов в текущем кластере
        int nHits = static_cast<int>(cluster.hits_size());

        // Определяем детектор по позиции кластера
        if (isClusterInEcal(cluster)) {
            stats.nHitsEcal += nHits;
            stats.nClustersEcal++;
        } else {
            stats.nHitsHcal += nHits;
            stats.nClustersHcal++;
        }
    }

    return stats;
}
