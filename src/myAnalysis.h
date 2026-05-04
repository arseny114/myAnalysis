// src/myAnalysis.h

#ifndef MY_ANALYSIS_H
#define MY_ANALYSIS_H

#include "GaudiKernel/Algorithm.h"
#include "TFile.h"
#include "TTree.h"
#include "k4FWCore/DataHandle.h"

#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include <string>
#include <vector>

/**
 * @brief Класс для хранения данных одного события
 */
class EventData {
public:
    EventData() { reset(); }
    void reset() {
        // PFO kinematics + index
        eventNumber = 0;
        pfoE.clear();
        pfoPx.clear();
        pfoPy.clear();
        pfoPz.clear();
        pfoIdx.clear();
        leptonConeEnergy.clear();
        isIsolatedLeptonFlag.clear();
        pfoNHitsEcal.clear();
        pfoNHitsHcal.clear();
        pfoNClustersEcal.clear();
        pfoNClustersHcal.clear();
        isLepton.clear();
        isChargedHadron.clear();
        particleType.clear();

        // Inclusive jets
        nJetsInclusive = 0;
        inclusiveJetConstituentsPfoIdx.clear();
        inclusiveJetPx.clear();
        inclusiveJetPy.clear();
        inclusiveJetPz.clear();
        inclusiveJetE.clear();
        inclusiveJetThrust.clear();
        inclusiveJetSize.clear();

        // Exclusive jets
        nJetsExclusive = 0;
        exclusiveJetConstituentsPfoIdx.clear();
        exclusiveJetPx.clear();
        exclusiveJetPy.clear();
        exclusiveJetPz.clear();
        exclusiveJetE.clear();
        exclusiveJetThrust.clear();
        exclusiveJetSize.clear();
    }

    int eventNumber = 0;
    std::vector<double> pfoE, pfoPx, pfoPy, pfoPz;

    std::vector<int> pfoIdx, pfoNHitsEcal, pfoNHitsHcal, pfoNClustersEcal, pfoNClustersHcal;
    std::vector<double> leptonConeEnergy;
    std::vector<int> isIsolatedLeptonFlag, particleType, isLepton, isChargedHadron;

    // Inclusive jets
    int nJetsInclusive = 0;
    std::vector<std::vector<int>> inclusiveJetConstituentsPfoIdx;
    std::vector<double> inclusiveJetPx, inclusiveJetPy, inclusiveJetPz, inclusiveJetE,
        inclusiveJetThrust, inclusiveJetSize;

    // Exclusive jets
    int nJetsExclusive = 0;
    std::vector<std::vector<int>> exclusiveJetConstituentsPfoIdx;
    std::vector<double> exclusiveJetPx, exclusiveJetPy, exclusiveJetPz, exclusiveJetE,
        exclusiveJetThrust, exclusiveJetSize;
};

/**
 * @brief Основной алгоритм анализа
 */
class myAnalysis : public Algorithm {
public:
    myAnalysis(const std::string &name, ISvcLocator *pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:
    // Настройки по умолчанию
    //
    // =========================================================================
    // Общие настройки
    // =========================================================================
    Gaudi::Property<std::string> myOutputFileName{this, "outputRootFile", "analysis_output.root"};
    Gaudi::Property<double> myCenterOfMassEnergy{this, "centerOfMassEnergy", 240.0};

    // =========================================================================
    // Настройки изоляции лептонов (ILC-style, по энергии в конусе)
    // =========================================================================
    Gaudi::Property<double> myCosConeAngle{this, "cosConeAngle", 0.98};

    Gaudi::Property<bool> myUseRectangularIsolation{this, "useRectangularIsolation", true};
    Gaudi::Property<double> myIsoMinTrackEnergy{this, "isoMinTrackEnergy", 15.0};
    Gaudi::Property<double> myIsoMaxTrackEnergy{this, "isoMaxTrackEnergy", 1e20};
    Gaudi::Property<double> myIsoMinConeEnergy{this, "isoMinConeEnergy", 0.0};
    Gaudi::Property<double> myIsoMaxConeEnergy{this, "isoMaxConeEnergy", 2.0};

    Gaudi::Property<bool> myUsePolynomialIsolation{this, "usePolynomialIsolation", false};
    Gaudi::Property<double> myIsoPolynomialA{this, "isoPolynomialA", 0.0};
    Gaudi::Property<double> myIsoPolynomialB{this, "isoPolynomialB", 20.0};
    Gaudi::Property<double> myIsoPolynomialC{this, "isoPolynomialC", -300.0};

    // =========================================================================
    // Настройки кластеризации джетов
    // =========================================================================
    Gaudi::Property<int> myNumberJets{this, "numberJets", 2};
    Gaudi::Property<double> myJetR{this, "jetR", 0.5};
    Gaudi::Property<double> myJetP{this, "jetP", 1.0};
    Gaudi::Property<double> myJetPtMin{this, "jetPtMin", 5.0};

    // =========================================================================
    // Параметры геометрии детектора (для определения ECAL/HCAL)
    // =========================================================================
    Gaudi::Property<double> ecalRMax{this, "ECALRMax", 2130.0,
                                     "Максимальный радиус барреля ECAL [мм]"};
    Gaudi::Property<double> ecalZMax{this, "ECALZMax", 3230.0, "Максимальная |Z| торца ECAL [мм]"};

    // =========================================================================
    // Настройки сбора статистики хитов
    // =========================================================================
    Gaudi::Property<bool> myCollectHitStats{
        this, "collectHitStats", true,
        "Собирать статистику хитов в ECAL/HCAL (true) или пропускать (false)"};

    // Данные текущего события
    EventData myEventData;

    // Выходной файл и дерево
    TFile *myOutputFile = nullptr;
    TTree *myOutputTree = nullptr;

    // Коллекция PFO
    DataHandle<edm4hep::ReconstructedParticleCollection> pfoCollHandler{
        "CyberPFOPID", Gaudi::DataHandle::Reader, this};
    const edm4hep::ReconstructedParticleCollection *myPfoCollPtr = nullptr;

    // =========================================================================
    // Вспомогательные структуры
    // =========================================================================
    /**
     * @brief Структура для хранения статистики хитов в детекторах
     */
    struct ClusterHitStats {
        int nHitsEcal = 0;     // Число хитов в ECAL
        int nHitsHcal = 0;     // Число хитов в HCAL
        int nClustersEcal = 0; // Число кластеров в ECAL
        int nClustersHcal = 0; // Число кластеров в HCAL

        int totalHits() const { return nHitsEcal + nHitsHcal; }
        int totalClusters() const { return nClustersEcal + nClustersHcal; }
    };

    // =========================================================================
    // Вспомогательные функции
    // =========================================================================
    bool isInvalidPFO(const edm4hep::ReconstructedParticle &pfo) const;
    bool pfoIsLepton(const edm4hep::ReconstructedParticle &pfo) const;
    bool pfoIsChargedHadron(const edm4hep::ReconstructedParticle &pfo) const;

    // Вспомогательные функции для изоляции лептонов (ILC-style)
    double getConeEnergy(const edm4hep::ReconstructedParticle &pfo,
                         const edm4hep::ReconstructedParticleCollection *allPfos) const;
    bool isIsolatedRectangular(double trackEnergy, double coneEnergy) const;
    bool isIsolatedPolynomial(double trackEnergy, double coneEnergy) const;
    bool isIsolatedLepton(const edm4hep::ReconstructedParticle &pfo,
                          const edm4hep::ReconstructedParticleCollection *allPfos) const;

    void fillJetData(const std::vector<fastjet::PseudoJet> &jets, bool isInclusive);

    bool getPfoCollection();

    /**
     * @brief Подсчитывает хиты в кластерах данного PFO с разделением на ECAL/HCAL
     * @param pfo Реконструированная частица
     * @return Структура со статистикой хитов и кластеров
     */
    ClusterHitStats countClusterHits(const edm4hep::ReconstructedParticle &pfo) const;

    /**
     * @brief Определяет, находится ли кластер в ECAL по его позиции
     * @param cluster Ссылка на кластер
     * @return true если кластер в ECAL, false если в HCAL
     */
    bool isClusterInEcal(const edm4hep::Cluster &cluster) const;
};

#endif // MY_ANALYSIS_H
