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
        pfoE.clear();
        pfoPx.clear();
        pfoPy.clear();
        pfoPz.clear();
        leptonConeEnergy.clear();
        isIsolatedLeptonFlag.clear();

        reconstructedJetConstituentsPfoIdx.clear();
        reconstructedJetPx.clear();
        reconstructedJetPy.clear();
        reconstructedJetPz.clear();
        reconstructedJetE.clear();
        reconstructedJetThrust.clear();

        isLepton.clear();
        isChargedHadron.clear();
        particleType.clear();

        jetSize.clear();

        pfoTotalE = pfoTotalPx = pfoTotalPy = pfoTotalPz = 0.0;
        invariantMassAllPFO = invariantMassJets = 0.0;
        recoilMassAllPFO = recoilMassJets = 0.0;

        skippedByJets = skippedByIsolatedLepton = 0;

        numberJetsInEvent = 0;
    }

    int eventNumber = 0;
    int numberJetsInEvent = 0;

    // Кинематика PFO
    std::vector<double> pfoE, pfoPx, pfoPy, pfoPz;
    double pfoTotalE = 0, pfoTotalPx = 0, pfoTotalPy = 0, pfoTotalPz = 0;

    // Изоляция лептонов (ветки для отладки)
    std::vector<double> leptonConeEnergy;  // Энергия в конусе для лептонов
    std::vector<int> isIsolatedLeptonFlag; // Флаг: 1 = лептон изолирован

    // Было ли пропущено событие из-за изоляции
    int skippedByIsolatedLepton = 0;

    // Тип частицы (PFO)
    std::vector<int> particleType;
    std::vector<int> isLepton;
    std::vector<int> isChargedHadron;

    // Джеты
    std::vector<std::vector<int>> reconstructedJetConstituentsPfoIdx;
    std::vector<double> reconstructedJetPx, reconstructedJetPy, reconstructedJetPz,
        reconstructedJetE;
    std::vector<double> reconstructedJetThrust;

    // Было ли пропущено событие из-за условий на джеты
    int skippedByJets = 0;

    // Количество PFO, которые входят в джет
    std::vector<int> jetSize;

    // Физические величины
    double invariantMassAllPFO = 0;
    double invariantMassJets = 0;
    double recoilMassAllPFO = 0;
    double recoilMassJets = 0;
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
    Gaudi::Property<double> myPfoEnergyMin{this, "pfoEnergyMin", 0.5};

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
    Gaudi::Property<bool> myUseInclusive{this, "useInclusive", true};
    Gaudi::Property<double> myJetR{this, "jetR", 0.5};
    Gaudi::Property<double> myJetP{this, "jetP", 1.0};
    Gaudi::Property<double> myJetPtMin{this, "jetPtMin", 5.0};
    Gaudi::Property<size_t> myMinConstPerJet{this, "minConstPerJet", 6};

    // =========================================================================
    // Настройки отбора событий
    // =========================================================================
    Gaudi::Property<bool> myApplyJetSelection{this, "applyJetSelection", true};
    Gaudi::Property<bool> myApplyIsolationSelection{this, "applyIsolationSelection", true};

    // Данные текущего события
    EventData myEventData;

    // Выходной файл и дерево
    TFile *myOutputFile = nullptr;
    TTree *myOutputTree = nullptr;

    // Коллекция PFO
    DataHandle<edm4hep::ReconstructedParticleCollection> pfoCollHandler{
        "CyberPFOPID", Gaudi::DataHandle::Reader, this};
    const edm4hep::ReconstructedParticleCollection *myPfoCollPtr = nullptr;

    // Вспомогательные функции
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

    void saveJetClusteringResults(const std::vector<fastjet::PseudoJet> &jets);

    bool getPfoCollection();
};

#endif // MY_ANALYSIS_H
