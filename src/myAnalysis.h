// src/myAnalysis.h

#ifndef MY_ANALYSIS_H
#define MY_ANALYSIS_H

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"
#include "TFile.h"
#include "TTree.h"

#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include <vector>
#include <string>

/** Минимальный поперечный импульс для расчёта изоляции (ГэВ) */
constexpr double MIN_PT_FOR_ISOLATION = 2.0;

/** Порог изоляции, ниже которого частица считается изолированной */
constexpr double ISOLATION_THRESHOLD = 0.1;

/** Минимальное количество частиц в джете для принятия события */
constexpr size_t MIN_CONSTITUENTS_PER_JET = 6;

/**
 * @brief Класс для хранения данных одного события
 */
class EventData
{
public:
    EventData() { reset(); }

    void reset()
    {
        pfoE.clear();   pfoPx.clear();   pfoPy.clear();   pfoPz.clear();
        relativeIsolation.clear();

        reconstructedJetConstituentsPfoIdx.clear();
        reconstructedJetPx.clear();
        reconstructedJetPy.clear();
        reconstructedJetPz.clear();
        reconstructedJetE.clear();
        reconstructedJetThrust.clear();

        pfoTotalE = pfoTotalPx = pfoTotalPy = pfoTotalPz = 0.0;
        invariantMassAllPFO = invariantMassDijets = 0.0;
        recoilMassAllPFO = recoilMassDijets = 0.0;
    }

    int eventNumber = 0;

    // PFO
    std::vector<double> pfoE, pfoPx, pfoPy, pfoPz;
    double pfoTotalE = 0, pfoTotalPx = 0, pfoTotalPy = 0, pfoTotalPz = 0;

    // Изоляция
    std::vector<double> relativeIsolation;

    // Джеты
    std::vector<std::vector<int>>   reconstructedJetConstituentsPfoIdx;
    std::vector<double> reconstructedJetPx, reconstructedJetPy, reconstructedJetPz, reconstructedJetE;
    std::vector<double> reconstructedJetThrust;

    // Физические величины
    double invariantMassAllPFO = 0;
    double invariantMassDijets = 0;
    double recoilMassAllPFO    = 0;
    double recoilMassDijets    = 0;
};

/**
 * @brief Основной алгоритм анализа
 */
class myAnalysis : public Algorithm
{
public:
    myAnalysis(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute()    override;
    StatusCode finalize()   override;

private:
    // Настраиваемые параметры
    Gaudi::Property<std::string> myJetClusteringAlgoName {this, "jetClusteringAlgoName",   "ee_kt_algorithm"};
    Gaudi::Property<int>         myNumberJets            {this, "numberJets",              2};
    Gaudi::Property<double>      myIsolationDeltaR       {this, "isolationDeltaR",         0.4};
    Gaudi::Property<std::string> myOutputFileName        {this, "outputRootFile",          "analysis_output.root"};
    Gaudi::Property<double>      myCenterOfMassEnergy    {this, "centerOfMassEnergy",      240.0};

    // Данные текущего события
    EventData myEventData;

    // Выходной файл и дерево
    TFile* myOutputFile = nullptr;
    TTree* myOutputTree = nullptr;

    // Коллекция PFO
    DataHandle<edm4hep::ReconstructedParticleCollection> pfoCollHandler{"CyberPFOPID", Gaudi::DataHandle::Reader, this};
    const edm4hep::ReconstructedParticleCollection* myPfoCollPtr = nullptr;

    // Вспомогательные функции
    bool isInvalidPFO(const edm4hep::ReconstructedParticle& pfo) const;

    void calculateIsolationForPFO(const edm4hep::ReconstructedParticle& pfo, double deltaR);

    void saveJetClusteringResults(const std::vector<fastjet::PseudoJet>& jets);

    bool getPfoCollection();
};

#endif // MY_ANALYSIS_H
