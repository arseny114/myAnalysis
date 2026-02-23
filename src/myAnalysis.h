// src/myAnalysis.h (updated)

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
        relativeIsolationForLeptons.clear();
        relativeIsolationForHadrons.clear();

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

        numberJetsInEvent = 0;
    }

    int eventNumber = 0;
    int numberJetsInEvent = 0;

    // Кинематика PFO
    std::vector<double> pfoE, pfoPx, pfoPy, pfoPz;
    double pfoTotalE = 0, pfoTotalPx = 0, pfoTotalPy = 0, pfoTotalPz = 0;

    // Изоляция
    std::vector<double> relativeIsolation;
    std::vector<double> relativeIsolationForLeptons;
    std::vector<double> relativeIsolationForHadrons;

    // Тип частицы (PFO) 
    std::vector<int> particleType;
    std::vector<int> isLepton;
    std::vector<int> isChargedHadron;

    // Джеты
    std::vector<std::vector<int>>   reconstructedJetConstituentsPfoIdx;
    std::vector<double> reconstructedJetPx, reconstructedJetPy, reconstructedJetPz, reconstructedJetE;
    std::vector<double> reconstructedJetThrust;

    // Количество PFO, которые входят в джет
    std::vector<int> jetSize;

    // Физические величины
    double invariantMassAllPFO = 0;
    double invariantMassJets = 0;
    double recoilMassAllPFO    = 0;
    double recoilMassJets    = 0;
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
    Gaudi::Property<int>         myNumberJets         {this, "numberJets",         2};
    Gaudi::Property<double>      myIsolationDeltaR    {this, "isolationDeltaR",    0.4};
    Gaudi::Property<std::string> myOutputFileName     {this, "outputRootFile",     "analysis_output.root"};
    Gaudi::Property<double>      myCenterOfMassEnergy {this, "centerOfMassEnergy", 240.0};

    // Новые свойства для ee_genkt_algorithm
    Gaudi::Property<double>      myJetR               {this, "jetR",               0.5};  // Радиус кластеризации R
    Gaudi::Property<double>      myJetP               {this, "jetP",               1.0};  // Параметр p (1 для kt-like)
    Gaudi::Property<double>      myJetPtMin           {this, "jetPtMin",           5.0};  // Для inclusive режима (GeV)
    Gaudi::Property<bool>        myUseInclusive       {this, "useInclusive",       true}; // Флаг: true для inclusive, false для exclusive
    Gaudi::Property<double>      myPfoEnergyMin       {this, "pfoEnergyMin",       0.5};  // Минимальная энергия PFO для джет кластеринга
    Gaudi::Property<double>      myMinPtForIsolation  {this, "minPtForIsolation",  2.0};  // Минимальный Pt для расчёта изоляции (ГэВ)
    Gaudi::Property<double>      myIsolationThreshold {this, "isolationThreshold", 0.1};  // Порог изоляции
    Gaudi::Property<size_t>      myMinConstPerJet     {this, "minConstPerJet",     6};    // Минимальное количество частиц в джете

    Gaudi::Property<bool>        myApplyJetSelection{this,   "applyJetSelection",  true}; // Включить отбор по джетам
    Gaudi::Property<bool>        myApplyIsolationSelection{this, "applyIsolationSelection", true}; // Включить отбор по изоляции

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
