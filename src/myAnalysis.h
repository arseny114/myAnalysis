#ifndef MY_ANALYSIS
#define MY_ANALYSIS

/* TODO вот это я все от балды пока что скопировал, надо потом разораться с инклудами */
#include "UTIL/ILDConf.h"
#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"
#include <random>
#include <vector>
#include <string>
#include "GaudiKernel/NTuple.h"
#include "TFile.h"
#include "TTree.h"

#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticleCollectionData.h"
#include "edm4hep/RecDqdx.h"
#include "edm4hep/RecDqdxCollection.h"
#include "edm4hep/ParticleIDCollection.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

/*
 * Вот эта штука это какое-то расширение edm4hep от команды CEPC
 * ее репозиторий - еhttps://code.ihep.ac.cn/cepc/externals/EDM4cepc
 */
#include "edm4cepc/RecTof.h"
#include "edm4cepc/RecTofCollection.h"

/* Дефайны для настроек кода */

/* 
 * Минимальный Pt при котором для этого 
 * PFO будет вычисляться изоляция 
 */
#define MIN_PT_MOMENTUM_FOR_ISOLATION     2         /* GeV */

#define MAX_DELTA_R_FOR_ISOLATION         0.40f     /* Radians */

/* 
 * Количество возможных гипотез для 
 * которых мы считаем хи-квадраты
 */
#define NUMBER_OF_PARTICLE_HYPOTHESES     5

/* 
 * Создаем вспомогательный класс для
 * хранения данных одного события 
 */
class EventData 
{
public:
  /* Конструктор, инициализирует все переменные */
  EventData() { reset(); };

  /* Функция для сборса всех переменных */
  void reset()
  {
    relativeIsolationCharged.clear();
    relativeIsolationNeutral.clear();

    pfoChi2TOF.clear();
    pfoChi2TPC.clear();

    pfoE.clear();
    pfoPx.clear();
    pfoPy.clear();
    pfoPz.clear();

    pfoTotalE = 0;
    pfoTotalPx = 0;
    pfoTotalPy = 0;
    pfoTotalPz = 0;

    reconstructedJetConstituentsPfoIdx.clear();
    reconstructedJetPx.clear();
    reconstructedJetPy.clear();
    reconstructedJetPz.clear();
    reconstructedJetE.clear();
    reconstructedJetThrust.clear();

    mcParticleE.clear();
    mcParticlePx.clear();
    mcParticlePy.clear();
    mcParticlePz.clear();
    mcParticlePDGid.clear();
    mcParticleCharge.clear();
    mcParticleGeneratorStatus.clear();
    mcParticleParentsSize.clear();
    mcJetPx.clear();
    mcJetPy.clear();
    mcJetPz.clear();
    mcJetE.clear();
    mcJetThrust.clear();

    mcTotalE = 0;
    mcTotalNeutrinosE = 0;
    mcTotalPx = 0;
    mcTotalPy = 0;
    mcTotalPz = 0;

    minDeltaR.clear();
    matchedMCIdx.clear();
  }

  /* Номер обрабатываемого события */
  int     eventNumber;

  /* 
   * Величина изоляции для каждого pfo объекта из этого события.
   *
   * Значения в первом векторе это изоляция (относительная) 
   * относительно заряженных соседей. Значения во втором векторе 
   * это изоляция (относительная) относительно нейтральных соседей.
   */
  std::vector<double> relativeIsolationCharged;
  std::vector<double> relativeIsolationNeutral;

  /* Значения хи квадратов PID, которые вычесленны для каждого PFO */
  std::vector<std::vector<double>> pfoChi2TOF;
  std::vector<std::vector<double>> pfoChi2TPC;

  /* Переменные для PFO коллекции  */
  std::vector<double> pfoE;
  std::vector<double> pfoPx;
  std::vector<double> pfoPy;
  std::vector<double> pfoPz;

  double pfoTotalE;
  double pfoTotalPx;
  double pfoTotalPy;
  double pfoTotalPz;

  /* Составляющие джетов, которые нашел джет кластеринг */
  std::vector<std::vector<int>> reconstructedJetConstituentsPfoIdx;
  std::vector<double> reconstructedJetPx;
  std::vector<double> reconstructedJetPy;
  std::vector<double> reconstructedJetPz;
  std::vector<double> reconstructedJetE;
  std::vector<double> reconstructedJetThrust;

  /* Переменные для MC коллекции */
  std::vector<double> mcParticleE;
  std::vector<double> mcParticlePx;
  std::vector<double> mcParticlePy;
  std::vector<double> mcParticlePz;
  std::vector<int> mcParticlePDGid;
  std::vector<int> mcParticleCharge;
  std::vector<int> mcParticleGeneratorStatus;
  std::vector<int> mcParticleParentsSize;
  std::vector<double> mcJetPx;
  std::vector<double> mcJetPy;
  std::vector<double> mcJetPz;
  std::vector<double> mcJetE;
  std::vector<double> mcJetThrust;

  double mcTotalE;
  double mcTotalNeutrinosE;
  double mcTotalPx;
  double mcTotalPy;
  double mcTotalPz;

  /* Соотвествие MC и PFO */
  std::vector<double> minDeltaR;
  std::vector<int> matchedMCIdx;
};

/* Gaudi::Property - настройки конфигурации (он как раз берет настройки из Python) */

/* 
 * NOTE: Как я понял для джет кластеринга используется библиотека FastJet:
 * github - https://github.com/scikit-hep/fastjet/
 * doc - https://fastjet.readthedocs.io/en/latest/ (питоновский интерфейс)
 * doc - http://fastjet.fr (с++ интерфейс, но сайт работает почему-то только с vpn)
 */

/*
 * NOTE: Можно еще почитать документацию по Gaudi на их сайте
 */

/*
 * NOTE: для того чтобы понять какие методы есть у тайпдефов выше нужно
 * почтать документацию edm4hep (видимо key4hep это тоже самое):
 * что-то похожее на doc - https://key4hep.github.io/key4hep-doc/main/how-tos/key4hep-tutorials/edm4hep_analysis/edm4hep_api_intro.html 
 * Еще есть вот это - https://edm4hep.web.cern.ch/classes.html
 */

/* 
 * Добавляем свой алгоритм анализа в Gaudi, чтобы он его вызвал
 * (как я понимаю тут происходит наследование от Algorithm, 
 * который является базовым классом Gaudi)
 */
class myAnalysis : public Algorithm
{
public:

  /* 
   * Конструктуор класса:
   * name - имя алгоритма
   * pSvcLocator - указатель на навигатор сервисов (он позволяет находить другие компоненты системы)
   */
  myAnalysis( const std::string& name, ISvcLocator* pSvcLocator );

  /* 
   * override это проверка компилятором того,
   * что такая функция есть у родительского класса.
   *
   * Как я понимаю по такому шаблону объявяются
   * все Gaudi алгоритмы.
   */
  StatusCode initialize() override;   /* Инициализация алгоритма, вызывается 1 раз */
  StatusCode execute() override;      /* Выполняется для каждого события (ее вызывает ApplicationMgr) */
  StatusCode finalize() override;     /* Выполняется 1 раз в конце, для отчистки */

private:
  /*
   * NOTE: m_PFOColHdl и тд это просто имена переменных в моем коде
   * NOTE: this это указатель на текущий объект (чтобы Gaudi знал кому принадлежит handle)
   */
  DataHandle<edm4hep::ReconstructedParticleCollection> pfoCollHandler{"CyberPFOPID", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::MCParticleCollection> mcCollHandler{"MCParticle", Gaudi::DataHandle::Reader, this};

  /* 
   * Это используется для идентефикации частиц: 
   * TOF - time if flight - изменение времени пролета частицы для определения массы
   * DQDX - dQ / dx - измерение потерь энергии в детекторе для идентефикации частиц
   */
  DataHandle<edm4hep::RecTofCollection> tofCollHandler{"RecTofCollection", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::RecDqdxCollection> dqdxCollHandler{"DndxTracks", Gaudi::DataHandle::Reader, this};

  /* 
   * Gaudi::Property это настройки из конфигурации (видимо из питоновского скрипта) 
   *
   * Они задаются вот в таком виде:
   *
   * Gaudi::Property<Тип> имя_переменной{ 
   * this, 
   * "Имя_В_Конфиге", 
   * Значение_По_Умолчанию 
   * };
   *
   * Вроде бы этим штукам не нужно declareProperty(), потому что они сами себя регают (тк они уже Property)
   */
  Gaudi::Property<std::string> 
    myJetClusteringAlgoName{this, "jetClusteringAlgoName", "ee_kt_algorithm"};
  Gaudi::Property<int> myNumberJets{this, "numberJets", 2};
  Gaudi::Property<double> myJetsR{this, "jetsR", 0.6};
  Gaudi::Property<std::string> 
    myOutputFileName{this, "outputRootFile", "myAnalysisOutput.root"};

  /* Все переменные обрабатываемого события хранятся тут */
  EventData myEventData;

  /* Выходной файл и выходное дерево */
  TFile   *myOutputFile;
  TTree   *myOutputTree;

  /* Обязательные коллекции */
  const edm4hep::ReconstructedParticleCollection *myPfoCollPtr;
  const edm4hep::MCParticleCollection *myMcPartCollPtr;

  /* Необязательные коллекции */
  const edm4hep::RecTofCollection     *myTofCollPtr;
  const edm4hep::RecDqdxCollection    *myDqdxCollPtr;

  /* Функция для проверки инвалидности элемента PFO */
  bool isInvalidPFO(const edm4hep::ReconstructedParticle& pfo) const;

  /* 
   * Функция для идентефикации частиц по данным из коллекции dqdx
   */
  void doDqdxParticleIdentification(const edm4hep::ReconstructedParticle& pfoElement);

  /* 
   * Функция для идентефикации частиц по коллекциям dqdx и tof 
   */
  void doTOFParticleIdentification(const edm4hep::ReconstructedParticle& pfoElement);

  /* 
   * Функция для вычисления изоляции PFO.
   */
  void calculateIsolationForPFO(const edm4hep::ReconstructedParticle& pfoElement, 
                                double deltaR);

  /* 
   * Функция для обработки результатов джет кластеринга.
   *
   * Она предназначена как для генераторных джетов, так и для реконструированных.
   */
  void processJetClusteringResults(const std::vector<fastjet::PseudoJet>& jets, 
                                   bool isReconstructed);

  /* 
   * Функция для сопоставления реконструкции с MC.
   * inputParticles - должен быть созан на исходных PFO.
   */
  void doGenMatch(const std::vector<fastjet::PseudoJet>& inputParticles,
                  const std::vector<fastjet::PseudoJet>& inputParticlesMC);

  /* Для сброса переменных класса */
  void reset() { myEventData.reset(); }

  /* Для получения коллекций */
  bool getPfoCollection();
  bool getMcCollection();
  bool getTofCollection();
  bool getDqdxCollection();
};

#endif                  /* MY_ANALYSIS */
