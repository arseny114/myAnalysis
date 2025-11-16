/*
 * Описание того что это
 */

#ifndef MY_ANALYSIS
#define MY_ANALYSIS

/* TODO вот это я все от балды пока что скопировал, надо потом разораться с инклудами */
#include "UTIL/ILDConf.h"
#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"
#include <random>
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
 * Попробую тут добавить тайпдефы чтобы 
 * просто
 * было легче читать C++ код
 * TODO выглядит сомнительно я бы удалил это
 * (потому что тут тайпдефы вообще на хендлеры, а не на коллекции)
 */
typedef DataHandle<edm4hep::ReconstructedParticleCollection> myRecPartColl;
typedef DataHandle<edm4hep::MCParticleCollection> myMCPartColl;
typedef DataHandle<edm4hep::RecTofCollection> myRecTofColl;
typedef DataHandle<edm4hep::RecDqdxCollection> myRecDqdxColl;

/* TODO это видимо просто значит, что у нас режим на чтение */
typedef Gaudi::DataHandle::Reader myGaudiReader;

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

  /* TODO надо понять какие мне нужны 
   *
   * NOTE: m_PFOColHdl и тд это просто имена переменных в моем коде
   * NOTE: this это указатель на текущий объект (чтобы Gaudi знал кому принадлежит handle)
   *
   * Я пока что возьму те же самые коллекции
   *
   * Тут надо сразу написать какие есть методы у этих хендлеров:
   * .get()
   *
   *
   *
   *
   */
  myRecPartColl pfoCollHandler{"CyberPFOPID", myGaudiReader, this};
  myMCPartColl mcCollHandler{"MCParticle", myGaudiReader, this};

  /* 
   * Это используется для идентефикации частиц: 
   * TOF - time if flight - изменение времени пролета частицы для определения массы
   * DQDX - dQ / dx - измерение потерь энергии в детекторе для идентефикации частиц
   */
  myRecTofColl tofCollHandler{"RecTofCollection", myGaudiReader, this};
  myRecDqdxColl dqdxCollHandler{"DndxTracks", myGaudiReader, this};

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
    myOutputRootFile{this, "outputRootFile", "myAnalysisOutput.root"}

  /* 
   * TODO
   * В MissingET есть еще некоторые функции ниже, я пока что их добавлять не буду 
   * потому что не оч пока понял зачем именно их в этом классе писать, почему не написать
   * их вне этого класса внутри myAnalysis.cpp 
   */

  /* Для извлечения коллекции (шаблонная функция должна быть реализована тут) */
  template<typename CollectionType>
  bool getCollection(const DataHandle<CollectionType>& handle, 
                     const CollectionType **collection, 
                     const int eventNumber);
  {
    try 
    {
      *collection = handle.get();
    }
    catch (const GaudiException& e)
    {
      info() << "Collection " << handle.fullKey() << 
            " is unavailable in event " << eventNumber << endmsg;
      return false;
    }

    if ((*collection)->empty())
    {
      info() << "Collection " << handle.fullKey() << 
            " is empty in event " << eventNumber << endmsg;
      return false;
    }

    return true;
  }

  /* Функция для проверки инвалидности элемента PFO */
  bool isInvalidPFO(const edm4hep::ReconstructedParticle& pfo) const 
  {
    return (std::isnan(pfo.getMomentum()[0]) || 
            std::isnan(pfo.getMomentum()[1]) || 
            std::isnan(pfo.getMomentum()[2]) || 
            pfo.getEnergy() == 0);
  }

  /* 
   * Функция для идентефикации частиц по данным из коллекции dqdx
   */
  void doDqdxParticleIdentification(const edm4hep::RecDqdxCollection *dqdxColl, 
                                    edm4hep::MutableReconstructedParticle &pfoElement)
  {
    /* 
     * Тут будут храниться значения chi^2 для каждой из гипотез:
     * (изнанчально заполнены значениями -1)
     * 0 - электрон
     * 1 - мюон
     * 2 - пион
     * 3 - каон
     * 4 - протон
     *
     * Я не уверен, что я правильно определил номера гипотез, но
     * тут можно найти некий конфиг файл в котором есть эти номера:
     *
     * https://code.ihep.ac.cn/cepc/externals/EDM4cepc/-/blob/main/
     * edm4cepc.yaml?ref_type=heads
     *
     * Помимо номеров гипотез в нем есть (как я понимаю) вообще все
     * члены TOF коллекции (которая добавлена именно разработчиками CEPCSW).
     */
    std::vector<double> chi2Tpc(NUMBER_OF_PARTICLE_HYPOTHESES, -1);

    /* 
     * Запускаем цикл по всем трекам для конкретного PFO.
     *
     * По какой-то причине в PFO может быть несколько
     * треков, это может быть вызвано либо ошибкой реконструкции,
     * либо же еще какими-то эффектами реконструкции (типа вторичных частиц наверное).
     */
    // TODO а точно ли у PFO может быть больше одного трека?
    // я просто не понимаю как мы возврашаем результаты этой функции тогда
    for (auto trk : pfoElement.getTracks()) /* Вроде возвращает 1 трек */
    {
      /* 
       * Цикл по всем элементам коллекции dqdx.
       *
       * TPC измеряет dE/dx (потери энергии на единицу пути), разные
       * частицы имеют разные характерные потери при заданном импульсе.
       * Для каждой гипотезы вычисляется хи-квадрат исходя из ожидаемых потерь.
       * (насколько хорошо полученные потери согласуются с ожидаемыми).
       *
       * Методы этой коллекции и в целом edm4hep можно поискать в nauch_maga/EDM4hep
       * (я скачал ее с сервера CEPC, надеюсь это нужная версия)
       *
       * Каждый элемент коллекции dqdx содержит связь с конкретным треком.
       */
      for (auto dqdxElement : *dqdxColl) 
      {
        if (dqdxElement.getTrack() == trk) /* Тоже возвращает один трек, но 
                                            * тут другое название метода почему-то */
          for (int i = 0; i < NUMBER_OF_PARTICLE_HYPOTHESES; i++)
            chi2Tpc[i] = dqdxElement.getHypotheses(i).chi2;
        
        break; /* Завершаем цикл если нашли нужный трек */
      }
    }

    /* 
     * Запысываем резуьтат для последнего трека PFO 
     * (FIXME нужно чтобы не только для последнего было) 
     */
    pfoChi2TPC.push_back(chi2Tpc);
  }


  /* 
   * Функция для идентефикации частиц по коллекциям dqdx и tof 
   */
  void doTOFParticleIdentification(const edm4hep::RecTofCollection *tofColl, 
                                   edm4hep::MutableReconstructedParticle &pfoElement)
  {
    /* 
     * Тут будут храниться значения chi^2 для каждой из гипотез:
     * (изнанчально заполнены значениями -1)
     * 0 - электрон
     * 1 - мюон
     * 2 - пион
     * 3 - каон
     * 4 - протон
     *
     * Я не уверен, что я правильно определил номера гипотез, но
     * тут можно найти некий конфиг файл в котором есть эти номера:
     *
     * https://code.ihep.ac.cn/cepc/externals/EDM4cepc/-/blob/main/
     * edm4cepc.yaml?ref_type=heads
     *
     * Помимо номеров гипотез в нем есть (как я понимаю) вообще все
     * члены TOF коллекции (которая добавлена именно разработчиками CEPCSW).
     */
    std::vector<double> chi2Tof(NUMBER_OF_PARTICLE_HYPOTHESES, -1);

    /* 
     * Запускаем цикл по всем трекам для конкретного PFO.
     *
     * По какой-то причине в PFO может быть несколько
     * треков, это может быть вызвано либо ошибкой реконструкции,
     * либо же еще какими-то эффектами реконструкции (типа вторичных частиц наверное).
     */
    // TODO а точно ли у PFO может быть больше одного трека?
    // я просто не понимаю как мы возврашаем результаты этой функции тогда
    for (auto trk : pfoElement.getTracks()) /* Вроде возвращает 1 трек */
    {
      /* 
       * Цикл по всем элементам коллекции tof (reconstructed time of flight). 
       *
       * Делает примерно тоже самое, что и цикл выше, но на основе измерений
       * времени пролета частиц.
       *
       * Хи-квадрат тут считаем вручную по формуле:
       * chi2 = [(t_ожидаемое - t_измеренное) / погрешность]^2
       */
      for (auto tofElement : *tofColl)
      {
        /* TODO тут странная формула для хи квадрат была, я поменял знаменатель */
        if (tofElement.getTrack() == trk)
          chi2Tof[i] = std::pow((tofElement.getTime() - tofElement.getTimeExp()[i]) / 
                                tofElement.getTimeExp()[i], 2);

        break; /* Завершаем цикл если нашли нужный трек */
      }
    }

    /* 
     * Запысываем резуьтат для последнего трека PFO 
     * (FIXME нужно чтобы не только для последнего было) 
     */
    pfoChi2TOF.push_back(chi2Tof);
  }

  /* 
   * Функция для вычисления изоляции PFO.
   */
  void calculateIsolationForPFO(edm4hep::ReconstructedParticleCollection *myPfoCollPtr,
                                edm4hep::MutableReconstructedParticle &pfoCollElement, 
                                EventData &eventData,
                                double deltaR)
  {
    /* Создаем 4-вектор для переданного PFO объекта */
    TLorentzVector pfo4Momentum;
    TLorentzVector pfo4MomentumNeighbour;

    /* Суммарный импульс */
    double totalMomentumCharge = 0; 
    double totalMomentumNeutral = 0;

    /* Заполняем 4 вектор для обрабатываемого pfo объекта */
    pfo4Momentum.SetPxPyPzE(pfoCollElement.getMomentum()[0], 
                            pfoCollElement.getMomentum()[1], 
                            pfoCollElement.getMomentum()[2], 
                            pfoCollElement.getEnergy());

    /* 
     * Для частиц у которых поперечный импульс
     * меньше 2 ГэВ мы считаем изоляцию равной -1.
     * (по сути просто не считаем)
     */
    if (pfo4Momentum.Pt() < MIN_PT_MOMENTUM_FOR_ISOLATION)
    {
      eventData.relativeIsolationCharged.push_back(-1);
      eventData.relativeIsolationNeutral.push_back(-1);
      return;
    }

    /* Все проверки пройдены, тогда вычисляем изоляцию */

    /* 
     * Запускаем цикл по всем pfo элеметам 
     * коллекции чтобы посчитать изоляцию
     */
    for (const auto& pfoCollElementNeighbour : *myPfoCollPtr)
    {
      /* 
       * Пропускаем невалидные + пропускаем 
       * совпадение с нашим обрабатываемым.
       */
      if (isInvalidPFO(pfoCollElementNeighbour) || 
          &pfoCollElementNeighbour == &pfoCollElement)
        continue;

      /* Также заполняем 4 вектор для текущего соседа */
      pfo4MomentumNeighbour.SetPxPyPzE(pfoCollElementNeighbour.getMomentum()[0], 
                                       pfoCollElementNeighbour.getMomentum()[1], 
                                       pfoCollElementNeighbour.getMomentum()[2], 
                                       pfoCollElementNeighbour.getEnergy());    

      /* 
       * Если эта соседняя частица входит 
       * в интересующий нас конус 
       */
      if (pfo4Momentum.DeltaR(pfo4MomentumNeighbour) < deltaR)
      {
        /* Если заряд соседа */
        if (pfoCollElementNeighbour.getCharge() == 0)
          totalMomentumCharge += pfoCollElementNeighbour.P();
        else
          totalMomentumNeutral += pfoCollElementNeighbour.P();
      }
    }

    /* 
     * Теперь нормируем изоляцию на импульс самого pfo 
     * и запысываем результат в eventData 
     */
    eventData.relativeIsolationCharged.push_back(totalMomentumCharge / pfo4Momentum.P());
    eventData.relativeIsolationNeutral.push_back(totalMomentumNeutral / pfo4Momentum.P());
  }

  /* 
   * Функция для обработки результатов джет кластеринга.
   *
   * Она предназначена как для генераторных джетов, так и для реконструированных.
   */
  void processJetClusteringResults(const std::vector<fastjet::PseudoJet>& jets, 
                                   bool isReconstructed,
                                   EventData& eventData)
  {
    /* Запускаем цикл по полученным джетам чтобы извлечь их параметры */
    for (int i = 0; i < eventData.jets.size(); i++) 
    {
      /* 
       * TODO сохранение всех переменных для джета типа энергии и тд,
       * пока что я не буду это сохранять так как пока не принял
       * решение какие перменные мы будем сохранят в EventData
       */
      
    }

  }

  /* 
   * TODO нужно добавить сюда все, что я хочу считать
   * Также нам необходимо поместить сюда все нужные нам переменные, 
   * которые должны быть доступны между вызовами функций
   */

  /* Номер обрабатываемого события */
  int     eventNumber;

  /* Суммарное количество обработанных частиц */
  int numberParticles;

  /* Все переменные обрабатываемого события хранятся тут */
  EventData myEventData;

  /* Выходной файл и выходное дерево */
  TFile   *outputFile;
  TTree   *outputTree;

  /* 
   * Значения хи квадратов PID, которые вычесленны для каждого PFO.
   * TODO В текущей реализации мы берем только последний трек для PFO, 
   * это нужно исправить
   *
   * TODO это мб надо перенести в EventData вместе с самими коллекциями
   */
  std::vector<std::vector<double>> pfoChi2TPC, pfoChi2TOF;
};

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
  /* TODO  сюда надо добавить все используемые перменные */
  void reset()
  {

  }

  /* TODO все перменные как публичные члены класса */

  /* 
   * Величина изоляции для каждого pfo объекта из этого события.
   *
   * Значения в первом векторе это изоляция (относительная) 
   * относительно заряженных соседей. Значения во втором векторе 
   * это изоляция (относительная) относительно нейтральных соседей.
   */
  std::vector<double> relativeIsolationCharged;
  std::vector<double> relativeIsolationNeutral;

  /* Суммарный 4 импульс всех pfo в событии */
  TLorentzVector totalPFO4Momentum;

  /* Джеты, которые нашел алгоритм джет кластеринга */
  vector<fastjet::PseudoJet> jets;

};


























#endif                  /* MY_ANALYSIS */
