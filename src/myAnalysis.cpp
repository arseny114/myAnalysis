#include "myAnalysis.h"

/* TODO вот это я все от балды пока что скопировал, надо потом разораться с инклудами */
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <math.h>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
// #include "fastjet/contrib/Nsubjettiness.hh" TODO уже было в комменте

#include <chrono>

#include "TLorentzVector.h"

/*
 * Какой-то макрос из Gaudi, который нужен для создания
 * фабрики, которая будет создавать объекты моего класса
 */
DECLARE_COMPONENT(myAnalysis)

/* Список инициализации для класса */
myAnalysis::myAnalysis(const std::string& name, ISvcLocator* pSvcLocator ) :
  Algorithm( name, pSvcLocator ) 
{ 
  /* Объявлем входные параметры из питоновского файла */
  declareProperty("jetClusteringAlgoName", myJetClusteringAlgoName, "jet Clustering Algo Name");
  declareProperty("numberJets", myNumberJets, "Number Jets");
  declareProperty("jetsR", myJetsR, "Jets R");
  declareProperty("outputRootFile", myOutputFileName, "output Root File");
}

/*
 * ================================================
 *
 * Функция инициализации
 *
 *        Вызывается в самом начале один раз
 *
 * ================================================
 */
StatusCode myAnalysis::initialize()
{
  /* Начинаем обработку с нулевого события */
  myEventData = EventData();
  myEventData.eventNumber = 0;

  /* Создаем и открываем выходной файл + дерево */
  myOutputFile = TFile::Open(myOutputFileName.value().c_str(), "RECREATE");
  myOutputTree = new TTree("outputTree", "outputTree");

  /* Создаем все ветки для заполнения данными из анализа */
  myOutputTree->Branch("pfoE", &myEventData.pfoE);
  myOutputTree->Branch("pfoPx", &myEventData.pfoPx);
  myOutputTree->Branch("pfoPy", &myEventData.pfoPy);
  myOutputTree->Branch("pfoPz", &myEventData.pfoPz);

  myOutputTree->Branch("pfoTotalE", &myEventData.pfoTotalE);
  myOutputTree->Branch("pfoTotalPx", &myEventData.pfoTotalPx);
  myOutputTree->Branch("pfoTotalPy", &myEventData.pfoTotalPy);
  myOutputTree->Branch("pfoTotalPz", &myEventData.pfoTotalPz);

  myOutputTree->Branch("relativeIsolationCharged", &myEventData.relativeIsolationCharged);
  myOutputTree->Branch("relativeIsolationNeutral", &myEventData.relativeIsolationNeutral);

  myOutputTree->Branch("pfoChi2TOF", &myEventData.pfoChi2TOF);
  myOutputTree->Branch("pfoChi2TPC", &myEventData.pfoChi2TPC);

  myOutputTree->Branch("reconstructedJetConstituentsPfoIdx", &myEventData.reconstructedJetConstituentsPfoIdx);
  myOutputTree->Branch("reconstructedJetPx", &myEventData.reconstructedJetPx);
  myOutputTree->Branch("reconstructedJetPy", &myEventData.reconstructedJetPy);
  myOutputTree->Branch("reconstructedJetPz", &myEventData.reconstructedJetPz);
  myOutputTree->Branch("reconstructedJetE", &myEventData.reconstructedJetE);
  myOutputTree->Branch("reconstructedJetThrust", &myEventData.reconstructedJetThrust);

  myOutputTree->Branch("mcParticleE", &myEventData.mcParticleE);
  myOutputTree->Branch("mcParticlePx", &myEventData.mcParticlePx);
  myOutputTree->Branch("mcParticlePy", &myEventData.mcParticlePy);
  myOutputTree->Branch("mcParticlePz", &myEventData.mcParticlePz);
  myOutputTree->Branch("mcParticlePDGid", &myEventData.mcParticlePDGid);
  myOutputTree->Branch("mcParticleCharge", &myEventData.mcParticleCharge);
  myOutputTree->Branch("mcParticleGeneratorStatus", &myEventData.mcParticleGeneratorStatus);
  myOutputTree->Branch("mcParticleParentsSize", &myEventData.mcParticleParentsSize);
  myOutputTree->Branch("mcJetPx", &myEventData.mcJetPx);
  myOutputTree->Branch("mcJetPy", &myEventData.mcJetPy);
  myOutputTree->Branch("mcJetPz", &myEventData.mcJetPz);
  myOutputTree->Branch("mcJetE", &myEventData.mcJetE);
  myOutputTree->Branch("mcJetThrust", &myEventData.mcJetThrust);

  myOutputTree->Branch("mcTotalE", &myEventData.mcTotalE);
  myOutputTree->Branch("mcTotalNeutrinosE", &myEventData.mcTotalNeutrinosE);
  myOutputTree->Branch("mcTotalPx", &myEventData.mcTotalPx);
  myOutputTree->Branch("mcTotalPy", &myEventData.mcTotalPy);
  myOutputTree->Branch("mcTotalPz", &myEventData.mcTotalPz);

  myOutputTree->Branch("minDeltaR", &myEventData.minDeltaR);
  myOutputTree->Branch("matchedMCIdx", &myEventData.matchedMCIdx);

  return StatusCode::SUCCESS;
}

/*
 * ================================================
 *
 * Функция обработки
 *
 *        Вызывается для каждого события
 *
 * ================================================
 */
StatusCode myAnalysis::execute()
{
  /* Локальные переменные нужные для работы метода */
  bool hasTof, hasDqdx;
  hasTof = hasDqdx = true;
  TLorentzVector totalPFO4Momentum = TLorentzVector(0., 0., 0., 0.);

  /* Итераторы цикла */
  int pfoIndex = 0;

  /* 
   * Для джет кластеринга создаем вектор, элементами которого 
   * будут являться PseudoJet. PseudoJet это класс fastjet,
   * который представляет 4-вектор.
   *
   * Конструктор у него такой: PseudoJet(px, py, pz, E)
   *
   * Этот вектор является 4-импульсами входных частиц
   * (которые мы будем далее обрабатывать).
   */
  std::vector<fastjet::PseudoJet> inputParticlesPseudoJet;

  /* Извлекаем обязательные коллекции (или ловим исключения) */
  if (getPfoCollection() == false || getMcCollection() == false)
    return StatusCode::SUCCESS;

  /* Извлекаем необязательные коллекции */
  hasTof = getTofCollection();
  hasDqdx = getDqdxCollection();

  /* Перед обработкой собтия чистим все переменные */
  myEventData.reset();

  /* 
   * === Часть 1. Обработка PFO коллекции === 
   *
   * Сохранение нужных переменных для каждого
   * PFO + идентефикация + изоляция.
   */

  /* 
   * Запускаем цикл по всем элементам коллекции myPfoCollPtr.
   * (по сути это цикл по всем реконструированным частицам).
   * Это какой-то C++ цикл, но суть в том, что на каждой
   * итерации этого цикла переменная pfoCollElement принимает
   * значение следующего элемента из коллекции myPfoCollPtr. 
   */
  for (const auto& pfoCollElement : *myPfoCollPtr) 
  {
    /* Проверяем на валидность, если надо пропускаем итерацию */
    if (isInvalidPFO(pfoCollElement))
      continue;

    /* Считываем все параметры элемента PFO */
    myEventData.pfoE.push_back(pfoCollElement.getEnergy());
    myEventData.pfoPx.push_back(pfoCollElement.getMomentum()[0]);
    myEventData.pfoPy.push_back(pfoCollElement.getMomentum()[1]);
    myEventData.pfoPz.push_back(pfoCollElement.getMomentum()[2]);

    /* 
     * Создаем новый элемент PseudoJet и добавляем его в конец 
     * вектора (с помощтю push_back) inputParticlesPseudoJet. 
     *
     * По сути мы тут переписали просто импульсы и энергию
     * из элемента pfo в элемент PseudoJet.
     */
    inputParticlesPseudoJet.push_back(fastjet::PseudoJet(pfoCollElement.getMomentum()[0],
                                                pfoCollElement.getMomentum()[1],
                                                pfoCollElement.getMomentum()[2],
                                                pfoCollElement.getEnergy()));

    /* 
     * Далее нам надо сохранить индекс добавленного pfo объекта
     * в поле user_index, которое предоставляется классом PseudoJet 
     * (см. документацию) (для восстановления сосдава джета)
     */
    inputParticlesPseudoJet.back().set_user_index(pfoIndex++);

    /* 
     * Вычисляем суммарный 4-импульс всех pfo в событии 
     * (тут я создаю временный объект и прибавляю его)
     */
    totalPFO4Momentum += TLorentzVector(pfoCollElement.getMomentum()[0], 
                                        pfoCollElement.getMomentum()[1], 
                                        pfoCollElement.getMomentum()[2], 
                                        pfoCollElement.getEnergy());

    /* Выполняем идентефикацию частиц по данным из Dqdx */
    if (hasDqdx)
      doDqdxParticleIdentification(pfoCollElement);

    /* Выполняем идентефикацию частиц по данным из TOF */
    if (hasTof)
      doTOFParticleIdentification(pfoCollElement);

    /* 
     * TODO Далее в оригинате тут считается энергия и количество
     * кластеров в электромагнитном и адронном калориметрах,
     * что тоже может помочь идентефициаровать частицу, но я
     * пока что не буду это считать, хотя это можно потом добавить
     */

    /* 
     * TODO строка, которая извлекает из реконструкции PDGID, 
     * но я пока не понял как надо добавить ее после того как 
     * я это пойму (сам код извлекает некий getType())
     */

    /* 
     * Вычисляем изоляцию частицы. 
     *
     * Делаем это только при условии, что поперечный импульс больше 2 ГэВ.
     */
    calculateIsolationForPFO(pfoCollElement, MAX_DELTA_R_FOR_ISOLATION);
  }

  /* Сохраняем суммарный импульс всех PFO в событии */
  myEventData.pfoTotalE = totalPFO4Momentum.E();
  myEventData.pfoTotalPx = totalPFO4Momentum.Px();
  myEventData.pfoTotalPy = totalPFO4Momentum.Py();
  myEventData.pfoTotalPz = totalPFO4Momentum.Pz();

  /* === Часть 2. Джет кластеринг === */

  /* Выбираем алгоритм для джет кластеринга */
  fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm);

  /* 
   * Запускаем алгоритм, извлекаем джеты отсортированные по pt 
   *
   * NOTE: в MissingET тут еще дополнительная проверка на количество
   * частиц в inputParticlesPseudoJet, но я пока ее убрал, потому
   * что в документации в fastjet ее не было, пока попробую так.
   *
   * NOTE: есть два доступных способа извлечения джетов, exclusive 
   * и inclusive. На ee ускорителе почему-то используется первый.
   * Про это более подробно написано в документации (см. пункт 3.3.2)
   *
   * Для exclusive есть два способа получения:
   * - С заданием порога (dcut) (не знаю что это, нужно смотреть документацю, 
   *   вроде что-то связанное с ограничением расстояния)
   * - С заданием желаемого числа джетов (пока что используем его)
   */
  fastjet::ClusterSequence cs(inputParticlesPseudoJet, jet_def);
  std::vector<fastjet::PseudoJet> reconstructedJets = 
                    fastjet::sorted_by_pt(cs.exclusive_jets(myNumberJets));

  /* Функция обрабатвает результаты джет кластеринга */
  processJetClusteringResults(reconstructedJets, true);

  /* === Часть 3. Обработка MC коллекции === */

  /* Суммарный импульс MC частиц в событии */
  TLorentzVector totalMC4Momentum = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector totalMCNeutrinos4Momentum = TLorentzVector(0., 0., 0., 0.);

  /* Для MC джет кластеринга */
  std::vector<fastjet::PseudoJet> inputParticlesPseudoJetMC;

  /* Запускаем цикл по MC коллекции*/
  for (const auto& mcCollElement : *myMcPartCollPtr) 
  {
    /* Считываем все параметры элемента mc */
    myEventData.mcParticleE.push_back(mcCollElement.getEnergy());
    myEventData.mcParticlePx.push_back(mcCollElement.getMomentum()[0]);
    myEventData.mcParticlePy.push_back(mcCollElement.getMomentum()[1]);
    myEventData.mcParticlePz.push_back(mcCollElement.getMomentum()[2]);
    myEventData.mcParticlePDGid.push_back(mcCollElement.getPDG());
    myEventData.mcParticleCharge.push_back(mcCollElement.getCharge());
    myEventData.mcParticleGeneratorStatus.push_back(mcCollElement.getGeneratorStatus());
    myEventData.mcParticleParentsSize.push_back(mcCollElement.parents_size());

    /* 
     * 1 - стабильная конечная частица (детектируемая)
     * 2 - промежуточная нестабильная (распадается)
     * (другие коды для всяких технических вещей)
     */
    if (mcCollElement.getGeneratorStatus() != 1) continue;

    /* Создаем 4 вектор для текущей частицы */
    TLorentzVector mc4Momentum = TLorentzVector(mcCollElement.getMomentum()[0], 
                                                mcCollElement.getMomentum()[1],
                                                mcCollElement.getMomentum()[2],
                                                mcCollElement.getEnergy());
    totalMC4Momentum += mc4Momentum; /* Считаем суммарный импульс */

    /* 
     * Если частица это нейтрино, то добавляем 
     * энергию к нейтрино для того чтобы 
     * посчитать недостающую энергию от нейтрино.
     */
    if (abs(mcCollElement.getPDG()) == 12 || 
        abs(mcCollElement.getPDG()) == 14 || 
        abs(mcCollElement.getPDG()) == 16)
    {
      totalMCNeutrinos4Momentum += mc4Momentum;
    } 
    else /* Добавляем частицу в псевдо джет для джет кластеринга MC */
    {
      inputParticlesPseudoJetMC.push_back(fastjet::PseudoJet(mc4Momentum.Px(), 
                                                             mc4Momentum.Py(), 
                                                             mc4Momentum.Pz(), 
                                                             mc4Momentum.E()));
    }
  }

  /* Записываем суммарные значения */
  myEventData.mcTotalE = totalMC4Momentum.E();
  myEventData.mcTotalNeutrinosE = totalMCNeutrinos4Momentum.E();
  myEventData.mcTotalPx = totalMC4Momentum.Px();
  myEventData.mcTotalPy = totalMC4Momentum.Py();
  myEventData.mcTotalPz = totalMC4Momentum.Pz();

  /* Теперь хотим сделать сопоставление реконструкции с MC */
  doGenMatch(inputParticlesPseudoJet, inputParticlesPseudoJetMC);

  /* === Часть 4. Джет кластеринг для MC === */

  fastjet::JetDefinition mc_jet_def(fastjet::ee_kt_algorithm);
  fastjet::ClusterSequence mc_cs(inputParticlesPseudoJetMC, mc_jet_def);
  std::vector<fastjet::PseudoJet> mcJets = 
                    fastjet::sorted_by_pt(mc_cs.exclusive_jets(myNumberJets));

  /* Функция обрабатвает результаты джет кластеринга */
  processJetClusteringResults(mcJets, false);

  /* Успешно завершаемся */
  myEventData.eventNumber++;
  myOutputTree->Fill();
  return StatusCode::SUCCESS;
}

/*
 * ================================================
 *
 * Функция завершения процессора
 *
 *        Вызывается в конце один раз
 *
 * ================================================
 */
StatusCode myAnalysis::finalize()
{
  /* Успешно завершаемся */
  myOutputFile->Write();
  myOutputFile->Close();
  return StatusCode::SUCCESS;
}

/*
 * ================================================
 *
 * Вспомогательные функции (utils)
 *
 * ================================================
 */

/* Функция для проверки инвалидности элемента PFO */
bool myAnalysis::isInvalidPFO(const edm4hep::ReconstructedParticle& pfo) const 
{
  return (std::isnan(pfo.getMomentum()[0]) || 
          std::isnan(pfo.getMomentum()[1]) || 
          std::isnan(pfo.getMomentum()[2]) || 
          pfo.getEnergy() == 0);
}

/* 
 * Функция для идентефикации частиц по данным из коллекции dqdx
 */
void myAnalysis::doDqdxParticleIdentification(const edm4hep::ReconstructedParticle& pfoElement)
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
    for (auto dqdxElement : *myDqdxCollPtr) 
    {
      if (dqdxElement.getTrack() == trk) /* Тоже возвращает один трек, но 
                                          * тут другое название метода почему-то */
      {
        for (size_t i = 0; i < NUMBER_OF_PARTICLE_HYPOTHESES; i++)
          chi2Tpc[i] = dqdxElement.getHypotheses(i).chi2;

        break; /* Завершаем цикл если нашли нужный трек */
      }
    }
  }

  /* 
   * Запысываем резуьтат для последнего трека PFO 
   * (FIXME нужно чтобы не только для последнего было) 
   */
  myEventData.pfoChi2TPC.push_back(chi2Tpc);
}

/* 
 * Функция для идентефикации частиц по коллекциям dqdx и tof 
 */
void myAnalysis::doTOFParticleIdentification(const edm4hep::ReconstructedParticle& pfoElement)
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
    for (auto tofElement : *myTofCollPtr)
    {
      /* 
       * TODO тут странная формула для хи квадрат была, я поменял знаменатель,
       * но возможно это неправильно 
       */
      if (tofElement.getTrack() == trk)
      {
        for (size_t i = 0; i < NUMBER_OF_PARTICLE_HYPOTHESES; i++)
          chi2Tof[i] = std::pow((tofElement.getTime() - tofElement.getTimeExp()[i]) / 
                                tofElement.getTimeExp()[i], 2);

        break; /* Завершаем цикл если нашли нужный трек */
      }
    }
  }

  /* 
   * Запысываем резуьтат для последнего трека PFO 
   * (FIXME нужно чтобы не только для последнего было) 
   */
  myEventData.pfoChi2TOF.push_back(chi2Tof);
}

/* 
 * Функция для вычисления изоляции PFO.
 */
void myAnalysis::calculateIsolationForPFO(const edm4hep::ReconstructedParticle& pfoElement,
                                          double deltaR)
{
  /* Создаем 4-вектор для переданного PFO объекта */
  TLorentzVector pfo4Momentum;
  TLorentzVector pfo4MomentumNeighbour;

  /* Суммарный импульс */
  double totalMomentumCharge = 0; 
  double totalMomentumNeutral = 0;

  /* Заполняем 4 вектор для обрабатываемого pfo объекта */
  pfo4Momentum.SetPxPyPzE(pfoElement.getMomentum()[0], 
                          pfoElement.getMomentum()[1], 
                          pfoElement.getMomentum()[2], 
                          pfoElement.getEnergy());

  /* 
   * Для частиц у которых поперечный импульс
   * меньше 2 ГэВ мы считаем изоляцию равной -1.
   * (по сути просто не считаем)
   */
  if (pfo4Momentum.Pt() < MIN_PT_MOMENTUM_FOR_ISOLATION)
  {
    myEventData.relativeIsolationCharged.push_back(-1);
    myEventData.relativeIsolationNeutral.push_back(-1);
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
        &pfoCollElementNeighbour == &pfoElement)
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
        totalMomentumNeutral += pfo4MomentumNeighbour.P();
      else
        totalMomentumCharge += pfo4MomentumNeighbour.P();
    }
  }

  /* 
   * Теперь нормируем изоляцию на импульс самого pfo 
   * и запысываем результат в eventData 
   */
  myEventData.relativeIsolationCharged.push_back(totalMomentumCharge / pfo4Momentum.P());
  myEventData.relativeIsolationNeutral.push_back(totalMomentumNeutral / pfo4Momentum.P());
}

/* 
 * Функция для обработки результатов джет кластеринга.
 *
 * Она предназначена как для генераторных джетов, так и для реконструированных.
 */
void myAnalysis::processJetClusteringResults(const std::vector<fastjet::PseudoJet>& jets, 
                                             bool isReconstructed)
{
  /* Запускаем цикл по полученным джетам чтобы извлечь их параметры */
  for (size_t i = 0; i < jets.size(); i++) 
  {
    /* Создаем 4импульс джета */
    TLorentzVector p4jet;
    p4jet.SetPxPyPzE(jets[i].px(), jets[i].py(), 
                     jets[i].pz(), jets[i].E());

    /* Извлекаем все составляющие джет псевдоджеты */
    std::vector<fastjet::PseudoJet> jetConstituents = jets[i].constituents();

    /* Запускаем цикл по всем составляющим джета */
    double jetThrust = 0.;                /* Мера коллимированности струи */
    double jetPSum = 0.;                  /* Сумма модулей импульсов составляющих джета */
    std::vector<int> jetConstituentsIdx;  /* Индексы PFO, которые составляют джет */
    for (size_t j = 0; j < jetConstituents.size(); j++)
    {
      /* Добавляем индекс составляющей джета в вектор его составляющих */
      jetConstituentsIdx.push_back(jetConstituents[j].user_index());

      /* Создаем 4импульс этой составляющей */
      TLorentzVector p4jetConstituents;
      p4jetConstituents.SetPxPyPzE(jetConstituents[j].px(), 
                                   jetConstituents[j].py(),
                                   jetConstituents[j].pz(), 
                                   jetConstituents[j].E());

      /* Вычисляем jetThrust */
      jetThrust += abs(p4jetConstituents.Px() * p4jet.Px() +
                       p4jetConstituents.Py() * p4jet.Py() +
                       p4jetConstituents.Pz() * p4jet.Pz()) / p4jet.P();
      jetPSum += p4jetConstituents.P(); /* Для нормировки jetThrust */

      /* 
       * Вычисляем Energy Correlation Functions (ECF) 
       * TODO решил пока не делать, нужно это обсудить с М.В.
       */
    }

    /* Делаем нормировку jetThrust */
    jetThrust = jetThrust / jetPSum;

    /* 
     * TODO также тут считаются N-subjettiness (τ₁, τ₂, τ₃),
     * решил пока что тоже не делать, потому что нужно это обсудить.
     */
    
    /* Записываем все параметры джета в eventData */
    if (isReconstructed) /* Если джеты из реконструкции */
    {
      myEventData.reconstructedJetConstituentsPfoIdx.push_back(jetConstituentsIdx);
      myEventData.reconstructedJetPx.push_back(jets[i].px());
      myEventData.reconstructedJetPy.push_back(jets[i].py());
      myEventData.reconstructedJetPz.push_back(jets[i].pz());
      myEventData.reconstructedJetE.push_back(jets[i].E());
      myEventData.reconstructedJetThrust.push_back(jetThrust);
    }
    else /* Если это MC */ 
    {
      myEventData.mcJetPx.push_back(jets[i].px());
      myEventData.mcJetPy.push_back(jets[i].py());
      myEventData.mcJetPz.push_back(jets[i].pz());
      myEventData.mcJetE.push_back(jets[i].E());
      myEventData.mcJetThrust.push_back(jetThrust);
    }
  }
}

/* 
 * Функция для сопоставления реконструкции с MC.
 * inputParticles - должен быть созан на исходных PFO.
 */
void myAnalysis::doGenMatch(const std::vector<fastjet::PseudoJet>& inputParticles, 
                            const std::vector<fastjet::PseudoJet>& inputParticlesMC)
{

  /* Запускаем цикл по всем реконструированным частицам */
  for (size_t i = 0; i < inputParticles.size(); i++) 
  {
    /* Получаем текущую частицу */
    fastjet::PseudoJet inputParticle = inputParticles[i];

    /* Индекс совпавшей частицы */
    int matchedMCIdx = -1;

    /* Минимальное расстояние */
    double minDeltaR = 999.0;

    /* Ищем ближайшую частицу */
    for (size_t j = 0; j < inputParticlesMC.size(); j++) 
    {
      /* Получаем текущую MC частицу */
      fastjet::PseudoJet inputParticleMC = inputParticlesMC[j];

      /* Записываем эту частицу, если расстояние меньше */
      if (inputParticle.delta_R(inputParticleMC) < minDeltaR)
      {
        matchedMCIdx = j;
        minDeltaR = inputParticle.delta_R(inputParticleMC);
      }
    }

    /* Записываем результаты */
    myEventData.minDeltaR.push_back(minDeltaR);
    myEventData.matchedMCIdx.push_back(matchedMCIdx);
  }
}

/*
 * Функции для извлечения коллекций
 */
bool myAnalysis::getPfoCollection()
{
  try 
  {
    myPfoCollPtr = pfoCollHandler.get();
  }
  catch (const GaudiException& e) 
  {
    info() << "Collection " << pfoCollHandler.fullKey() 
           << " is unavailable in event " << myEventData.eventNumber << endmsg;
    return false;
  }
  
  if (myPfoCollPtr->empty()) 
  {
    info() << "Collection " << pfoCollHandler.fullKey() 
           << " is empty in event " << myEventData.eventNumber << endmsg;
    return false;
  }
  
  return true;
}

bool myAnalysis::getMcCollection()
{
  try 
  {
    myMcPartCollPtr = mcCollHandler.get();
  }
  catch (const GaudiException& e) 
  {
    info() << "Collection " << mcCollHandler.fullKey() 
           << " is unavailable in event " << myEventData.eventNumber << endmsg;
    return false;
  }
  
  if (myMcPartCollPtr->empty()) 
  {
    info() << "Collection " << mcCollHandler.fullKey() 
           << " is empty in event " << myEventData.eventNumber << endmsg;
    return false;
  }
  
  return true;
}

bool myAnalysis::getTofCollection()
{
  try 
  {
    myTofCollPtr = tofCollHandler.get();
    return true;
  }
  catch (const GaudiException& e) 
  {
    info() << "Collection " << tofCollHandler.fullKey() 
           << " is unavailable in event " << myEventData.eventNumber << endmsg;
    return false;
  }
}

bool myAnalysis::getDqdxCollection()
{
  try 
  {
    myDqdxCollPtr = dqdxCollHandler.get();
    return true;
  }
  catch (const GaudiException& e) 
  {
      info() << "Collection " << dqdxCollHandler.fullKey() 
             << " is unavailable in event " << myEventData.eventNumber << endmsg;
      return false;
  }
}
