/*
 * Описание того что это
 * 
 * Описание терминов:
 * PFO - по сути это частица, которую восстановил алгоритм
 *       реконструкции, он не знает что это именно, но записал параметры.
 */





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
// #include "fastjet/contrib/Nsubjettiness.hh"

#include <chrono>

#include "TLorentzVector.h"

/*
 * Какой-то макрос из Gaudi, который нужен для создания
 * фабрики, которая будет создавать объекты моего класса
 */
DECLARE_COMPONENT(myAnalysis)

/*
 * Вызываем конструктор класса
 * Вот это : Algorithm( name, pSvcLocator ) вызывает конструктор родительского 
 * класса Algorithm перед выполнением тела нашего конструктора.
 * name - Имя нашего алгоритма
 * pSvcLocator - указатель на навигатор сервисов
 *
 * Тут (внутри тела конструктора) еще можно делать некие declareProperty(), но я пока что не буду.
 * Как я понял это может понадобиться если мы хотим чтобы можно было менять входные данные из python.
 */
myAnalysis::myAnalysis(const string& name, ISvcLocator* pSvcLocator ) : Algorithm( name, pSvcLocator ) { }

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
  eventNumber = 0;

  /* Суммарное количество частиц */
  numberParticles = 0;

  /* Создаем и открываем выходной файл + дерево */
  myOutputFile = TFile::Open(m_outputFile.value().c_str(), "RECREATE");
  myOutputTree = new TTree("outputTree", "outputTree");

  /* Создаем все ветки для заполнения данными из анализа */
  /* TODO тут надо это создать в зависимости от того что я буду считать */

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

  /* Итераторы цикла */
  int pfoIndex = 0;

  /* Обязательные коллекции */
  edm4hep::ReconstructedParticleCollection *myPfoCollPtr = NULL;
  edm4hep::MCParticleCollection *myMcPartCollPtr = NULL;

  /* Необязательные коллекции */
  edm4hep::RecTofCollection     *myTofCollPtr = NULL;
  edm4hep::RecDqdxCollection    *myDqdxCollPtr = NULL;

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
  vector<fastjet::PseudoJet> inputParticlesPseudoJet;

  /* Извлекаем обязательные коллекции (или ловим исключения) */
  if (getCollection(pfoCollHandler, &myPfoCollPtr, eventNumber) == false || 
      getCollection(mcCollHandler, &myMcPartCollPtr, eventNumber) == false)
  {
    return StatusCode::SUCCESS;
  }

  /* Извлекаем необязательные коллекции */
  hasTof = getCollection(tofCollHandler, &myTofCollPtr, eventNumber);
  hasDqdx = getCollection(dqdxCollHandler, &myDqdxCollPtr, eventNumber);

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

    /* 
     * Создаем новый элемент PseudoJet и добавляем его в конец 
     * вектора (с помощтю push_back) inputParticlesPseudoJet. 
     *
     * По сути мы тут переписали просто импульсы и энергию
     * из элемента pfo в элемент PseudoJet.
     */
    inputParticlesPseudoJet.push_back(PseudoJet(pfoCollElement.getMomentum()[0],
                                                pfoCollElement.getMomentum()[1],
                                                pfoCollElement.getMomentum()[2],
                                                pfoCollElement.getEnergy()));

    /* 
     * Далее нам надо сохранить индекс добавленного pfo объекта
     * в поле user_index, которое предоставляется классом PseudoJet 
     * (см. документацию) (TODO для чего не очень понятно, но видимо 
     * потом для восстановления сосдава джета)
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
      doDqdxParticleIdentification(myDqdxCollPtr, pfoCollElement);

    /* Выполняем идентефикацию частиц по данным из TOF */
    if (hasTof)
      doTOFParticleIdentification(myDqdxCollPtr, pfoCollElement);

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
    calculateIsolationForPFO(myPfoCollPtr, pfoCollElement,
                             myEventData, MAX_DELTA_R_FOR_ISOLATION);
  }

  /* 
   * TODO Далее идет сохранение суммарного импульса pfo
   * по всем компонентам в выходное дерево (просто сохраняют 
   * значения из totalPFO4Momentum). Я пока что не буду это 
   * делать так как не решил еще полностью какие переменные
   * мы будем сохранять для каждого события
   */

  /* === Часть 2. Джет кластеринг === */

  /* Выбираем алгоритм для джет кластеринга */
  JetDefinition jet_def(ee_kt_algorithm);

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
  ClusterSequence cs(inputParticlesPseudoJet, jet_def);
  myEventData.jets = sorted_by_pt(cs.exclusive.jets(myNumberJets));

  /* Функция обрабатвает результаты джет кластеринга */
  processJetClusteringResults(jets, false, myEventData);





  

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
  
  return StatusCode::SUCCESS;
}

/*
 * ================================================
 *
 * Вспомогательные функции (utils)
 *
 * ================================================
 */



























