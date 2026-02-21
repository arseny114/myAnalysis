// src/myAnalysis.cpp (updated)

/**
 * @file myAnalysis.cpp
 * @brief Реализация алгоритма myAnalysis для обработки событий в CEPCSW
 *
 * Алгоритм выполняет следующие основные задачи:
 * - чтение реконструированных частиц (PFO)
 * - расчёт изоляции частиц
 * - отсев событий с изолированными лептонами, фотонами и заряженными адронами
 * - кластеризация на джеты с помощью fastjet (теперь инклюзивно по умолчанию)
 * - отсев событий, в которых хотя бы один джет содержит менее 6 частиц
 * - расчёт инвариантной массы и массы отдачи двумя способами
 */

#include "myAnalysis.h"

// Стандартные библиотеки
#include <cmath>
#include <chrono>

// Компоненты Gaudi
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

// Интерфейсы детектора и единицы измерения
#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"

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
myAnalysis::myAnalysis(const std::string& name, ISvcLocator* pSvcLocator)
    : Algorithm(name, pSvcLocator)
{
    // Старые свойства
    declareProperty("numberJets",              myNumberJets,              "Желаемое количество джетов");
    declareProperty("isolationDeltaR",         myIsolationDeltaR,         "Радиус конуса для расчёта изоляции (ΔR)");
    declareProperty("outputRootFile",          myOutputFileName,          "Имя выходного ROOT-файла");
    declareProperty("centerOfMassEnergy",      myCenterOfMassEnergy,      "Полная энергия в системе центра масс (ГэВ)");

    // Новые свойства для ee_genkt_algorithm
    declareProperty("jetR",                    myJetR,                    "Радиус джета R для ee_genkt_algorithm");
    declareProperty("jetP",                    myJetP,                    "Параметр p для ee_genkt_algorithm (1 для kt-like)");
    declareProperty("jetPtMin",                myJetPtMin,                "Минимальный pT джета для inclusive режима (GeV)");
    declareProperty("useInclusive",            myUseInclusive,            "Использовать inclusive кластеризацию (true) или exclusive (false)");
    declareProperty("pfoEnergyMin",            myPfoEnergyMin,            "Минимальная энергия PFO для кластеризации джетов (GeV)");
    declareProperty("minPtForIsolation",       myMinPtForIsolation,       "Минимальный поперечный импульс для расчёта изоляции (ГэВ)");
    declareProperty("isolationThreshold",      myIsolationThreshold,      "Порог изоляции, ниже которого частица считается изолированной");
    declareProperty("minConstPerJet",          myMinConstPerJet,          "Минимальное количество частиц в джете для принятия события");

    declareProperty("applyJetSelection",       myApplyJetSelection,       "Применять отбор по джетам (true) или сохранять все (false)");
    declareProperty("applyIsolationSelection", myApplyIsolationSelection, "Применять отбор по изоляции (true) или сохранять все (false)");
}

/**
 * @brief Инициализация алгоритма (вызывается один раз в начале)
 * 
 * Создаёт структуру данных события, открывает выходной файл,
 * создаёт дерево и все необходимые ветви.
 */
StatusCode myAnalysis::initialize()
{
    myEventData = EventData();
    myEventData.eventNumber = 0;

    // Открываем файл в режиме пересоздания
    myOutputFile = TFile::Open(myOutputFileName.value().c_str(), "RECREATE");
    if (!myOutputFile || myOutputFile->IsZombie()) {
        error() << "Не удалось создать выходной файл: " << myOutputFileName << endmsg;
        return StatusCode::FAILURE;
    }

    myOutputTree = new TTree("outputTree", "Results");

    // Ветки для отдельных PFO
    myOutputTree->Branch("pfoE",   &myEventData.pfoE);
    myOutputTree->Branch("pfoPx",  &myEventData.pfoPx);
    myOutputTree->Branch("pfoPy",  &myEventData.pfoPy);
    myOutputTree->Branch("pfoPz",  &myEventData.pfoPz);

    // Суммарные характеристики всех PFO в событии
    myOutputTree->Branch("pfoTotalE",  &myEventData.pfoTotalE);
    myOutputTree->Branch("pfoTotalPx", &myEventData.pfoTotalPx);
    myOutputTree->Branch("pfoTotalPy", &myEventData.pfoTotalPy);
    myOutputTree->Branch("pfoTotalPz", &myEventData.pfoTotalPz);

    // Изоляция
    myOutputTree->Branch("relativeIsolation", &myEventData.relativeIsolation);

    // Какие типы частиц были найдены
    myOutputTree->Branch("particleType", &myEventData.particleType);
    myOutputTree->Branch("isLepton", &myEventData.isLepton);
    myOutputTree->Branch("isChargedHadron", &myEventData.isChargedHadron);

    // Размеры джетов и количество найденных джетов в событии
    myOutputTree->Branch("jetSize", &myEventData.jetSize);
    myOutputTree->Branch("numberJetsInEvent", &myEventData.numberJetsInEvent);

    // Характеристики реконструированных джетов
    myOutputTree->Branch("reconstructedJetConstituentsPfoIdx", &myEventData.reconstructedJetConstituentsPfoIdx);
    myOutputTree->Branch("reconstructedJetPx", &myEventData.reconstructedJetPx);
    myOutputTree->Branch("reconstructedJetPy", &myEventData.reconstructedJetPy);
    myOutputTree->Branch("reconstructedJetPz", &myEventData.reconstructedJetPz);
    myOutputTree->Branch("reconstructedJetE",  &myEventData.reconstructedJetE);
    myOutputTree->Branch("reconstructedJetThrust", &myEventData.reconstructedJetThrust);

    // Физические величины, которые нас интересуют
    myOutputTree->Branch("invariantMassAllPFO", &myEventData.invariantMassAllPFO);
    myOutputTree->Branch("invariantMassJets", &myEventData.invariantMassJets);
    myOutputTree->Branch("recoilMassAllPFO",    &myEventData.recoilMassAllPFO);
    myOutputTree->Branch("recoilMassJets",    &myEventData.recoilMassJets);

    return StatusCode::SUCCESS;
}

/**
 * @brief Основная функция обработки каждого события
 */
StatusCode myAnalysis::execute()
{
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
    for (const auto& pfo : *myPfoCollPtr)
    {
        // Пропускаем невалидные
        if (isInvalidPFO(pfo)) continue;

        // Увеличиваем индекс при обработке нового PFO
        pfoIndex++;

        // Сохраняем кинематику
        myEventData.pfoE.push_back(pfo.getEnergy());
        myEventData.pfoPx.push_back(pfo.getMomentum()[0]);
        myEventData.pfoPy.push_back(pfo.getMomentum()[1]);
        myEventData.pfoPz.push_back(pfo.getMomentum()[2]);

        // Суммируем полный четырёхимпульс события
        totalPFO4Momentum += TLorentzVector(
            pfo.getMomentum()[0], pfo.getMomentum()[1],
            pfo.getMomentum()[2], pfo.getEnergy()
        );

        // Если частица слишком мягкая, то не будем считать для нее
        // изоляцию и она не будет участвовать в джет кластеринге
        if (pfo.getEnergy() < myPfoEnergyMin.value()) 
        {
            // Для слишком мягких частиц изоляцию не считаем 
            myEventData.relativeIsolation.push_back(-1.0);
            continue;
        }

        // Готовим частицу для кластеризации
        // emplace_back() это эффективный способ добавления в конец вектора
        // (он также создает объект PseudoJet с переданными параметрами импульса)
        inputParticles.emplace_back(
            pfo.getMomentum()[0], pfo.getMomentum()[1],
            pfo.getMomentum()[2], pfo.getEnergy()
        );
        // Сохраняем для восстановления состава джета
        inputParticles.back().set_user_index(pfoIndex);

        // Считаем изоляцию
        calculateIsolationForPFO(pfo, myIsolationDeltaR);
    }

    // Сохраняем суммарный четырёхимпульс всех PFO в событии
    myEventData.pfoTotalE  = totalPFO4Momentum.E();
    myEventData.pfoTotalPx = totalPFO4Momentum.Px();
    myEventData.pfoTotalPy = totalPFO4Momentum.Py();
    myEventData.pfoTotalPz = totalPFO4Momentum.Pz();

    // ── 2. Проверка на изолированные объекты ────────────────────────────

    // Проходим по всем сохранённым PFO (реконструированным частицам) события
    for (size_t i = 0; i < myEventData.relativeIsolation.size(); ++i)
    {
        // Восстанавливаем четырёхимпульс i-й частицы из сохранённых значений
        // (мы их ранее записали в векторы pfoPx, pfoPy, pfoPz, pfoE)
        TLorentzVector p4(
            myEventData.pfoPx[i],
            myEventData.pfoPy[i],
            myEventData.pfoPz[i],
            myEventData.pfoE[i]
        );

        // Пропускаем частицы, для которых изоляция не вычислялась
        if (myEventData.relativeIsolation[i] == -1.0) continue;

        // Получаем тип частицы
        int pfoType   = (*myPfoCollPtr)[i].getType();

        int absPdg    = std::abs(pfoType);                  // модуль типа
        double relIso = myEventData.relativeIsolation[i];   // ранее посчитанная изоляция

        // Классифицируем частицу по типу
        bool isLepton        = (absPdg == 11 || absPdg == 13);                        // e⁻/e⁺, μ⁻/μ⁺
        bool isChargedHadron = (absPdg == 2212 || absPdg == 321 || absPdg == 211);    // заряженные адроны

        // Запоминмаем тип этой частицы и что это была за частица
        myEventData.particleType.push_back(pfoType);
        myEventData.isLepton.push_back(isLepton ? 1 : 0);
        myEventData.isChargedHadron.push_back(isChargedHadron ? 1 : 0);

        // Условие отбора события:
        // Если частица изолирована (мало энергии/импульса в конусе вокруг неё)
        // И при этом она относится к одной из «нежелательных» категорий, то 
        // отбрасываем это событие
        if (myApplyIsolationSelection.value()) // отбрасываем события, только если включен режим отбрасывания
        {
            if (relIso > 0 && relIso < myIsolationThreshold.value() && (isLepton || isChargedHadron))
            {
                return StatusCode::SUCCESS;
            }
        }
    }

    // ── 3. Кластеризация джетов ───────────────────────────────────
    fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, myJetR.value(), myJetP.value());
    fastjet::ClusterSequence cs(inputParticles, jet_def);

    std::vector<fastjet::PseudoJet> jets;

    // Выбор режима: exclusive или inclusive
    if (myUseInclusive.value()) 
    {
        // Inclusive режим: джеты с pT > myJetPtMin, может быть >2 джетов, с неприсвоенными частицами
        jets = fastjet::sorted_by_pt(cs.inclusive_jets(myJetPtMin.value()));

        // Записываем колличество найденных джетов в событии
        myEventData.numberJetsInEvent = jets.size();

        if (myApplyJetSelection.value()) // отбрасываем события, только если включен режим отбрасывания
        {
            // Отсев, если нашлось не 2 джета (TODO: этот отсев нужно переделать)
            if (jets.size() != 2) 
            {
                return StatusCode::SUCCESS;
            }
        }
    } 
    else 
    {
        // Exclusive режим: ровно myNumberJets джетов, все частицы присвоены
        jets = fastjet::sorted_by_pt(cs.exclusive_jets(myNumberJets.value()));
    }

    // Проверка минимального количества частиц во всех джетах
    for (const auto& jet : jets) 
    {
        // Запоминаем размер джета
        myEventData.jetSize.push_back(jet.constituents().size());

        if (myApplyJetSelection.value()) // отбрасываем события, только если включен режим отбрасывания
        {
            // Отсев, если в размер джетов меньше требуемого
            if (jet.constituents().size() < myMinConstPerJet.value()) 
            {
                return StatusCode::SUCCESS;
            }
        }
    }

    // Сохраняем информацию о джетах
    saveJetClusteringResults(jets);

    // ── 4. Расчёт инвариантных масс и масс отдачи ───────────────────────
    myEventData.invariantMassAllPFO = totalPFO4Momentum.M();
    if (myEventData.invariantMassAllPFO > myCenterOfMassEnergy) // отбрасываем ошибки реконструкции (TODO: почему они возникли?)
      return StatusCode::SUCCESS;

    // Создаём четырёхимпульсы ведущих джетов
    // Для invariantMassJets и recoilMassJets теперь суммируем все джеты
    TLorentzVector summedJets(0., 0., 0., 0.);
    for (const auto& jet : jets) {
        summedJets += TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.E());
    }

    // Инвариантная масса системы джетов (масса всех джетов вместе)
    // Это одна из ключевых наблюдаемых величин в анализе
    myEventData.invariantMassJets = summedJets.M();
    if (myEventData.invariantMassJets > myCenterOfMassEnergy) // отбрасываем ошибки реконструкции
          return StatusCode::SUCCESS;

    // Энергия в системе центра масс (обычно 240 ГэВ для CEPC при Z-фабрике)
    double sqrts = myCenterOfMassEnergy;

    // ── Расчёт массы отдачи по всем PFO (реконструированная полная система) ──
    
    // Энергия отдачи = полная энергия события минус энергия всех PFO
    double recoilE_all = sqrts - myEventData.pfoTotalE;

    // Квадрат полного импульса системы всех PFO
    double recoilP2_all = std::pow(myEventData.pfoTotalPx, 2) +
                          std::pow(myEventData.pfoTotalPy, 2) +
                          std::pow(myEventData.pfoTotalPz, 2);

    // Масса отдачи = √(E_recoil² - p_recoil²)
    // В идеале должна быть близка к массе Z-бозона (~91 ГэВ) в процессе e⁺e⁻ → ZH
    myEventData.recoilMassAllPFO = std::sqrt(recoilE_all * recoilE_all - recoilP2_all);

    // ── Расчёт массы отдачи по джетам (система всех джетов) ──

    // Энергия отдачи = полная энергия минус энергия системы джетов
    double recoilE_jets = sqrts - summedJets.E();
    double recoilP2_jets = std::pow(summedJets.Px(), 2) + std::pow(summedJets.Py(), 2) + std::pow(summedJets.Pz(), 2);
    myEventData.recoilMassJets = std::sqrt(recoilE_jets * recoilE_jets - recoilP2_jets);

    // Заполняем дерево и увеличиваем счётчик событий
    myEventData.eventNumber++;
    myOutputTree->Fill();

    return StatusCode::SUCCESS;
}

/**
 * @brief Завершение работы алгоритма (вызывается один раз в конце)
 */
StatusCode myAnalysis::finalize()
{
    if (myOutputFile)
    {
        myOutputFile->Write();
        myOutputFile->Close();
    }
    return StatusCode::SUCCESS;
}

/**
 * @brief Проверка, является ли PFO невалидным
 */
bool myAnalysis::isInvalidPFO(const edm4hep::ReconstructedParticle& pfo) const
{
    const auto& mom = pfo.getMomentum();
    return (std::isnan(mom[0]) || std::isnan(mom[1]) || std::isnan(mom[2]) ||
            pfo.getEnergy() <= 0);
}

/**
 * @brief Вычисляет относительную изоляцию данной PFO-частицы
 *
 * Относительная изоляция определяется как сумма модулей импульсов всех других
 * частиц (кроме самой рассматриваемой), находящихся внутри конуса с раствором
 * deltaR вокруг данной частицы, нормированная на модуль импульса самой частицы.
 *
 * Используется суммарная изоляция (учитываются и заряженные, и нейтральные частицы).
 * Частицы с поперечным импульсом Pt меньше порогового значения myMinPtForIsolation.value()
 * считаются не подлежащими изоляционному анализу и получают значение -1.
 *
 * @param pfo Ссылка на объект реконструированной частицы (PFO), для которой
 *            производится расчёт изоляции
 * @param deltaR Угловой радиус конуса изоляции (в единицах ΔR = √(Δη² + Δφ²))
 */
void myAnalysis::calculateIsolationForPFO(const edm4hep::ReconstructedParticle& pfo,
                                          double deltaR)
{
    // Формируем четырёхимпульс рассматриваемой частицы
    TLorentzVector thisP4(
        pfo.getMomentum()[0], pfo.getMomentum()[1],
        pfo.getMomentum()[2], pfo.getEnergy()
    );

    // Если поперечный импульс слишком мал, то изоляцию не считаем
    // (обычно такие частицы не представляют интереса для анализа изоляции)
    if (thisP4.Pt() < myMinPtForIsolation.value())
    {
        myEventData.relativeIsolation.push_back(-1.);
        return;
    }

    // Сумма модулей импульсов всех частиц, попавших в конус изоляции
    double sumP_inCone = 0.0;

    // Проходим по всей коллекции реконструированных частиц события
    for (const auto& other : *myPfoCollPtr)
    {
        // Пропускаем саму рассматриваемую частицу
        if (&other == &pfo) continue;

        // Пропускаем невалидные/повреждённые объекты
        if (isInvalidPFO(other)) continue;

        // Формируем четырёхимпульс соседней частицы
        TLorentzVector otherP4(
            other.getMomentum()[0], other.getMomentum()[1],
            other.getMomentum()[2], other.getEnergy()
        );

        // Проверяем, попадает ли соседняя частица в конус изоляции
        if (thisP4.DeltaR(otherP4) < deltaR)
        {
            // Добавляем модуль импульса соседней частицы к сумме
            sumP_inCone += otherP4.P();
        }
    }

    // Относительная изоляция = сумма импульсов в конусе / импульс центральной частицы
    // Значение сохраняется в вектор для последующей записи в дерево
    myEventData.relativeIsolation.push_back(sumP_inCone / thisP4.P());
}

/**
 * @brief Сохраняет основные характеристики найденных джетов в структуру данных события
 *
 * Для каждого джета функция:
 *  - сохраняет индексы всех составляющих его частиц (PFO) через user_index
 *  - вычисляет параметр "thrust" джета — меру коллимированности (насколько хорошо
 *    импульсы составляющих направлены вдоль оси джета)
 *  - сохраняет четырёхимпульс джета (px, py, pz, E)
 *
 * Thrust здесь считается как средневзвешенный косинус угла между импульсом каждой
 * составляющей и направлением джета (нормированный на модуль импульса джета).
 * Значение thrust близкое к 1 — джет очень коллимированный (струя узкая),
 * близкое к 0 — составляющие разлетаются в разные стороны.
 */
void myAnalysis::saveJetClusteringResults(const std::vector<fastjet::PseudoJet>& jets)
{
    // Проходим по всем найденным джетам
    for (const auto& jet : jets)
    {
        // Вектор для хранения индексов исходных PFO, входящих в данный джет
        std::vector<int> constituentsIdx;

        // Переменные для расчёта thrust
        double thrust = 0.0;     // суммарный проекционный вклад
        double pSum   = 0.0;     // сумма модулей импульсов всех составляющих

        // Получаем список всех частиц, входящих в этот джет
        auto cons = jet.constituents();

        // Проходим по каждой составляющей джета
        for (const auto& c : cons)
        {
            // Сохраняем индекс исходной PFO-частицы (был установлен при создании PseudoJet)
            constituentsIdx.push_back(c.user_index());

            // Четырёхимпульсы составляющей и самого джета
            TLorentzVector c4(c.px(), c.py(), c.pz(), c.E());
            TLorentzVector j4(jet.px(), jet.py(), jet.pz(), jet.E());

            // Вклад в thrust: |p_const · p_jet| / |p_jet|
            // Это проекция импульса частицы на направление джета, нормированная
            // c4.Vect() - трехмерный вектор импульса частицы
            // c4.Vect().Dot(j4.Vect()) - скалярное произведение двух векторов импульсов (берем от него abs)
            thrust += std::abs(c4.Vect().Dot(j4.Vect())) / j4.P();

            // Суммируем модули импульсов для последующей нормировки
            pSum += c4.P();
        }

        // Нормируем thrust на суммарный импульс составляющих
        // Защита от деления на ноль (хотя в реальных джетах pSum почти никогда не бывает нулевым)
        if (pSum > 1e-9)
            thrust /= pSum;
        else
            thrust = 0.0;

        // Сохраняем все характеристики текущего джета в структуру данных события
        myEventData.reconstructedJetConstituentsPfoIdx.push_back(std::move(constituentsIdx));
        myEventData.reconstructedJetPx.push_back(jet.px());
        myEventData.reconstructedJetPy.push_back(jet.py());
        myEventData.reconstructedJetPz.push_back(jet.pz());
        myEventData.reconstructedJetE.push_back(jet.E());
        myEventData.reconstructedJetThrust.push_back(thrust);
    }
}

/**
 * @brief Получение коллекции PFO из хранилища событий
 */
bool myAnalysis::getPfoCollection()
{
    try
    {
        myPfoCollPtr = pfoCollHandler.get();
    }
    catch (const GaudiException& e)
    {
        info() << "Коллекция PFO недоступна в событии " << myEventData.eventNumber << endmsg;
        return false;
    }

    if (myPfoCollPtr->empty())
    {
        info() << "Коллекция PFO пуста в событии " << myEventData.eventNumber << endmsg;
        return false;
    }

    return true;
}
