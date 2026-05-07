// include/zh_invisible_analysis.h
#ifndef ZH_INVISIBLE_ANALYSIS_H
#define ZH_INVISIBLE_ANALYSIS_H

#include <string>

// =============================================================================
// НАСТРОЙКИ ОТБОРА
// =============================================================================

// Вето на изолированные лептоны
#define APPLY_LEPTON_VETO true

// Вето на высокоэнергетические фотоны
#define APPLY_PHOTON_VETO true
#define PHOTON_ENERGY_CUT_GEV 30.0 // Порог энергии фотона для вето

// Требование ровно 2 джета в инклюзивном режиме
#define REQUIRE_EXACTLY_TWO_INCLUSIVE_JETS true

// Требование минимального числа конституентов в каждом из двух джетов
#define REQUIRE_MIN_CONSTITUENTS_PER_JET true
#define MIN_CONSTITUENTS_PER_JET 6

// Окно инвариантной массы диджета (ГэВ)
#define APPLY_DIJET_MASS_WINDOW true
#define DIJET_MASS_WINDOW_MIN_GEV 65.0
#define DIJET_MASS_WINDOW_MAX_GEV 110.0

// Окно массы отдачи (ГэВ)
#define APPLY_RECOIL_MASS_WINDOW true
#define RECOIL_MASS_WINDOW_MIN_GEV 110.0
#define RECOIL_MASS_WINDOW_MAX_GEV 165.0

// =============================================================================
// ПАРАМЕТРЫ ФИЗИКИ И ГИСТОГРАММ
// =============================================================================

const std::string TREE_NAME = "outputTree";

// Массы частиц и энергия столкновения
const double MZ_GEV = 91.2;
const double MH_GEV = 125.26;
const double SQRT_S_GEV = 240.0;

// PDG коды частиц
const int PDG_PHOTON = 22;
const int PDG_ELECTRON = 11;
const int PDG_MUON = 13;

const int MASS_BINS = 200;
const double MASS_MIN_GEV = 0.0;
const double MASS_MAX_GEV = 250.0;

const int RECOIL_BINS = 200;
const double RECOIL_MIN_GEV = 0.0;
const double RECOIL_MAX_GEV = 250.0;

// Параметры для 2D гистограммы: E_photon vs M_recoil
const int PHOTON_E_BINS = 100;
const double PHOTON_E_MIN_GEV = PHOTON_ENERGY_CUT_GEV;
const double PHOTON_E_MAX_GEV = SQRT_S_GEV;

const int RECOIL_MASS_2D_BINS = 200;
const double RECOIL_MASS_2D_MIN_GEV = 0.0;
const double RECOIL_MASS_2D_MAX_GEV = 250.0;

// Параметры для гистограмм угловых распределений
const int THETA_BINS = 100;
const double THETA_MIN_RAD = 0.0;
const double THETA_MAX_RAD = 3.14159; // π

const int DELTA_R_BINS = 100;
const double DELTA_R_MIN = 0.0;
const double DELTA_R_MAX = 5.0;

// =============================================================================
// НАСТРОЙКИ ВЫВОДА И ЛОГИРОВАНИЯ
// =============================================================================

const std::string OUTPUT_BASE_DIR = "../pdf_results";

const int LOG_INTERVAL_EVENTS = 1000;
const bool LOG_PERCENTAGE = true;
const bool PRINT_CUT_STATISTICS = true;

// =============================================================================
// СТРУКТУРА ДЛЯ СТАТИСТИКИ ПРОХОЖДЕНИЯ КАТОВ
// =============================================================================

struct CutStatistics {
    Long64_t totalEvents = 0;
    Long64_t afterLeptonVeto = 0;
    Long64_t afterPhotonVeto = 0;
    Long64_t afterJetCount = 0;
    Long64_t afterConstituents = 0;
    Long64_t afterDijetMassWindow = 0;
    Long64_t afterRecoilMassWindow = 0;
    Long64_t finalSelected = 0;

    void print(const std::string &processName) const {
        if (!PRINT_CUT_STATISTICS)
            return;
        std::cout << "\n═══════════════════════════════════════════════════" << std::endl;
        std::cout << "Статистика отбора для процесса: " << processName << std::endl;
        std::cout << "═══════════════════════════════════════════════════" << std::endl;
        std::cout << "Всего событий:                    " << totalEvents << std::endl;
        std::cout << "После lepton veto:                " << afterLeptonVeto << " ("
                  << 100.0 * afterLeptonVeto / totalEvents << "%)" << std::endl;
        std::cout << "После требования 2 джетов:        " << afterJetCount << " ("
                  << 100.0 * afterJetCount / totalEvents << "%)" << std::endl;
        std::cout << "После требования конституентов:   " << afterConstituents << " ("
                  << 100.0 * afterConstituents / totalEvents << "%)" << std::endl;
        std::cout << "После photon veto:                " << afterPhotonVeto << " ("
                  << 100.0 * afterPhotonVeto / totalEvents << "%)" << std::endl;
        if (APPLY_DIJET_MASS_WINDOW) {
            std::cout << "После окна массы диджета (65-110 ГэВ): " << afterDijetMassWindow << " ("
                      << 100.0 * afterDijetMassWindow / totalEvents << "%)" << std::endl;
        }
        if (APPLY_RECOIL_MASS_WINDOW) {
            std::cout << "После окна массы отдачи (110-165 ГэВ): " << afterRecoilMassWindow << " ("
                      << 100.0 * afterRecoilMassWindow / totalEvents << "%)" << std::endl;
        }
        std::cout << "───────────────────────────────────────────────────" << std::endl;
        std::cout << "ФИНАЛЬНО ОТОБРАНО:                  " << finalSelected << " ("
                  << 100.0 * finalSelected / totalEvents << "%)" << std::endl;
        std::cout << "═══════════════════════════════════════════════════\n" << std::endl;
    }
};

#endif // ZH_INVISIBLE_ANALYSIS_H
