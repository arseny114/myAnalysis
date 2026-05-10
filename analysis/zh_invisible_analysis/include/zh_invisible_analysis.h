// include/zh_invisible_analysis.h
#ifndef ZH_INVISIBLE_ANALYSIS_H
#define ZH_INVISIBLE_ANALYSIS_H

#include <string>

// =============================================================================
// НАСТРОЙКИ ОТБОРА
// =============================================================================

// Вето на изолированные лептоны
#define APPLY_LEPTON_VETO true
const double LEPTON_ISO_COS_CONE_ANGLE = 0.985; // cosConeAngle
const double LEPTON_ISO_MIN_TRACK_E_GEV = 5.0;  // isoMinTrackEnergy
const double LEPTON_ISO_MAX_TRACK_E_GEV = 1e20; // isoMaxTrackEnergy
const double LEPTON_ISO_MIN_CONE_E_GEV = 0.0;   // isoMinConeEnergy
const double LEPTON_ISO_MAX_CONE_E_GEV = 2.0;   // isoMaxConeEnergy

// Вето на высокоэнергетические фотоны
#define APPLY_HIGH_E_PHOTON_VETO true
#define PHOTON_ENERGY_CUT_GEV 30.0 // Порог энергии фотона для вето

// Вето на изолированные фотоны
#define APPLY_ISOLATED_PHOTON_VETO true
const double PHOTON_ISO_MIN_ENERGY_GEV = 5.0;
const double PHOTON_ISO_MAX_CONE_ENERGY_GEV = 2.0;
const double PHOTON_ISO_COS_CONE_ANGLE = 0.985;

// Требование ровно 2 джета в инклюзивном режиме
#define REQUIRE_EXACTLY_TWO_INCLUSIVE_JETS true

// Требование минимального числа конституентов в каждом из двух джетов
#define REQUIRE_MIN_CONSTITUENTS_PER_JET true
#define MIN_CONSTITUENTS_PER_JET 6

// Кат на MET (Missing Transverse Energy)
#define APPLY_MET_CUT true
const double MET_CUT_MIN_GEV = 20.0;

// Окно инвариантной массы диджета (ГэВ)
#define APPLY_DIJET_MASS_WINDOW false
#define DIJET_MASS_WINDOW_MIN_GEV 75.0
#define DIJET_MASS_WINDOW_MAX_GEV 105.0

// Окно массы отдачи (ГэВ)
#define APPLY_RECOIL_MASS_WINDOW false
#define RECOIL_MASS_WINDOW_MIN_GEV 105.0
#define RECOIL_MASS_WINDOW_MAX_GEV 145.0

// Кат на полярный угол системы двух джетов (Z-бозона)
#define APPLY_COS_THETA_Z_CUT true
const double COS_THETA_Z_CUT = 0.98;

// Настройки эллиптического ката на M_jj vs M_recoil
#define APPLY_ELLIPSE_CUT true

// Параметры эллипса: ((x-x0)cosθ + (y-y0)sinθ)²/a² + (-(x-x0)sinθ + (y-y0)cosθ)²/b² <= 1
const double ELLIPSE_CX_GEV = 87.5;           // Центр по M_inv
const double ELLIPSE_CY_GEV = 130.0;          // Центр по M_recoil
const double ELLIPSE_A_GEV = 27.95;           // Большая полуось
const double ELLIPSE_B_GEV = 5.59;            // Малая полуось
const double ELLIPSE_THETA = std::atan(-2.0); // Угол поворота

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

const int COS_THETA_JET_BINS = 100;
const double COS_THETA_JET_MIN = -1.0;
const double COS_THETA_JET_MAX = 1.0;

// Параметры для гистограмм потерянной поперечной энергии
const int MET_PFO_BINS = 100;
const double MET_PFO_MIN = 0.0;
const double MET_PFO_MAX = 150.0;

const int MET_JET_BINS = 100;
const double MET_JET_MIN = 0.0;
const double MET_JET_MAX = 150.0;

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
    Long64_t afterHighEPhotonVeto = 0;
    Long64_t afterIsoPhotonVeto = 0;
    Long64_t afterJetCount = 0;
    Long64_t afterConstituents = 0;
    Long64_t afterMetCut = 0;
    Long64_t afterDijetMassWindow = 0;
    Long64_t afterCosThetaZCut = 0;
    Long64_t afterRecoilMassWindow = 0;
    Long64_t afterEllipseCut = 0;
    Long64_t finalSelected = 0;

    void print(const std::string &processName) const {
        if (!PRINT_CUT_STATISTICS)
            return;
        std::cout << "\n═══════════════════════════════════════════════════" << std::endl;
        std::cout << "Статистика отбора для процесса: " << processName << std::endl;
        std::cout << "═══════════════════════════════════════════════════" << std::endl;
        std::cout << "Всего событий:                    " << totalEvents << std::endl;
        if (APPLY_LEPTON_VETO) {
            std::cout << "После lepton veto:                " << afterLeptonVeto << " ("
                      << 100.0 * afterLeptonVeto / totalEvents << "%)" << std::endl;
        }
        if (REQUIRE_EXACTLY_TWO_INCLUSIVE_JETS) {
            std::cout << "После требования 2 джетов:        " << afterJetCount << " ("
                      << 100.0 * afterJetCount / totalEvents << "%)" << std::endl;
        }
        if (REQUIRE_MIN_CONSTITUENTS_PER_JET) {
            std::cout << "После требования конституентов:   " << afterConstituents << " ("
                      << 100.0 * afterConstituents / totalEvents << "%)" << std::endl;
        }
        if (APPLY_MET_CUT) {
            std::cout << "После MET:                        " << afterMetCut << " ("
                      << 100.0 * afterMetCut / totalEvents << "%)" << std::endl;
        }
        if (APPLY_HIGH_E_PHOTON_VETO) {
            std::cout << "После high-E photon veto:         " << afterHighEPhotonVeto << " ("
                      << 100.0 * afterHighEPhotonVeto / totalEvents << "%)" << std::endl;
        }
        if (APPLY_ISOLATED_PHOTON_VETO) {
            std::cout << "После isolated photon veto:       " << afterIsoPhotonVeto << " ("
                      << 100.0 * afterIsoPhotonVeto / totalEvents << "%)" << std::endl;
        }
        if (APPLY_DIJET_MASS_WINDOW) {
            std::cout << "После окна массы диджета: " << afterDijetMassWindow << " ("
                      << 100.0 * afterDijetMassWindow / totalEvents << "%)" << std::endl;
        }
        if (APPLY_COS_THETA_Z_CUT) {
            std::cout << "После |cos#theta_{Z}| < " << COS_THETA_Z_CUT << ":  " << afterCosThetaZCut
                      << " (" << 100.0 * afterCosThetaZCut / totalEvents << "%)" << std::endl;
        }
        if (APPLY_RECOIL_MASS_WINDOW) {
            std::cout << "После окна массы отдачи: " << afterRecoilMassWindow << " ("
                      << 100.0 * afterRecoilMassWindow / totalEvents << "%)" << std::endl;
        }
        if (APPLY_ELLIPSE_CUT) {
            std::cout << "После эллиптического cut'а:      " << afterEllipseCut << " ("
                      << 100.0 * afterEllipseCut / totalEvents << "%)" << std::endl;
        }
        std::cout << "───────────────────────────────────────────────────" << std::endl;
        std::cout << "ФИНАЛЬНО ОТОБРАНО:                  " << finalSelected << " ("
                  << 100.0 * finalSelected / totalEvents << "%)" << std::endl;
        std::cout << "═══════════════════════════════════════════════════\n" << std::endl;
    }
};

#endif // ZH_INVISIBLE_ANALYSIS_H
