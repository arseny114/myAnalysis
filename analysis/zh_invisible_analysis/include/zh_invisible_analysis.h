// include/zh_invisible_analysis.h
#ifndef ZH_INVISIBLE_ANALYSIS_H
#define ZH_INVISIBLE_ANALYSIS_H

#include <iomanip>
#include <iostream>
#include <string>

// =============================================================================
// НАСТРОЙКИ ОТБОРА
// =============================================================================

// =============================================================================
// ПРЕДОТБОРЫ
// =============================================================================

// --- 1. Veto на изолированные лептоны ---
#define APPLY_PRE_LEPTON_VETO true
const double LEPTON_ISO_COS_CONE_ANGLE = 0.985; // cosConeAngle для изоляции
const double LEPTON_ISO_MIN_TRACK_E_GEV = 5.0;  // Мин. энергия трека-кандидата
const double LEPTON_ISO_MAX_TRACK_E_GEV = 1e20; // Макс. энергия трека
const double LEPTON_ISO_MIN_CONE_E_GEV = 0.0;   // Мин. энергия в конусе
const double LEPTON_ISO_MAX_CONE_E_GEV = 2.0;   // Макс. энергия в конусе (порог изоляции)

// --- 2. Veto на высокоэнергетические фотоны ---
#define APPLY_PRE_HIGH_E_PHOTON_VETO true
#define PHOTON_ENERGY_CUT_GEV 30.0 // Порог энергии фотона для вето

// --- 3. Veto на изолированные фотоны ---
#define APPLY_PRE_ISOLATED_PHOTON_VETO true
const double PHOTON_ISO_MIN_ENERGY_GEV = 5.0;      // Мин. энергия фотона-кандидата
const double PHOTON_ISO_MAX_CONE_ENERGY_GEV = 2.0; // Макс. энергия в конусе вокруг фотона
const double PHOTON_ISO_COS_CONE_ANGLE = 0.985;    // cosConeAngle для фотонной изоляции

// --- 4. Требование ровно 2 инклюзивных джета ---
#define APPLY_PRE_TWO_JETS_REQUIREMENT true

// --- 5. Требование минимального числа конституентов в каждом джете ---
#define APPLY_PRE_CONSTITUENTS_REQUIREMENT true
#define MIN_CONSTITUENTS_PER_JET 6 // Мин. число PFO в каждом из двух джетов

// =============================================================================
// ОСНОВНЫЕ ОТБОРЫ
// =============================================================================

// --- 6. Кат на MET (Missing Transverse Energy) ---
#define APPLY_MAIN_MET_CUT true
const double MET_CUT_MIN_GEV = 30.0;
const double MET_CUT_MAX_GEV = 60.0;

// --- 7. Кат на deltaPhi между джетами ---
#define APPLY_MAIN_DELTA_PHI_CUT true
const double DELTA_PHI_CUT_MAX = 3.0;

// --- 8. Кат на полярный угол системы двух джетов ---
#define APPLY_MAIN_COS_THETA_Z_CUT true
const double COS_THETA_Z_CUT = 0.7;

// --- 9. Окно инвариантной массы диджета (ГэВ) ---
#define APPLY_MAIN_DIJET_MASS_WINDOW true
#define DIJET_MASS_WINDOW_MIN_GEV 78.0
#define DIJET_MASS_WINDOW_MAX_GEV 100.0

// --- 10. Кат на Pmiss ---
#define APPLY_MAIN_PMISS_CUT true
const double PMISS_CUT_MIN_GEV = 30.0;
const double PMISS_CUT_MAX_GEV = 70.0;

// --- 11. Окно массы отдачи (ГэВ) ---
#define APPLY_MAIN_RECOIL_MASS_WINDOW false
#define RECOIL_MASS_WINDOW_MIN_GEV 100.0
#define RECOIL_MASS_WINDOW_MAX_GEV 165.0

// --- 12. Эллиптический кат на плоскости M_jj vs M_recoil ---
#define APPLY_MAIN_ELLIPSE_CUT false
// Параметры эллипса: ((x-x₀)cosθ + (y-y₀)sinθ)²/a² + (-(x-x₀)sinθ + (y-y₀)cosθ)²/b² ≤ 1
const double ELLIPSE_CX_GEV = 85.0;  // Центр по M_inv (ГэВ)
const double ELLIPSE_CY_GEV = 132.5; // Центр по M_recoil (ГэВ)
const double ELLIPSE_A_GEV = 25.00;  // Большая полуось (ГэВ)
const double ELLIPSE_B_GEV = 7.0;    // Малая полуось (ГэВ)
const double ELLIPSE_THETA = -55.0;  // Угол поворота (градусы)

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

const int RECOIL_STACK_BINS = 70;
const double RECOIL_STACK_MIN_GEV = RECOIL_MASS_WINDOW_MIN_GEV;
const double RECOIL_STACK_MAX_GEV = RECOIL_MASS_WINDOW_MAX_GEV;
const double RECOIL_STACK_MIN_Y = 1e-2;
const double RECOIL_STACK_MAX_Y = 30000;

// Параметры для гистограмм угловых распределений
const int COS_THETA_Z_BINS = 100;
const double COS_THETA_Z_MIN = -1.0;
const double COS_THETA_Z_MAX = 1.0;

const int DELTA_R_BINS = 100;
const double DELTA_R_MIN = 0.0;
const double DELTA_R_MAX = 5.0;

const int COS_THETA_JET_BINS = 100;
const double COS_THETA_JET_MIN = -1.0;
const double COS_THETA_JET_MAX = 1.0;

// Параметры для угловых распределений разниц между джетами
const int DELTA_THETA_BINS = 100;
const double DELTA_THETA_MIN = 0.0;
const double DELTA_THETA_MAX = M_PI;

const int DELTA_PHI_BINS = 100;
const double DELTA_PHI_MIN = 0.0;
const double DELTA_PHI_MAX = M_PI;

// Параметры для гистограмм потерянной поперечной энергии
const int MET_PFO_BINS = 100;
const double MET_PFO_MIN = 0.0;
const double MET_PFO_MAX = 150.0;

const int MET_JET_BINS = 100;
const double MET_JET_MIN = 0.0;
const double MET_JET_MAX = 150.0;

// Параметры для гистограмм недостающего 3-импульса (Pmiss)
const int PMISS_BINS = 100;
const double PMISS_MIN_GEV = 0.0;
const double PMISS_MAX_GEV = 150.0;

const int COS_THETA_PMISS_BINS = 100;
const double COS_THETA_PMISS_MIN = -1.0;
const double COS_THETA_PMISS_MAX = 1.0;

// Параметры для гистограммы энергии диджета
const int DIJET_ENERGY_BINS = 200;
const double DIJET_ENERGY_MIN_GEV = 0.0;
const double DIJET_ENERGY_MAX_GEV = 300.0;

// =============================================================================
// НАСТРОЙКИ ВЫВОДА И ЛОГИРОВАНИЯ
// =============================================================================

const std::string OUTPUT_BASE_DIR = "../pdf_results";

const int LOG_INTERVAL_EVENTS = 100000;
const bool LOG_PERCENTAGE = true;
const bool PRINT_CUT_STATISTICS = true;

// =============================================================================
// СТРУКТУРА ДЛЯ СТАТИСТИКИ ПРОХОЖДЕНИЯ КАТОВ
// =============================================================================

// Считает число символов (не байт) в UTF-8 строке
inline size_t utf8_char_count(const std::string &s) {
    size_t count = 0;
    for (unsigned char c : s) {
        if ((c & 0xC0) != 0x80)
            ++count;
    }
    return count;
}

// Дополняет строку пробелами справа до нужной ширины в символах
inline std::string utf8_ljust(const std::string &s, size_t width) {
    size_t len = utf8_char_count(s);
    if (len >= width)
        return s;
    return s + std::string(width - len, ' ');
}

struct CutStatistics {
    Long64_t totalEvents = 0;

    // Статистика по предотборам
    Long64_t afterPreLeptonVeto = 0;
    Long64_t afterPreHighEPhotonVeto = 0;
    Long64_t afterPreIsoPhotonVeto = 0;
    Long64_t afterPreJetCount = 0;
    Long64_t afterPreConstituents = 0;

    // Основные отборы
    Long64_t afterMetCut = 0;
    Long64_t afterDeltaPhiCut = 0;
    Long64_t afterCosThetaZCut = 0;
    Long64_t afterDijetMassWindow = 0;
    Long64_t afterPmissCut = 0;
    Long64_t afterRecoilMassWindow = 0;
    Long64_t afterEllipseCut = 0;
    Long64_t finalSelected = 0;

    void print(const std::string &processName) const {
        if (!PRINT_CUT_STATISTICS)
            return;

        // Лямбда для вывода строки с корректным выравниванием UTF-8
        auto printRow = [this](const std::string &name, Long64_t passed, Long64_t prev) {
            double pct_prev = (prev > 0) ? 100.0 * passed / prev : 0.0;
            double pct_total = (totalEvents > 0) ? 100.0 * passed / totalEvents : 0.0;
            std::cout << utf8_ljust(name, 55) << std::right << std::setw(8) << passed << "  ("
                      << std::fixed << std::setprecision(2) << std::setw(6) << pct_prev << "% prev)"
                      << "  (" << std::fixed << std::setprecision(2) << std::setw(6) << pct_total
                      << "% tot)\n";
        };

        const std::string sepWide(100, '=');
        const std::string sepNarrow(100, '-');

        std::cout << "\n" << sepWide << "\n";
        std::cout << "Статистика отбора для процесса: " << processName << "\n";
        std::cout << sepWide << "\n\n";

        // ───────────────── ПРЕДОТБОРЫ ─────────────────
        std::cout << "\n" << sepNarrow << "\n";
        std::cout << "ПРЕДОТБОРЫ:\n" << sepNarrow << "\n";
        std::cout << utf8_ljust("  Всего событий:", 55) << std::right << std::setw(8) << totalEvents
                  << "\n";

        Long64_t current = totalEvents;

        if (APPLY_PRE_LEPTON_VETO) {
            printRow("Veto на изолированный лептон:", afterPreLeptonVeto, current);
            current = afterPreLeptonVeto;
        }

        if (APPLY_PRE_HIGH_E_PHOTON_VETO) {
            printRow("Veto на энергичный фотон (>30 GeV):", afterPreHighEPhotonVeto, current);
            current = afterPreHighEPhotonVeto;
        }

        if (APPLY_PRE_ISOLATED_PHOTON_VETO) {
            printRow("Veto на изолированный фотон:", afterPreIsoPhotonVeto, current);
            current = afterPreIsoPhotonVeto;
        }

        if (APPLY_PRE_TWO_JETS_REQUIREMENT) {
            printRow("Требование ровно 2 джетов:", afterPreJetCount, current);
            current = afterPreJetCount;
        }

        if (APPLY_PRE_CONSTITUENTS_REQUIREMENT) {
            printRow("Требование >= 6 конституентов/джет:", afterPreConstituents, current);
            current = afterPreConstituents;
        }

        // ───────────────── ОСНОВНЫЕ ОТБОРЫ ─────────────────
        std::cout << "\n" << sepNarrow << "\n";
        std::cout << "ОСНОВНЫЕ ОТБОРЫ:\n" << sepNarrow << "\n";

        if (APPLY_MAIN_MET_CUT) {
            printRow("MET cut:", afterMetCut, current);
            current = afterMetCut;
        }

        if (APPLY_MAIN_DELTA_PHI_CUT) {
            printRow("Cut on #Delta#phi: ", afterDeltaPhiCut, current);
            current = afterDeltaPhiCut;
        }

        if (APPLY_MAIN_COS_THETA_Z_CUT) {
            printRow("|cos#theta_{Z}:", afterCosThetaZCut, current);
            current = afterCosThetaZCut;
        }

        if (APPLY_MAIN_DIJET_MASS_WINDOW) {
            printRow("Окно массы диджета:", afterDijetMassWindow, current);
            current = afterDijetMassWindow;
        }

        if (APPLY_MAIN_PMISS_CUT) {
            printRow("PMISS cut:", afterPmissCut, current);
            current = afterPmissCut;
        }

        if (APPLY_MAIN_RECOIL_MASS_WINDOW) {
            printRow("Окно массы отдачи:", afterRecoilMassWindow, current);
            current = afterRecoilMassWindow;
        }

        if (APPLY_MAIN_ELLIPSE_CUT) {
            printRow("Эллиптический cut (Mjj vs Mrecoil):", afterEllipseCut, current);
            current = afterEllipseCut;
        }

        // ───────────────── ИТОГ ─────────────────
        double finalPct = (totalEvents > 0) ? 100.0 * current / totalEvents : 0.0;
        std::cout << "\n" << sepNarrow << "\n";
        std::cout << utf8_ljust("ФИНАЛЬНО ОТОБРАНО:", 55) << std::right << std::setw(8) << current
                  << "  (" << std::fixed << std::setprecision(2) << std::setw(6) << finalPct
                  << "%)\n";
        std::cout << sepWide << "\n\n";
    }
};

// =============================================================================
// СТРУКТУРА ДЛЯ СТАТИСТИКИ ИЗОЛИРОВАННЫХ ЭЛЕКТРОНОВ
// =============================================================================
struct IsoElectronStats {
    Long64_t total = 0;
    Long64_t barrel = 0; // |cos(theta)| < 0.7
    Long64_t endcap = 0; // |cos(theta)| >= 0.7

    void print() const {
        if (total == 0)
            return;
        std::cout << "\n── Статистика изолированных электронов ──" << std::endl;
        std::cout << utf8_ljust("  Всего изолированных e:", 40) << total << std::endl;
        std::cout << utf8_ljust("  В барреле (|cos#theta| < 0.7):", 40) << barrel << " ("
                  << std::fixed << std::setprecision(2) << 100.0 * barrel / total << "%)"
                  << std::endl;
        std::cout << utf8_ljust("  В эндкапе (|cos#theta| >= 0.7):", 40) << endcap << " ("
                  << std::fixed << std::setprecision(2) << 100.0 * endcap / total << "%)"
                  << std::endl;
        std::cout << "──────────────────────────────────────────" << std::endl;
    }
};

// =============================================================================
// НАСТРОЙКИ ДЛЯ МНОГОФАЙЛОВОГО АНАЛИЗА И СТЕКА
// =============================================================================

const double LUMINOSITY_FB1 = 5050.0;

// Структура для описания процесса
struct ProcessInfo {
    std::string legendName;
    double weight = 1.0;
    Color_t color = kBlack;
    int fillStyle = 1001;

    // Конструктор для удобной инициализации
    ProcessInfo(std::string leg = "", double w = 1.0, Color_t c = kBlack, int fs = 1001)
        : legendName(std::move(leg)), weight(w), color(c), fillStyle(fs) {}
};

// База весов и стилей по имени файла
std::map<std::string, ProcessInfo> getProcessDatabase() {
    std::map<std::string, ProcessInfo> db;

    db["merged_E240_qqHX.root"] = ProcessInfo{"qqHX", 0.8111, kRed - 3, 1001};
    db["merged_E240_qq.root"] = ProcessInfo{"qq", 513.7456, kGray + 2, 1001};
    db["merged_E240_4f_sw_sl0qq.root"] = ProcessInfo{"4f_sw_sl0qq", 32.8556, kAzure - 4, 1001};
    db["merged_E240_4f_sze_sl0uu.root"] = ProcessInfo{"4f_sze_sl0uu", 4.1155, kGreen + 2, 1001};
    db["merged_E240_4f_sznu_sl0nu_up.root"] =
        ProcessInfo{"4f_sznu_sl0nu_up", 1.8633, kOrange - 3, 1001};
    db["merged_E240_4f_ww_h0cuxx.root"] = ProcessInfo{"4f_ww_h0cuxx", 47.6774, kMagenta - 2, 1001};
    db["merged_E240_4f_zz_h0dtdt.root"] = ProcessInfo{"4f_zz_h0dtdt", 6.4182, kCyan + 1, 1001};
    db["merged_E240_qqHinvi.root"] = ProcessInfo{"qqHinvi (signal)", 0.0080, kRed + 1, 3005};

    return db;
}

// Фиксированный порядок процессов в стеке: от фона к сигналу
const std::vector<std::string> RECOIL_STACK_ORDER = {
    "E240_4f_sznu_sl0nu_up", // nunuZ
    "E240_4f_sw_sl0qq",      // sw
    "E240_4f_ww_h0cuxx",     // ww
    "E240_4f_sze_sl0uu",     // sze
    "E240_4f_zz_h0dtdt",     // zz
    "E240_qq",               // qq background
    "E240_qqHX",             // qqHX
    "E240_qqHinvi"           // signal, всегда последним
};

#endif // ZH_INVISIBLE_ANALYSIS_H
