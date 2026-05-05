#ifndef PFO_CHARGED_NEUTRAL_COUNT_H
#define PFO_CHARGED_NEUTRAL_COUNT_H

#include <set>
#include <string>

// ───────────────────────────────────────────────────────────────
// Имена ветвей дерева
// ───────────────────────────────────────────────────────────────
const std::string TREE_NAME = "outputTree";
const std::string BRANCH_PARTICLE_TYPE = "particleType";
const std::string BRANCH_PFO_ENERGY = "pfoE"; // нужен только для чтения, но не для cuts

// ───────────────────────────────────────────────────────────────
// Параметры гистограмм
// ───────────────────────────────────────────────────────────────
const std::string HIST_TITLE = "Number of PFOs per Event";
const std::string X_AXIS_TITLE = "Number of PFOs";
const std::string Y_AXIS_TITLE = "Events / bin";

const int NUM_BINS = 100;
const int LEFT_BOARD = 0;
const int RIGHT_BOARD = 160;

const int HIST_COLOR_CHARGED = kBlue + 2;
const int HIST_COLOR_NEUTRAL = kRed + 2;

// ───────────────────────────────────────────────────────────────
// Настройки отображения
// ───────────────────────────────────────────────────────────────
const bool DRAW_MEAN_LINES = true;
const bool PRINT_TYPE_STATISTICS = true;
const bool SORT_STATS_BY_COUNT = true;

// ───────────────────────────────────────────────────────────────
// Классификация частиц по PDG-кодам
// ───────────────────────────────────────────────────────────────
// Заряженные частицы (проверяем по абсолютному значению)
// 11 = e±, 13 = μ±, 211 = π±, 321 = K±, 2212 = p/p̄
const std::set<int> CHARGED_PDG_CODES = {11, 13, 211, 321, 2212};

// Нейтральные частицы (проверяем точное значение)
// 22 = γ, 111 = π⁰, 130 = K⁰L, 310 = K⁰S, 2112 = n, 221 = η, 331 = η'
const std::set<int> NEUTRAL_PDG_CODES = {22, 111, 130, 310, 2112, 221, 331};

// Флаг: считать частицы, не попавшие ни в один список, как нейтральные
const bool TREAT_UNKNOWN_AS_NEUTRAL = false;

// ───────────────────────────────────────────────────────────────
// Настройки логирования и вывода
// ───────────────────────────────────────────────────────────────
const int LOG_INTERVAL_EVENTS = 1000;
const bool LOG_PERCENTAGE = true;

// ───────────────────────────────────────────────────────────────
// Базовая директория для результатов
// ───────────────────────────────────────────────────────────────
const std::string OUTPUT_BASE_DIR = "../pdf_results";

#endif // PFO_CHARGED_NEUTRAL_COUNT_H
