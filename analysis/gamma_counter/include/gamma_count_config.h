// include/gamma_count_config.h
#ifndef GAMMA_COUNT_CONFIG_H
#define GAMMA_COUNT_CONFIG_H

#include <string>

// ───────────────────────────────────────────────────────────────
// Параметры дерева и ветвей
// ───────────────────────────────────────────────────────────────
const std::string TREE_NAME = "outputTree";

// Ветви с 4-импульсами всех PFO
const std::string BRANCH_PFO_E = "pfoE";
const std::string BRANCH_PFO_PX = "pfoPx";
const std::string BRANCH_PFO_PY = "pfoPy";
const std::string BRANCH_PFO_PZ = "pfoPz";

// Ветвь с типом частицы (PDG-код)
const std::string BRANCH_PARTICLE_TYPE = "particleType";

// PDG-код фотона
const int PDG_PHOTON = 22;

// ───────────────────────────────────────────────────────────────
// Энергетические пороги для подсчёта фотонов
// ───────────────────────────────────────────────────────────────
const double THRESHOLD_LOW = 30.0;  // ГэВ
const double THRESHOLD_HIGH = 50.0; // ГэВ

// ───────────────────────────────────────────────────────────────
// Параметры гистограммы спектра
// ───────────────────────────────────────────────────────────────
const double SPECTRUM_MIN = 0.0;   // ГэВ
const double SPECTRUM_MAX = 250.0; // ГэВ
const int SPECTRUM_BINS = 250;     // 1 ГэВ на бин

// ───────────────────────────────────────────────────────────────
// Директория для результатов
// ───────────────────────────────────────────────────────────────
const std::string OUTPUT_BASE_DIR = "../pdf_results";

#endif // GAMMA_COUNT_CONFIG_H
