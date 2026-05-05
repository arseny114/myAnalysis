// include/jet_inv_and_reco_mass_histograms.h
#ifndef JET_MASS_HISTOGRAMS_H
#define JET_MASS_HISTOGRAMS_H

#include <string>

// ───────────────────────────────────────────────────────────────
// Параметры гистограмм
// ───────────────────────────────────────────────────────────────
const std::string TREE_NAME = "outputTree";

const int NUM_BINS = 200;
const double MASS_MIN = 0.0;
const double MASS_MAX = 250.0;

const double MZ = 91.2;
const double MH = 125.26;
const double CENTER_OF_MASS_ENERGY = 240.0;

// ───────────────────────────────────────────────────────────────
// Настройки анализа
// ───────────────────────────────────────────────────────────────
const bool APPLY_LEPTON_VETO = true;
const bool APPLY_Z_MASS_WINDOW = false;

// Cuts для exclusive джетов
const double MIN_JET_ENERGY_EXCL = 5.0;
const double MZ_WINDOW_MIN_EXCL = 70.0;
const double MZ_WINDOW_MAX_EXCL = 110.0;

// Cuts для inclusive джетов
const double MIN_JET_ENERGY_INCL = 2.0;
const double MZ_WINDOW_MIN_INCL = 70.0;
const double MZ_WINDOW_MAX_INCL = 110.0;

// ───────────────────────────────────────────────────────────────
// Настройки логирования
// ───────────────────────────────────────────────────────────────
const int LOG_INTERVAL_EVENTS = 1000;
const bool LOG_PERCENTAGE = true;

// ───────────────────────────────────────────────────────────────
// Директория для результатов
// ───────────────────────────────────────────────────────────────
const std::string OUTPUT_BASE_DIR = "../pdf_results";

#endif // JET_MASS_HISTOGRAMS_H
