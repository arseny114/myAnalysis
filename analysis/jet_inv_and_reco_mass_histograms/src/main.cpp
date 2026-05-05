// Скрипт для построения распределения инвариантной массы и массы отдачи джетов в событии и
// наложения условий отбора.
// Запуск: ./pfo_charged_neutral_count <input_file.root> [-o output_dir]

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "../include/jet_inv_and_reco_mass_histograms.h"

namespace fs = std::filesystem;

// ───────────────────────────────────────────────────────────────
// Вспомогательная функция для логирования прогресса
// ───────────────────────────────────────────────────────────────
void logProgress(Long64_t current, Long64_t total, const std::string &prefix = "") {
    static auto startTime = std::chrono::high_resolution_clock::now();

    if (current % LOG_INTERVAL_EVENTS == 0 || current == total) {
        auto now = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();

        std::cout << "[" << prefix << "] ";
        if (LOG_PERCENTAGE && total > 0) {
            double pct = 100.0 * current / total;
            std::cout << std::fixed;
            std::cout << "Событие " << current << "/" << total << " (" << pct << "%)";
        } else {
            std::cout << "Событие " << current << "/" << total;
        }
        std::cout << " | Прошло: " << elapsed << " с";

        if (current > 0 && total > 0) {
            double rate = static_cast<double>(current) / elapsed;
            double remaining = (total - current) / rate;
            std::cout << " | Осталось: ~" << static_cast<int>(remaining) << " с";
        }
        std::cout << std::endl;
        std::cout.flush();
    }
}

// ───────────────────────────────────────────────────────────────
// Функция расчёта массы отдачи
// ───────────────────────────────────────────────────────────────
double calculateRecoilMass(const TLorentzVector &system, double sqrtS) {
    double recoilE = sqrtS - system.E();
    double recoilP2 =
        std::pow(system.Px(), 2) + std::pow(system.Py(), 2) + std::pow(system.Pz(), 2);
    double recoilMass2 = recoilE * recoilE - recoilP2;
    return (recoilMass2 > 0) ? std::sqrt(recoilMass2) : 0.0;
}

// ───────────────────────────────────────────────────────────────
// Функция проверки veto на изолированные лептоны
// ───────────────────────────────────────────────────────────────
bool hasIsolatedLepton(const std::vector<int> *flags) {
    if (!flags)
        return false;
    for (int flag : *flags) {
        if (flag == 1)
            return true;
    }
    return false;
}

// ───────────────────────────────────────────────────────────────
// Функция суммирования 4-импульсов джетов с cut по энергии
// ───────────────────────────────────────────────────────────────
TLorentzVector sumJets(const std::vector<double> *jetE, const std::vector<double> *jetPx,
                       const std::vector<double> *jetPy, const std::vector<double> *jetPz,
                       double minEnergy) {
    TLorentzVector sum(0., 0., 0., 0.);
    for (size_t i = 0; i < jetE->size(); ++i) {
        if (jetE->at(i) >= minEnergy) {
            sum += TLorentzVector(jetPx->at(i), jetPy->at(i), jetPz->at(i), jetE->at(i));
        }
    }
    return sum;
}

// ───────────────────────────────────────────────────────────────
// Функция отрисовки гистограммы с линиями масс
// ───────────────────────────────────────────────────────────────
void drawHistogramWithMarkers(TH1F *hist, const std::string &canvasTitle, const std::string &xTitle,
                              const std::string &outputFile, double markMass = -1,
                              const std::string &markLabel = "", Color_t markColor = kRed) {
    TCanvas *c = new TCanvas(canvasTitle.c_str(), canvasTitle.c_str(), 800, 600);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("Events");
    hist->Draw("HIST");

    c->cd();
    c->Update();

    if (markMass > 0) {
        double ymin = gPad->GetUymin();
        double ymax = gPad->GetUymax();

        TLine *line = new TLine(markMass, ymin, markMass, ymax);
        line->SetLineColor(markColor);
        line->SetLineWidth(3);
        line->SetLineStyle(7);
        line->Draw();

        TLatex *label = new TLatex(markMass + 5, ymax * 0.85, markLabel.c_str());
        label->SetTextColor(markColor);
        label->SetTextSize(0.035);
        label->Draw();
    }

    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;
    delete c;
}

// ───────────────────────────────────────────────────────────────
// Извлечение имени процесса из пути к файлу
// Ожидает формат: merged_<PROCESS_NAME>.root
// ───────────────────────────────────────────────────────────────
std::string extractProcessName(const std::string &filepath) {
    fs::path p(filepath);
    std::string filename = p.filename().string();

    // Удаляем префикс "merged_"
    const std::string prefix = "merged_";
    if (filename.find(prefix) == 0) {
        filename = filename.substr(prefix.length());
    }

    // Удаляем суффикс ".root"
    const std::string suffix = ".root";
    if (filename.size() >= suffix.size() &&
        filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
        filename = filename.substr(0, filename.size() - suffix.size());
    }

    return filename; // вернёт "E240_qqHX"
}

// ───────────────────────────────────────────────────────────────
// Вывод справки по использованию
// ───────────────────────────────────────────────────────────────
void printUsage(const char *progName) {
    std::cout << "Использование:\n"
              << "  " << progName << " <input_root_file> [options]\n\n"
              << "Обязательные аргументы:\n"
              << "  input_root_file    Путь к входному ROOT-файлу с деревом анализа\n\n"
              << "Опции:\n"
              << "  -h, --help         Показать эту справку\n"
              << "  -o, --output-dir   Базовая директория для результатов "
                 "(по умолчанию: ../pdf_results)\n"
              << "                     Внутри создаётся поддиректория с именем "
                 "процесса\n\n"
              << "Примеры:\n"
              << "  " << progName << " merged_E240_qqHX.root\n"
              << "  " << progName << " /path/to/merged_E240_qqHinvi.root -o ./plots\n";
}

// ───────────────────────────────────────────────────────────────
int main(int argc, char *argv[]) {
    // ─────────────────────────────────────────────────────────────
    // Парсинг аргументов командной строки
    // ─────────────────────────────────────────────────────────────
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string inputRootFile;
    std::string outputBaseDir = OUTPUT_BASE_DIR;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -o/--output-dir требуется указать путь к "
                             "директории"
                          << std::endl;
                return 1;
            }
            outputBaseDir = argv[++i];
        } else if (arg[0] != '-') {
            inputRootFile = arg;
        } else {
            std::cerr << "Неизвестный аргумент: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (inputRootFile.empty()) {
        std::cerr << "Ошибка: не указан входной файл" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    // ─────────────────────────────────────────────────────────────
    // Извлечение имени процесса и подготовка путей
    // ─────────────────────────────────────────────────────────────
    std::string processName = extractProcessName(inputRootFile);
    std::cout << "Процесс: " << processName << std::endl;
    std::cout << "Входной файл: " << inputRootFile << std::endl;

    // Формируем путь к поддиректории процесса и создаём её
    fs::path processOutputDir = fs::path(outputBaseDir) / processName;
    try {
        fs::create_directories(processOutputDir);
        std::cout << "Директория результатов: " << fs::absolute(processOutputDir) << std::endl;
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Ошибка при создании директории " << processOutputDir << ": " << e.what()
                  << std::endl;
        return 1;
    }

    // Формируем имена выходных файлов:
    // базовая_директория/процесс/имя_процесса.pdf
    auto makeOutputPath = [&](const std::string &basename) -> std::string {
        return (processOutputDir / (basename + "_" + processName + ".pdf")).string();
    };

    const std::string OUTPUT_PDF_INV_MASS_EXCL = makeOutputPath("invariant_mass_jets_exclusive");
    const std::string OUTPUT_PDF_RECOIL_MASS_EXCL = makeOutputPath("recoil_mass_jets_exclusive");
    const std::string OUTPUT_PDF_INV_MASS_INCL = makeOutputPath("invariant_mass_jets_inclusive");
    const std::string OUTPUT_PDF_RECOIL_MASS_INCL = makeOutputPath("recoil_mass_jets_inclusive");

    // ─────────────────────────────────────────────────────────────
    // Инициализация ROOT
    // ─────────────────────────────────────────────────────────────
    gStyle->SetOptStat(1);
    auto startTime = std::chrono::high_resolution_clock::now();

    // ─────────────────────────────────────────────────────────────
    // Открытие входного файла
    // ─────────────────────────────────────────────────────────────
    std::cout << "Открытие файла: " << inputRootFile << std::endl;
    TFile *inputFile = TFile::Open(inputRootFile.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Ошибка: не удалось открыть файл " << inputRootFile << std::endl;
        return 1;
    }

    TTree *tree = dynamic_cast<TTree *>(inputFile->Get(TREE_NAME.c_str()));
    if (!tree) {
        std::cerr << "Ошибка: дерево " << TREE_NAME << " не найдено" << std::endl;
        inputFile->Close();
        return 1;
    }
    std::cout << "Дерево найдено. Всего событий: " << tree->GetEntries() << std::endl;

    // ─────────────────────────────────────────────────────────────
    // Установка ветвей для exclusive джетов
    // ─────────────────────────────────────────────────────────────
    std::vector<double> *exclJetE = nullptr, *exclJetPx = nullptr;
    std::vector<double> *exclJetPy = nullptr, *exclJetPz = nullptr;
    tree->SetBranchAddress("exclusiveJetE", &exclJetE);
    tree->SetBranchAddress("exclusiveJetPx", &exclJetPx);
    tree->SetBranchAddress("exclusiveJetPy", &exclJetPy);
    tree->SetBranchAddress("exclusiveJetPz", &exclJetPz);

    // ─────────────────────────────────────────────────────────────
    // Установка ветвей для inclusive джетов
    // ─────────────────────────────────────────────────────────────
    std::vector<double> *inclJetE = nullptr, *inclJetPx = nullptr;
    std::vector<double> *inclJetPy = nullptr, *inclJetPz = nullptr;
    tree->SetBranchAddress("inclusiveJetE", &inclJetE);
    tree->SetBranchAddress("inclusiveJetPx", &inclJetPx);
    tree->SetBranchAddress("inclusiveJetPy", &inclJetPy);
    tree->SetBranchAddress("inclusiveJetPz", &inclJetPz);

    // ─────────────────────────────────────────────────────────────
    // Ветви для lepton veto
    // ─────────────────────────────────────────────────────────────
    std::vector<int> *isIsolatedLeptonFlag = nullptr;
    if (APPLY_LEPTON_VETO) {
        tree->SetBranchAddress("isIsolatedLeptonFlag", &isIsolatedLeptonFlag);
    }

    // ─────────────────────────────────────────────────────────────
    // Создание гистограмм (4 штуки)
    // ─────────────────────────────────────────────────────────────
    TH1F *hInvMassExcl =
        new TH1F("hInvMassExcl", "Invariant Mass (Exclusive Jets);Invariant Mass [GeV];Events",
                 NUM_BINS, MASS_MIN, MASS_MAX);
    TH1F *hRecoilMassExcl =
        new TH1F("hRecoilMassExcl", "Recoil Mass (Exclusive Jets);Recoil Mass [GeV];Events",
                 NUM_BINS, MASS_MIN, MASS_MAX);

    TH1F *hInvMassIncl =
        new TH1F("hInvMassIncl", "Invariant Mass (Inclusive Jets);Invariant Mass [GeV];Events",
                 NUM_BINS, MASS_MIN, MASS_MAX);
    TH1F *hRecoilMassIncl =
        new TH1F("hRecoilMassIncl", "Recoil Mass (Inclusive Jets);Recoil Mass [GeV];Events",
                 NUM_BINS, MASS_MIN, MASS_MAX);

    // ─────────────────────────────────────────────────────────────
    // Основной цикл по событиям
    // ─────────────────────────────────────────────────────────────
    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessedExcl = 0, nProcessedIncl = 0;

    std::cout << "\nНачало обработки событий..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        logProgress(i + 1, nEntries, "Processing");

        // ── Common cut: veto на изолированные лептоны ─────────────
        if (APPLY_LEPTON_VETO && hasIsolatedLepton(isIsolatedLeptonFlag)) {
            continue;
        }

        // ═════════════════════════════════════════════════════════
        // EXCLUSIVE JETS обработка
        // ═════════════════════════════════════════════════════════
        if (exclJetE->size() >= 2) {
            TLorentzVector summedExcl =
                sumJets(exclJetE, exclJetPx, exclJetPy, exclJetPz, MIN_JET_ENERGY_EXCL);
            double invMassExcl = summedExcl.M();
            double recoilMassExcl = calculateRecoilMass(summedExcl, CENTER_OF_MASS_ENERGY);

            bool passZWindow = !APPLY_Z_MASS_WINDOW || (invMassExcl >= MZ_WINDOW_MIN_EXCL &&
                                                        invMassExcl <= MZ_WINDOW_MAX_EXCL);

            if (passZWindow) {
                hInvMassExcl->Fill(invMassExcl);
                hRecoilMassExcl->Fill(recoilMassExcl);
                nProcessedExcl++;
            }
        }

        // ═════════════════════════════════════════════════════════
        // INCLUSIVE JETS обработка
        // ═════════════════════════════════════════════════════════
        if (inclJetE->size() >= 1) {
            TLorentzVector summedIncl =
                sumJets(inclJetE, inclJetPx, inclJetPy, inclJetPz, MIN_JET_ENERGY_INCL);
            double invMassIncl = summedIncl.M();
            double recoilMassIncl = calculateRecoilMass(summedIncl, CENTER_OF_MASS_ENERGY);

            bool passZWindow = !APPLY_Z_MASS_WINDOW || (invMassIncl >= MZ_WINDOW_MIN_INCL &&
                                                        invMassIncl <= MZ_WINDOW_MAX_INCL);

            if (passZWindow) {
                hInvMassIncl->Fill(invMassIncl);
                hRecoilMassIncl->Fill(recoilMassIncl);
                nProcessedIncl++;
            }
        }
    }

    // ─────────────────────────────────────────────────────────────
    // Итоговая статистика
    // ─────────────────────────────────────────────────────────────
    auto endTime = std::chrono::high_resolution_clock::now();
    auto totalSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    std::cout << "\n═══════════════════════════════════════════════════" << std::endl;
    std::cout << "Обработка завершена!" << std::endl;
    std::cout << "Всего событий: " << nEntries << std::endl;
    std::cout << "Прошло времени: " << totalSec << " с (" << totalSec / 60.0 << " мин)"
              << std::endl;
    std::cout << "───────────────────────────────────────────────────" << std::endl;
    std::cout << "EXCLUSIVE jets:" << std::endl;
    std::cout << "  Обработано событий: " << nProcessedExcl << std::endl;
    std::cout << "  hInvMassExcl entries: " << hInvMassExcl->GetEntries() << std::endl;
    std::cout << "  hRecoilMassExcl entries: " << hRecoilMassExcl->GetEntries() << std::endl;
    std::cout << "───────────────────────────────────────────────────" << std::endl;
    std::cout << "INCLUSIVE jets:" << std::endl;
    std::cout << "  Обработано событий: " << nProcessedIncl << std::endl;
    std::cout << "  hInvMassIncl entries: " << hInvMassIncl->GetEntries() << std::endl;
    std::cout << "  hRecoilMassIncl entries: " << hRecoilMassIncl->GetEntries() << std::endl;
    std::cout << "═══════════════════════════════════════════════════" << std::endl;

    // ─────────────────────────────────────────────────────────────
    // Отрисовка и сохранение гистограмм
    // ─────────────────────────────────────────────────────────────
    std::cout << "\nОтрисовка гистограмм..." << std::endl;

    drawHistogramWithMarkers(hInvMassExcl, "cInvExcl", "Invariant Mass [GeV]",
                             OUTPUT_PDF_INV_MASS_EXCL, MZ, "M_{Z} = 91.2 GeV", kRed);
    drawHistogramWithMarkers(hRecoilMassExcl, "cRecoilExcl", "Recoil Mass [GeV]",
                             OUTPUT_PDF_RECOIL_MASS_EXCL, MH, "M_{H} = 125 GeV", kRed);
    drawHistogramWithMarkers(hInvMassIncl, "cInvIncl", "Invariant Mass [GeV]",
                             OUTPUT_PDF_INV_MASS_INCL, MZ, "M_{Z} = 91.2 GeV", kBlue);
    drawHistogramWithMarkers(hRecoilMassIncl, "cRecoilIncl", "Recoil Mass [GeV]",
                             OUTPUT_PDF_RECOIL_MASS_INCL, MH, "M_{H} = 125 GeV", kBlue);

    // ─────────────────────────────────────────────────────────────
    // Очистка памяти
    // ─────────────────────────────────────────────────────────────
    delete hInvMassExcl;
    delete hRecoilMassExcl;
    delete hInvMassIncl;
    delete hRecoilMassIncl;
    inputFile->Close();
    delete inputFile;

    std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
              << std::endl;
    return 0;
}
