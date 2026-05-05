// Скрипт для построения распределения числа заряженных и нейтральных PFO в событии.
// Запуск: ./pfo_charged_neutral_count <input_file.root> [-o output_dir]

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../include/pfo_charged_neutral_count.h"

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

        if (current > 0 && total > 0 && elapsed > 0) {
            double rate = static_cast<double>(current) / elapsed;
            double remaining = (total - current) / rate;
            std::cout << " | Осталось: ~" << static_cast<int>(remaining) << " с";
        }
        std::cout << std::endl;
        std::cout.flush();
    }
}

// ───────────────────────────────────────────────────────────────
// Классификация частиц
// ───────────────────────────────────────────────────────────────
bool isCharged(int pdgCode) {
    int absCode = std::abs(pdgCode);
    return CHARGED_PDG_CODES.find(absCode) != CHARGED_PDG_CODES.end();
}

bool isNeutral(int pdgCode) {
    if (NEUTRAL_PDG_CODES.find(pdgCode) != NEUTRAL_PDG_CODES.end())
        return true;
    if (TREAT_UNKNOWN_AS_NEUTRAL && !isCharged(pdgCode))
        return true;
    return false;
}

// ───────────────────────────────────────────────────────────────
// Вывод статистики по типам частиц
// ───────────────────────────────────────────────────────────────
void printTypeStatistics(const std::map<int, int> &typeCounter, Long64_t totalParticles) {
    if (typeCounter.empty()) {
        std::cout << "\n[Статистика типов] Нет данных" << std::endl;
        return;
    }

    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << "Встреченные типы частиц (всего: " << totalParticles << "):" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    std::vector<std::pair<int, int>> sortedTypes(typeCounter.begin(), typeCounter.end());
    if (SORT_STATS_BY_COUNT) {
        std::sort(sortedTypes.begin(), sortedTypes.end(),
                  [](const auto &a, const auto &b) { return a.second > b.second; });
    }

    std::map<int, std::string> typeLabels = {{22, "γ"},   {111, "π⁰"},  {130, "K⁰L"}, {310, "K⁰S"},
                                             {2112, "n"}, {11, "e±"},   {13, "μ±"},   {211, "π±"},
                                             {321, "K±"}, {2212, "p/p̄"}};

    for (const auto &[type, count] : sortedTypes) {
        int absType = std::abs(type);
        std::string label = typeLabels.count(absType) ? typeLabels[absType] : "";
        std::cout << "type = " << std::setw(5) << type;
        if (!label.empty()) {
            std::cout << " (" << std::setw(4) << std::left << label << ")";
        } else {
            std::cout << "               ";
        }
        std::cout << "  |  count = " << count << std::endl;
    }

    int chargedTotal = 0, neutralTotal = 0, unknownTotal = 0;
    for (const auto &[type, count] : typeCounter) {
        if (isCharged(type))
            chargedTotal += count;
        else if (isNeutral(type))
            neutralTotal += count;
        else
            unknownTotal += count;
    }

    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Заряженные:  " << chargedTotal << std::endl;
    std::cout << "Нейтральные: " << neutralTotal << std::endl;
    if (unknownTotal > 0)
        std::cout << "Прочие:      " << unknownTotal << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
}

// ───────────────────────────────────────────────────────────────
// Извлечение имени процесса из пути к файлу
// ───────────────────────────────────────────────────────────────
std::string extractProcessName(const std::string &filepath) {
    fs::path p(filepath);
    std::string filename = p.filename().string();

    const std::string prefix = "merged_";
    if (filename.find(prefix) == 0) {
        filename = filename.substr(prefix.length());
    }

    const std::string suffix = ".root";
    if (filename.size() >= suffix.size() &&
        filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
        filename = filename.substr(0, filename.size() - suffix.size());
    }

    return filename;
}

// ───────────────────────────────────────────────────────────────
// Вывод справки
// ───────────────────────────────────────────────────────────────
void printUsage(const char *progName) {
    std::cout << "Использование:\n"
              << "  " << progName << " <input_root_file> [options]\n\n"
              << "Обязательные аргументы:\n"
              << "  input_root_file    Путь к входному ROOT-файлу с деревом анализа\n\n"
              << "Опции:\n"
              << "  -h, --help         Показать эту справку\n"
              << "  -o, --output-dir   Базовая директория для результатов "
              << "(по умолчанию: ../pdf_results)\n"
              << "                     Внутри создаётся поддиректория с именем "
                 "процесса\n\n"
              << "Примеры:\n"
              << "  " << progName << " merged_E240_qqHX.root\n"
              << "  " << progName << " /path/to/merged_E240_qqHinvi.root -o ./plots\n";
}

// ───────────────────────────────────────────────────────────────
int main(int argc, char *argv[]) {
    // ─────────────────────────────────────────────────────────────
    // Парсинг аргументов
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
    // Подготовка путей
    // ─────────────────────────────────────────────────────────────
    std::string processName = extractProcessName(inputRootFile);
    std::cout << "Процесс: " << processName << std::endl;
    std::cout << "Входной файл: " << inputRootFile << std::endl;

    fs::path processOutputDir = fs::path(outputBaseDir) / processName;
    try {
        fs::create_directories(processOutputDir);
        std::cout << "Директория результатов: " << fs::absolute(processOutputDir) << std::endl;
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Ошибка при создании директории " << processOutputDir << ": " << e.what()
                  << std::endl;
        return 1;
    }

    std::string outputPdf =
        (processOutputDir / ("pfo_charged_neutral_" + processName + ".pdf")).string();

    // ─────────────────────────────────────────────────────────────
    // Открытие файла
    // ─────────────────────────────────────────────────────────────
    auto startTime = std::chrono::high_resolution_clock::now();

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
    // Установка ветвей
    // ─────────────────────────────────────────────────────────────
    std::vector<int> *particleType = nullptr;
    std::vector<double> *pfoE = nullptr; // читаем, но не используем для cuts
    tree->SetBranchAddress(BRANCH_PARTICLE_TYPE.c_str(), &particleType);
    tree->SetBranchAddress(BRANCH_PFO_ENERGY.c_str(), &pfoE);

    // ─────────────────────────────────────────────────────────────
    // Создание гистограмм
    // ───────────────────────────────────────────────────────────────
    TH1D *hCharged = new TH1D("hCharged", "Charged PFOs", NUM_BINS, LEFT_BOARD, RIGHT_BOARD);
    TH1D *hNeutral = new TH1D("hNeutral", "Neutral PFOs", NUM_BINS, LEFT_BOARD, RIGHT_BOARD);

    hCharged->SetLineColor(HIST_COLOR_CHARGED);
    hCharged->SetLineWidth(3);
    hCharged->SetFillColorAlpha(HIST_COLOR_CHARGED, 0.25);

    hNeutral->SetLineColor(HIST_COLOR_NEUTRAL);
    hNeutral->SetLineWidth(3);
    hNeutral->SetFillColorAlpha(HIST_COLOR_NEUTRAL, 0.25);

    hCharged->GetXaxis()->SetTitle(X_AXIS_TITLE.c_str());
    hCharged->GetYaxis()->SetTitle(Y_AXIS_TITLE.c_str());
    hCharged->SetTitle(HIST_TITLE.c_str());

    // ─────────────────────────────────────────────────────────────
    // Основной цикл
    // ─────────────────────────────────────────────────────────────
    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessed = 0;

    std::map<int, int> typeCounter;
    Long64_t totalParticles = 0;

    std::cout << "\nНачало обработки событий (без фильтрации по энергии)..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        logProgress(i + 1, nEntries, "Processing");

        if (!particleType)
            continue;

        int nCharged = 0;
        int nNeutral = 0;

        // Проходим по всем частицам
        for (size_t j = 0; j < particleType->size(); ++j) {
            int pfoType = particleType->at(j);

            // Статистика по типам
            typeCounter[pfoType]++;
            totalParticles++;

            // Классификация для гистограмм
            if (isCharged(pfoType)) {
                ++nCharged;
            } else if (isNeutral(pfoType)) {
                ++nNeutral;
            }
            // Частицы, не попавшие в категории, не учитываются в гистограммах
            // (но считаются в totalParticles и typeCounter)
        }

        hCharged->Fill(nCharged);
        hNeutral->Fill(nNeutral);
        nProcessed++;
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
    std::cout << "Заряженные:  среднее = " << hCharged->GetMean() << " ± " << hCharged->GetRMS()
              << std::endl;
    std::cout << "Нейтральные: среднее = " << hNeutral->GetMean() << " ± " << hNeutral->GetRMS()
              << std::endl;
    std::cout << "═══════════════════════════════════════════════════" << std::endl;

    if (PRINT_TYPE_STATISTICS) {
        printTypeStatistics(typeCounter, totalParticles);
    }

    // ─────────────────────────────────────────────────────────────
    // Отрисовка
    // ─────────────────────────────────────────────────────────────
    std::cout << "\nОтрисовка гистограмм..." << std::endl;

    TCanvas *c1 = new TCanvas("c1", HIST_TITLE.c_str(), 900, 650);
    c1->cd();
    c1->SetLogy();

    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadRightMargin(0.07);
    gStyle->SetPadTopMargin(0.07);

    hCharged->Draw("HIST");
    hNeutral->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.65, 0.70, 0.88, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hCharged, "Charged PFOs", "l");
    legend->AddEntry(hNeutral, "Neutral PFOs", "l");
    legend->Draw();

    if (DRAW_MEAN_LINES) {
        TLine *lineCharged =
            new TLine(hCharged->GetMean(), 0, hCharged->GetMean(), hCharged->GetMaximum() * 1.05);
        lineCharged->SetLineColor(HIST_COLOR_CHARGED);
        lineCharged->SetLineWidth(2);
        lineCharged->SetLineStyle(9);
        lineCharged->Draw();

        TLine *lineNeutral =
            new TLine(hNeutral->GetMean(), 0, hNeutral->GetMean(), hNeutral->GetMaximum() * 1.05);
        lineNeutral->SetLineColor(HIST_COLOR_NEUTRAL);
        lineNeutral->SetLineWidth(2);
        lineNeutral->SetLineStyle(9);
        lineNeutral->Draw();
    }

    TLatex *latexCms = new TLatex(0.15, 0.92, "#sqrt{s} = 240 GeV");
    latexCms->SetNDC();
    latexCms->SetTextSize(0.035);
    latexCms->Draw();

    c1->SaveAs(outputPdf.c_str());
    std::cout << "График сохранён: " << outputPdf << std::endl;

    // ─────────────────────────────────────────────────────────────
    // Очистка
    // ─────────────────────────────────────────────────────────────
    delete hCharged;
    delete hNeutral;
    delete legend;
    delete c1;
    inputFile->Close();
    delete inputFile;

    std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
              << std::endl;
    return 0;
}
