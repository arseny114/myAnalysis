// Скрипт для подсчёта фотонов (PDG=22) и построения спектра их энергии.
// Запуск: ./gamma_counter <input_file.root>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "../include/gamma_count_config.h"

namespace fs = std::filesystem;

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
              << "  input_root_file    Путь к входному ROOT-файлу\n\n"
              << "Опции:\n"
              << "  -h, --help         Показать эту справку\n"
              << "  -o, --output-dir   Директория для результатов (по умолчанию: "
              << OUTPUT_BASE_DIR << ")\n\n"
              << "Пример:\n"
              << "  " << progName << " merged_E240_qqHX.root\n";
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

    std::string inputRootFile, outputBaseDir = OUTPUT_BASE_DIR;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -o нужен путь\n";
                return 1;
            }
            outputBaseDir = argv[++i];
        } else if (arg[0] != '-') {
            inputRootFile = arg;
        } else {
            std::cerr << "Неизвестный аргумент: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    if (inputRootFile.empty()) {
        std::cerr << "Ошибка: не указан входной файл\n";
        return 1;
    }

    // ─────────────────────────────────────────────────────────────
    // Подготовка путей
    // ─────────────────────────────────────────────────────────────
    std::string processName = extractProcessName(inputRootFile);
    fs::path processOutputDir = fs::path(outputBaseDir) / processName;
    try {
        fs::create_directories(processOutputDir);
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Ошибка создания директории: " << e.what() << "\n";
        return 1;
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    std::cout << "Входной файл: " << inputRootFile << "\n";

    // ─────────────────────────────────────────────────────────────
    // Открытие файла и дерева
    // ─────────────────────────────────────────────────────────────
    TFile *inputFile = TFile::Open(inputRootFile.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Ошибка: не удалось открыть " << inputRootFile << "\n";
        return 1;
    }
    TTree *tree = dynamic_cast<TTree *>(inputFile->Get(TREE_NAME.c_str()));
    if (!tree) {
        std::cerr << "Ошибка: дерево " << TREE_NAME << " не найдено\n";
        return 1;
    }
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Всего событий: " << nEntries << "\n";

    // ─────────────────────────────────────────────────────────────
    // Установка ветвей
    // ─────────────────────────────────────────────────────────────
    std::vector<double> *pfoE = nullptr, *pfoPx = nullptr, *pfoPy = nullptr, *pfoPz = nullptr;
    std::vector<int> *particleType = nullptr;
    tree->SetBranchAddress(BRANCH_PFO_E.c_str(), &pfoE);
    tree->SetBranchAddress(BRANCH_PFO_PX.c_str(), &pfoPx);
    tree->SetBranchAddress(BRANCH_PFO_PY.c_str(), &pfoPy);
    tree->SetBranchAddress(BRANCH_PFO_PZ.c_str(), &pfoPz);
    tree->SetBranchAddress(BRANCH_PARTICLE_TYPE.c_str(), &particleType);

    // ─────────────────────────────────────────────────────────────
    // Гистограмма и счётчики
    // ─────────────────────────────────────────────────────────────
    TH1F *hSpectrum = new TH1F("hSpectrum", "Photon Energy Spectrum;E_{#gamma} [GeV];Events",
                               SPECTRUM_BINS, SPECTRUM_MIN, SPECTRUM_MAX);
    Long64_t totalLow = 0, totalHigh = 0, withLow = 0, withHigh = 0;

    // ─────────────────────────────────────────────────────────────
    // Цикл по событиям
    // ─────────────────────────────────────────────────────────────
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        int cLow = 0, cHigh = 0;
        if (pfoE && particleType) {
            for (size_t j = 0; j < pfoE->size(); ++j) {
                if (particleType->at(j) != PDG_PHOTON)
                    continue;
                double E = pfoE->at(j);
                hSpectrum->Fill(E);
                if (E > THRESHOLD_LOW) {
                    cLow++;
                    totalLow++;
                }
                if (E > THRESHOLD_HIGH) {
                    cHigh++;
                    totalHigh++;
                }
            }
        }
        if (cLow > 0)
            withLow++;
        if (cHigh > 0)
            withHigh++;
    }

    // ─────────────────────────────────────────────────────────────
    // Вывод таблицы (printf для надёжного выравнивания)
    // ─────────────────────────────────────────────────────────────
    auto endTime = std::chrono::high_resolution_clock::now();
    int totalSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    printf("\n");
    printf("================================================================\n");
    printf("           СТАТИСТИКА ПО ФОТОНАМ (PDG=%d)\n", PDG_PHOTON);
    printf("================================================================\n");
    printf("Событий всего: %lld | Время: %d с\n", nEntries, totalSec);
    printf("----------------------------------------------------------------\n");
    printf("Порог   | Событий с gamma >= 1 | Доля (%%) | Всего gamma | Ср./событие\n");
    printf("----------------------------------------------------------------\n");
    printf("%6s  | %16lld | %8.2f | %9lld | %11.3f\n", "30 ГэВ", withLow,
           100.0 * withLow / nEntries, totalLow, 1.0 * totalLow / nEntries);
    printf("%6s  | %16lld | %8.2f | %9lld | %11.3f\n", "50 ГэВ", withHigh,
           100.0 * withHigh / nEntries, totalHigh, 1.0 * totalHigh / nEntries);
    printf("================================================================\n");

    // ─────────────────────────────────────────────────────────────
    // Отрисовка спектра
    // ─────────────────────────────────────────────────────────────
    std::string pdfPath = (processOutputDir / ("spectrum_" + processName + ".pdf")).string();

    TCanvas *c1 = new TCanvas("c1", "Photon Spectrum", 800, 600);
    c1->SetLogy();
    c1->SetGrid();

    // Рисуем гистограмму уже с logy
    hSpectrum->Draw("HIST");

    // Даём ROOT полностью пересчитать оси
    c1->Update();

    // Теперь берём реальные пределы осей
    double ymin = hSpectrum->GetMinimum();
    double ymax = hSpectrum->GetMaximum();

    TLine *lineLow = new TLine(THRESHOLD_LOW, ymin, THRESHOLD_LOW, ymax);
    lineLow->SetLineColor(kRed);
    lineLow->SetLineStyle(7);
    lineLow->SetLineWidth(2);

    TLine *lineHigh = new TLine(THRESHOLD_HIGH, ymin, THRESHOLD_HIGH, ymax);
    lineHigh->SetLineColor(kBlue);
    lineHigh->SetLineStyle(7);
    lineHigh->SetLineWidth(2);

    lineLow->Draw("SAME");
    lineHigh->Draw("SAME");

    // Легенда
    TLegend *leg = new TLegend(0.5, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(lineLow, Form("E > %.0f GeV", THRESHOLD_LOW), "l");
    leg->AddEntry(lineHigh, Form("E > %.0f GeV", THRESHOLD_HIGH), "l");
    leg->Draw();

    c1->SaveAs(pdfPath.c_str());
    std::cout << "Спектр сохранён: " << pdfPath << "\n";

    // ─────────────────────────────────────────────────────────────
    // Очистка
    // ─────────────────────────────────────────────────────────────
    delete hSpectrum;
    delete lineLow;
    delete lineHigh;
    delete leg;
    delete c1;
    inputFile->Close();
    delete inputFile;

    return 0;
}
