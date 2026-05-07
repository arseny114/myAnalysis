// Анализ для выделения процесса ee -> ZH -> (Z -> qq; H -> invisible)
//
// Применяет последовательность катов:
// 1. Veto на изолированные лептоны
// 2. Veto на высокоэнергетические фотоны
// 3. Требование ровно 2 инклюзивных джета
// 4. Требование минимального числа конституентов в каждом джете (через inclusiveJetSize)
// 5. Окно по инвариантной массе диджетов
// 5. Окно по массе отдачи
//
// Строит гистограммы:
// - Инвариантная масса двух джетов
// - Масса отдачи двух джетов
// - 2D распределение: инвариантная масса vs масса отдачи

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPad.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "../include/zh_invisible_analysis.h"

namespace fs = std::filesystem;

// =============================================================================
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
// =============================================================================

// Логирование прогресса обработки
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

        std::cout << std::endl;
        std::cout.flush();
    }
}

// Расчёт массы отдачи по формуле: M_recoil^2 = (sqrt(s) - E_sys)^2 - |p_sys|^2
double calculateRecoilMass(const TLorentzVector &system, double sqrtS) {
    double recoilE = sqrtS - system.E();
    double recoilP2 =
        system.Px() * system.Px() + system.Py() * system.Py() + system.Pz() * system.Pz();
    double recoilMass2 = recoilE * recoilE - recoilP2;
    return (recoilMass2 > 0) ? std::sqrt(recoilMass2) : 0.0;
}

// Вычисление псевдобыстроты η из четырёхвектора
// η = 0.5 * ln[(|p| + pz) / (|p| - pz)]
double calculatePseudorapidity(const TLorentzVector &vec) {
    double p = vec.P();
    double pz = vec.Pz();
    if (p <= 1e-9)
        return 0.0;
    if (std::abs(pz) >= p)
        return (pz > 0) ? 10.0 : -10.0;
    return 0.5 * std::log((p + pz) / (p - pz));
}

// Вычисление полярного угла θ из четырёхвектора
// θ = arccos(pz / |p|)
double calculatePolarAngle(const TLorentzVector &vec) {
    double p = vec.P();
    double pz = vec.Pz();
    if (p <= 1e-9)
        return 0.0;
    double cosTheta = pz / p;
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // защита от численных ошибок
    return std::acos(cosTheta);
}

// Вычисление ΔR между двумя четырёхвекторами
// ΔR = sqrt( (Δη)² + (Δφ)² )
double calculateDeltaR(const TLorentzVector &v1, const TLorentzVector &v2) {
    double eta1 = calculatePseudorapidity(v1);
    double eta2 = calculatePseudorapidity(v2);
    double phi1 = v1.Phi();
    double phi2 = v2.Phi();

    // Учёт периодичности φ: Δφ ∈ [-π, π]
    double dPhi = phi1 - phi2;
    while (dPhi > M_PI)
        dPhi -= 2 * M_PI;
    while (dPhi < -M_PI)
        dPhi += 2 * M_PI;

    double dEta = eta1 - eta2;
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

// Проверка наличия изолированных лептонов в событии
bool hasIsolatedLepton(const std::vector<int> *flags) {
    if (!flags)
        return false;
    for (int flag : *flags) {
        if (flag == 1)
            return true;
    }
    return false;
}

// Проверка наличия фотонов с энергией выше порога
bool hasHighEnergyPhoton(const std::vector<int> *particleTypes,
                         const std::vector<double> *particleEnergies, double energyCut) {
    if (!particleTypes || !particleEnergies)
        return false;
    if (particleTypes->size() != particleEnergies->size())
        return false;

    for (size_t i = 0; i < particleTypes->size(); ++i) {
        if (std::abs(particleTypes->at(i)) == PDG_PHOTON && particleEnergies->at(i) > energyCut) {
            return true;
        }
    }
    return false;
}

// Отрисовка 1D гистограммы с линиями маркерами
void drawHistogram1D(TH1F *hist, const std::string &canvasTitle, const std::string &xTitle,
                     const std::string &outputFile, double markLine = -1,
                     const std::string &markLabel = "", Color_t markColor = kRed,
                     int lineWidth = 2) {
    TCanvas *c = new TCanvas(canvasTitle.c_str(), canvasTitle.c_str(), 900, 700);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);

    gStyle->SetOptStat(1111);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitle("Events");
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->SetLineWidth(lineWidth);
    hist->Draw("HIST");

    if (markLine > 0) {
        double ymin = hist->GetMinimum();
        double ymax = hist->GetMaximum();
        if (gPad->GetLogy()) {
            ymin = std::max(0.1, ymin * 0.1);
            ymax *= 10;
        }

        TLine *line = new TLine(markLine, ymin, markLine, ymax);
        line->SetLineColor(markColor);
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        line->Draw();

        TLatex *label = new TLatex(markLine + 3, ymax * 0.9, markLabel.c_str());
        label->SetTextColor(markColor);
        label->SetTextSize(0.035);
        label->SetTextAlign(12);
        label->Draw();
    }

    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;
    delete c;
}

// Отрисовка 2D гистограммы с контурами масс
void drawHistogram2D(TH2F *hist, const std::string &canvasTitle, const std::string &xTitle,
                     const std::string &yTitle, const std::string &outputFile, double markX = -1,
                     double markY = -1, const std::string &markLabelX = "",
                     const std::string &markLabelY = "") {
    TCanvas *c = new TCanvas(canvasTitle.c_str(), canvasTitle.c_str(), 900, 800);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.15);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitle("Events");
    hist->GetZaxis()->SetTitleSize(0.045);

    hist->Draw("COLZ");

    if (markX > 0) {
        double ymin = hist->GetYaxis()->GetXmin();
        double ymax = hist->GetYaxis()->GetXmax();
        TLine *lineX = new TLine(markX, ymin, markX, ymax);
        lineX->SetLineColor(kRed);
        lineX->SetLineWidth(2);
        lineX->SetLineStyle(kDashed);
        lineX->Draw();

        TLatex *labelX = new TLatex(markX + 5, ymax * 0.95, markLabelX.c_str());
        labelX->SetTextColor(kRed);
        labelX->SetTextSize(0.03);
        labelX->SetTextAlign(12);
        labelX->Draw();
    }

    if (markY > 0) {
        double xmin = hist->GetXaxis()->GetXmin();
        double xmax = hist->GetXaxis()->GetXmax();
        TLine *lineY = new TLine(xmin, markY, xmax, markY);
        lineY->SetLineColor(kBlue);
        lineY->SetLineWidth(2);
        lineY->SetLineStyle(kDashed);
        lineY->Draw();

        TLatex *labelY = new TLatex(xmax * 0.75, markY + 3, markLabelY.c_str());
        labelY->SetTextColor(kBlue);
        labelY->SetTextSize(0.03);
        labelY->SetTextAlign(12);
        labelY->Draw();
    }

    std::string entriesText = "Entries: " + std::to_string(static_cast<int>(hist->GetEntries()));
    TLatex *entriesLabel = new TLatex(0.15, 0.92, entriesText.c_str());
    entriesLabel->SetNDC();
    entriesLabel->SetTextSize(0.035);
    entriesLabel->Draw();

    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;
    delete c;
}

// Извлечение имени процесса из пути к файлу
// Формат файла: merged_E240_qqHX.root
std::string extractProcessName(const std::string &filepath) {
    fs::path p(filepath);
    std::string filename = p.filename().string();

    // Убираем префикс "merged_"
    const std::string prefix = "merged_";
    if (filename.find(prefix) == 0) {
        filename = filename.substr(prefix.length());
    }

    // Убираем суффикс ".root"
    const std::string suffix = ".root";
    if (filename.size() >= suffix.size() &&
        filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
        filename = filename.substr(0, filename.size() - suffix.size());
    }

    return filename;
}

// Вывод справки по использованию
void printUsage(const char *progName) {
    std::cout << "Анализ для выделения ZH -> qq + invisible\n\n"
              << "Использование:\n"
              << "  " << progName << " <input_root_file> [options]\n\n"
              << "Обязательные аргументы:\n"
              << "  input_root_file    Путь к входному ROOT-файлу с деревом анализа\n\n"
              << "Опции:\n"
              << "  -h, --help         Показать эту справку\n"
              << "  -o, --output-dir   Базовая директория для результатов "
                 "(по умолчанию: ../pdf_results)\n"
              << "  -e, --energy-cut   Порог энергии фотона для veto (ГэВ, по умолчанию: 30)\n"
              << "  -c, --min-const    Мин. число конституентов в джете (по умолчанию: 6)\n"
              << "\nПримеры:\n"
              << "  " << progName << " merged_E240_qqHX.root\n"
              << "  " << progName << " /path/to/file.root -o ./plots -e 50 -c 8\n";
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main(int argc, char *argv[]) {
    // Парсинг аргументов командной строки
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string inputRootFile;
    std::string outputBaseDir = OUTPUT_BASE_DIR;
    double photonEnergyCut = PHOTON_ENERGY_CUT_GEV;
    int minConstituents = MIN_CONSTITUENTS_PER_JET;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -o/--output-dir требуется указать путь" << std::endl;
                return 1;
            }
            outputBaseDir = argv[++i];
        } else if (arg == "-e" || arg == "--energy-cut") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -e/--energy-cut требуется указать значение"
                          << std::endl;
                return 1;
            }
            photonEnergyCut = std::stod(argv[++i]);
        } else if (arg == "-c" || arg == "--min-const") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -c/--min-const требуется указать значение" << std::endl;
                return 1;
            }
            minConstituents = std::stoi(argv[++i]);
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

    // Извлечение имени процесса и подготовка путей
    std::string processName = extractProcessName(inputRootFile);
    std::cout << "Процесс: " << processName << std::endl;
    std::cout << "Входной файл: " << inputRootFile << std::endl;

    fs::path processOutputDir = fs::path(outputBaseDir) / processName;
    try {
        fs::create_directories(processOutputDir);
        std::cout << "Директория результатов: " << fs::absolute(processOutputDir) << std::endl;
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Ошибка при создании директории: " << e.what() << std::endl;
        return 1;
    }

    // Формирование имён выходных файлов
    auto makeOutputPath = [&](const std::string &basename) -> std::string {
        return (processOutputDir / (basename + "_" + processName + ".pdf")).string();
    };

    const std::string OUTPUT_INV_MASS = makeOutputPath("inv_mass_2jets");
    const std::string OUTPUT_RECOIL_MASS = makeOutputPath("recoil_mass_2jets");
    const std::string OUTPUT_2D_CORR = makeOutputPath("inv_vs_recoil_2d");
    const std::string OUTPUT_THETA_Z = makeOutputPath("theta_Z_polar_angle");
    const std::string OUTPUT_DELTA_R = makeOutputPath("deltaR_jet1_jet2");

    // Инициализация ROOT
    gStyle->SetOptStat(1111);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(false);
    auto startTime = std::chrono::high_resolution_clock::now();

    // Открытие входного файла
    std::cout << "\nОткрытие файла: " << inputRootFile << std::endl;
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

    // Установка ветвей для inclusive джетов
    std::vector<double> *inclJetE = nullptr, *inclJetPx = nullptr;
    std::vector<double> *inclJetPy = nullptr, *inclJetPz = nullptr;
    std::vector<double> *inclJetSize = nullptr;

    tree->SetBranchAddress("inclusiveJetE", &inclJetE);
    tree->SetBranchAddress("inclusiveJetPx", &inclJetPx);
    tree->SetBranchAddress("inclusiveJetPy", &inclJetPy);
    tree->SetBranchAddress("inclusiveJetPz", &inclJetPz);
    tree->SetBranchAddress("inclusiveJetSize", &inclJetSize);

    // Ветви для lepton veto и photon veto
    std::vector<int> *isIsolatedLeptonFlag = nullptr;
    std::vector<int> *particleType = nullptr;
    std::vector<double> *pfoE = nullptr;

    if (APPLY_LEPTON_VETO) {
        tree->SetBranchAddress("isIsolatedLeptonFlag", &isIsolatedLeptonFlag);
    }
    if (APPLY_PHOTON_VETO) {
        tree->SetBranchAddress("particleType", &particleType);
        tree->SetBranchAddress("pfoE", &pfoE);
    }

    // Создание гистограмм
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant Mass of Two Jets;M_{jj} [GeV];Events",
                              MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV);

    TH1F *hRecoilMass =
        new TH1F("hRecoilMass", "Recoil Mass Against Two Jets;M_{recoil} [GeV];Events", RECOIL_BINS,
                 RECOIL_MIN_GEV, RECOIL_MAX_GEV);

    TH2F *h2D_Correlation = new TH2F(
        "h2D_Correlation", "Invariant Mass vs Recoil Mass;M_{jj} [GeV];M_{recoil} [GeV]", MASS_BINS,
        MASS_MIN_GEV, MASS_MAX_GEV, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);

    TH1F *hThetaZ =
        new TH1F("hThetaZ", "Polar Angle of Z Boson (Two-Jet System);#theta_{Z} [rad];Events",
                 THETA_BINS, THETA_MIN_RAD, THETA_MAX_RAD);

    TH1F *hDeltaR = new TH1F("hDeltaR",
                             "Distance #Delta R Between Two Jets;#Delta R = #sqrt{#Delta#eta^{2} + "
                             "#Delta#phi^{2}};Events",
                             DELTA_R_BINS, DELTA_R_MIN, DELTA_R_MAX);

    // Статистика прохождения катов
    CutStatistics stats;

    // Основной цикл по событиям
    Long64_t nEntries = tree->GetEntries();
    std::cout << "\nНачало обработки событий..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        logProgress(i + 1, nEntries, "Processing");
        stats.totalEvents++;

        // Cut 1: Veto на изолированные лептоны
        if (APPLY_LEPTON_VETO && hasIsolatedLepton(isIsolatedLeptonFlag)) {
            continue;
        }
        stats.afterLeptonVeto++;

        // Cut 2: Veto на высокоэнергетические фотоны
        if (APPLY_PHOTON_VETO && hasHighEnergyPhoton(particleType, pfoE, photonEnergyCut)) {
            continue;
        }
        stats.afterPhotonVeto++;

        // Cut 3: Требование ровно 2 инклюзивных джета
        if (REQUIRE_EXACTLY_TWO_INCLUSIVE_JETS && inclJetE->size() != 2) {
            continue;
        }
        stats.afterJetCount++;

        // Cut 4: Требование минимального числа конституентов в каждом из двух джетов
        if (REQUIRE_MIN_CONSTITUENTS_PER_JET) {
            int nConst1 = static_cast<int>(inclJetSize->at(0));
            int nConst2 = static_cast<int>(inclJetSize->at(1));

            if (nConst1 < minConstituents || nConst2 < minConstituents) {
                continue;
            }
        }
        stats.afterConstituents++;

        // Расчёт инвариантной и массы отдачи для двух джетов
        TLorentzVector jet1(inclJetPx->at(0), inclJetPy->at(0), inclJetPz->at(0), inclJetE->at(0));
        TLorentzVector jet2(inclJetPx->at(1), inclJetPy->at(1), inclJetPz->at(1), inclJetE->at(1));

        TLorentzVector dijet = jet1 + jet2;
        double invMass = dijet.M();
        double recoilMass = calculateRecoilMass(dijet, SQRT_S_GEV);

        // Расчёт θ и ΔR
        double thetaZ = calculatePolarAngle(dijet);  // полярный угол системы двух джетов
        double deltaR = calculateDeltaR(jet1, jet2); // расстояние между джетами

        // Cut 5: Окно инвариантной массы диджета
        if (APPLY_DIJET_MASS_WINDOW) {
            if (invMass < DIJET_MASS_WINDOW_MIN_GEV || invMass > DIJET_MASS_WINDOW_MAX_GEV) {
                continue;
            }
        }
        stats.afterDijetMassWindow++;

        // Cut 6: Окно массы отдачи
        if (APPLY_RECOIL_MASS_WINDOW) {
            if (recoilMass < RECOIL_MASS_WINDOW_MIN_GEV ||
                recoilMass > RECOIL_MASS_WINDOW_MAX_GEV) {
                continue;
            }
        }
        stats.afterRecoilMassWindow++;

        // Заполнение гистограмм
        hInvMass->Fill(invMass);
        hRecoilMass->Fill(recoilMass);
        h2D_Correlation->Fill(invMass, recoilMass);
        hThetaZ->Fill(thetaZ);
        hDeltaR->Fill(deltaR);

        stats.finalSelected++;
    }

    // Итоговая статистика
    auto endTime = std::chrono::high_resolution_clock::now();
    auto totalSec = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    std::cout << "\n=======================================================" << std::endl;
    std::cout << "Обработка завершена!" << std::endl;
    std::cout << "Всего событий: " << nEntries << std::endl;
    std::cout << "Прошло времени: " << totalSec << " с (" << totalSec / 60.0 << " мин)"
              << std::endl;

    stats.print(processName);

    // Отрисовка и сохранение гистограмм (только PDF)
    std::cout << "Отрисовка гистограмм..." << std::endl;

    drawHistogram1D(hInvMass, "cInvMass", "M_{jj} [GeV]", OUTPUT_INV_MASS, MZ_GEV,
                    "M_{Z} = 91.2 GeV", kRed, 2);

    drawHistogram1D(hRecoilMass, "cRecoilMass", "M_{recoil} [GeV]", OUTPUT_RECOIL_MASS, MH_GEV,
                    "M_{H} = 125.3 GeV", kBlue, 2);

    drawHistogram2D(h2D_Correlation, "c2D_Correlation", "M_{jj} [GeV]", "M_{recoil} [GeV]",
                    OUTPUT_2D_CORR, MZ_GEV, MH_GEV, "M_{Z}", "M_{H}");

    drawHistogram1D(hThetaZ, "cThetaZ", "#theta_{Z} [rad]", OUTPUT_THETA_Z, -1, "", kGreen + 2, 2);

    drawHistogram1D(hDeltaR, "cDeltaR", "#Delta R", OUTPUT_DELTA_R, -1, "", kMagenta, 2);

    // Очистка памяти
    delete hInvMass;
    delete hRecoilMass;
    delete h2D_Correlation;
    delete hThetaZ;
    delete hDeltaR;
    inputFile->Close();
    delete inputFile;

    std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
              << std::endl;
    return 0;
}
