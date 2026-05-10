// Анализ для выделения процесса ee -> ZH -> (Z -> qq; H -> invisible)
//
// Применяет последовательность катов:
// 1. Veto на изолированные лептоны
// 2. Veto на высокоэнергетические фотоны
// 3. Требование ровно 2 инклюзивных джета
// 4. Требование минимального числа конституентов в каждом джете (через inclusiveJetSize)
// 5. Окно по инвариантной массе диджетов
// 6. Окно по массе отдачи
//
// Строит гистограммы:
// - Инвариантная масса двух джетов
// - Масса отдачи двух джетов
// - Полярный угол системы двух джетов
// - Расстояние deltaR между джетами
// - 2D распределение: инвариантная масса vs масса отдачи
// - 2D распределение: E_photon(>PHOTON_ENERGY_CUT_GEV) vs M_recoil

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
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

// Проверка, находится ли точка (x,y) внутри повёрнутого эллипса
// Параметры: центр (cx,cy), полуоси (a,b), угол поворота theta (радианы)
bool isInsideEllipse(double x, double y, double cx, double cy, double a, double b, double theta) {
    // 1. Перенос начала координат в центр эллипса
    double dx = x - cx;
    double dy = y - cy;

    // 2. Поворот системы координат на угол -theta (чтобы совместить с осями эллипса)
    double cosT = std::cos(theta);
    double sinT = std::sin(theta);
    double xRot = dx * cosT + dy * sinT;
    double yRot = -dx * sinT + dy * cosT;

    // 3. Проверка канонического уравнения эллипса: (x'/a)² + (y'/b)² <= 1
    double value = (xRot * xRot) / (a * a) + (yRot * yRot) / (b * b);
    return (value <= 1.0);
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

// Вычисление энергии в конусе вокруг произвольной PFO (исключая саму частицу)
double calculateConeEnergy(size_t centerIdx, const std::vector<double> *pfoE,
                           const std::vector<double> *pfoPx, const std::vector<double> *pfoPy,
                           const std::vector<double> *pfoPz, double cosConeCut) {
    if (!pfoE || !pfoPx || !pfoPy || !pfoPz)
        return 0.0;
    if (centerIdx >= pfoE->size())
        return 0.0;

    double coneE = 0.0;
    double px1 = pfoPx->at(centerIdx);
    double py1 = pfoPy->at(centerIdx);
    double pz1 = pfoPz->at(centerIdx);
    double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);

    if (p1 < 1e-9)
        return 0.0; // защита от деления на ноль

    for (size_t i = 0; i < pfoE->size(); ++i) {
        if (i == centerIdx)
            continue; // исключаем саму центральную частицу

        double px2 = pfoPx->at(i);
        double py2 = pfoPy->at(i);
        double pz2 = pfoPz->at(i);
        double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
        if (p2 < 1e-9)
            continue;

        // Вычисляем косинус угла между импульсами
        double cosTheta = (px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2);
        // Защита от численных ошибок (косинус может выйти за [-1, 1])
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

        if (cosTheta >= cosConeCut) {
            coneE += pfoE->at(i);
        }
    }
    return coneE;
}

// Проверка наличия изолированных фотонов
bool hasIsolatedPhoton(const std::vector<int> *particleTypes, const std::vector<double> *pfoE,
                       const std::vector<double> *pfoPx, const std::vector<double> *pfoPy,
                       const std::vector<double> *pfoPz, double minEnergy, double cosConeCut,
                       double maxConeEnergy) {
    if (!particleTypes || !pfoE || !pfoPx || !pfoPy || !pfoPz)
        return false;
    if (particleTypes->size() != pfoE->size())
        return false;

    for (size_t i = 0; i < particleTypes->size(); ++i) {
        // Проверяем, что это фотон и его энергия выше порога
        if (std::abs(particleTypes->at(i)) == PDG_PHOTON && pfoE->at(i) > minEnergy) {
            // Вычисляем энергию в конусе вокруг фотона
            double coneE = calculateConeEnergy(i, pfoE, pfoPx, pfoPy, pfoPz, cosConeCut);

            // Фотон считается изолированным, если энергия в конусе меньше порога
            if (coneE < maxConeEnergy) {
                return true;
            }
        }
    }
    return false;
}

// Проверка изоляции одного лептона по алгоритму Gaudi (прямоугольные критерии)
bool isLeptonIsolatedROOT(size_t idx, const std::vector<int> *types,
                          const std::vector<double> *energies, const std::vector<double> *px,
                          const std::vector<double> *py, const std::vector<double> *pz) {
    if (!types || !energies || !px || !py || !pz || idx >= types->size())
        return false;

    int pdg = std::abs(types->at(idx));
    if (pdg != 11 && pdg != 13)
        return false; // только e/mu

    double trackE = energies->at(idx);
    // Прямоугольные критерии по энергии трека
    if (trackE < LEPTON_ISO_MIN_TRACK_E_GEV || trackE > LEPTON_ISO_MAX_TRACK_E_GEV)
        return false;

    // Энергия в конусе (используем уже существующую функцию из вашего кода)
    double coneE = calculateConeEnergy(idx, energies, px, py, pz, LEPTON_ISO_COS_CONE_ANGLE);
    // Прямоугольные критерии по энергии в конусе
    if (coneE < LEPTON_ISO_MIN_CONE_E_GEV || coneE > LEPTON_ISO_MAX_CONE_E_GEV)
        return false;

    return true; // Лептон считается изолированным
}

// Проверка наличия хотя бы одного изолированного лептона в событии
bool hasIsolatedLeptonROOT(const std::vector<int> *types, const std::vector<double> *energies,
                           const std::vector<double> *px, const std::vector<double> *py,
                           const std::vector<double> *pz) {
    if (!types || !energies)
        return false;
    for (size_t i = 0; i < types->size(); ++i) {
        if (isLeptonIsolatedROOT(i, types, energies, px, py, pz))
            return true;
    }
    return false;
}

// Детальная таблица изоляции лептонов в событии
void printLeptonIsolationDebugTable(size_t eventNum, const std::vector<int> *types,
                                    const std::vector<double> *energies,
                                    const std::vector<double> *px, const std::vector<double> *py,
                                    const std::vector<double> *pz,
                                    const std::vector<int> *gaudiFlags) {
    if (!types || !energies || !px || !py || !pz || !gaudiFlags)
        return;
    if (gaudiFlags->size() != types->size())
        return;

    std::cout << "\n[DEBUG TABLE] Event " << eventNum << " | Lepton Isolation Check:\n";
    std::cout << std::left << std::setw(5) << "Idx" << std::setw(6) << "PDG" << std::setw(10)
              << "E(GeV)" << std::setw(8) << "Gaudi" << std::setw(12) << "ConeE(GeV)"
              << std::setw(8) << "ROOT"
              << "\n";
    std::cout << std::string(49, '-') << "\n";

    for (size_t i = 0; i < types->size(); ++i) {
        int pdg = std::abs(types->at(i));
        if (pdg != 11 && pdg != 13)
            continue; // Только электроны и мюоны

        double E = energies->at(i);
        int gFlag = gaudiFlags->at(i);
        double coneE = calculateConeEnergy(i, energies, px, py, pz, LEPTON_ISO_COS_CONE_ANGLE);

        // Прямое повторение логики isLeptonIsolatedROOT
        bool rootIso = (E >= LEPTON_ISO_MIN_TRACK_E_GEV && E <= LEPTON_ISO_MAX_TRACK_E_GEV &&
                        coneE >= LEPTON_ISO_MIN_CONE_E_GEV && coneE <= LEPTON_ISO_MAX_CONE_E_GEV);

        std::cout << std::left << std::setw(5) << i << std::setw(6) << pdg << std::setw(10)
                  << std::fixed << std::setprecision(2) << E << std::setw(8) << gFlag
                  << std::setw(12) << std::fixed << std::setprecision(3) << coneE << std::setw(8)
                  << (rootIso ? "ISO" : "clean") << "\n";
    }
    std::cout << std::string(49, '-') << "\n\n";
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

// Отрисовка 2D гистограммы с опциональным эллипсом
void drawHistogram2D(TH2F *hist, const std::string &canvasTitle, const std::string &xTitle,
                     const std::string &yTitle, const std::string &outputFile, double markX = -1,
                     double markY = -1, const std::string &markLabelX = "",
                     const std::string &markLabelY = "", double ellipseCx = -1,
                     double ellipseCy = -1, double ellipseA = -1, double ellipseB = -1,
                     double ellipseTheta = 0, bool drawEllipse = false) {
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

    // Отрисовка эллипса
    if (drawEllipse && ellipseA > 0 && ellipseB > 0) {
        const int nPts = 120;
        std::vector<double> ex(nPts + 1), ey(nPts + 1);
        for (int i = 0; i <= nPts; ++i) {
            double phi = 2.0 * M_PI * i / nPts;
            // Параметризация в собственной СК эллипса
            double xl = ellipseA * std::cos(phi);
            double yl = ellipseB * std::sin(phi);
            // Поворот и перенос в лабораторную СК
            ex[i] = ellipseCx + xl * std::cos(ellipseTheta) - yl * std::sin(ellipseTheta);
            ey[i] = ellipseCy + xl * std::sin(ellipseTheta) + yl * std::cos(ellipseTheta);
        }
        TPolyLine *ellipse = new TPolyLine(nPts + 1, ex.data(), ey.data());
        ellipse->SetLineColor(kGreen);
        ellipse->SetLineWidth(3);
        ellipse->SetLineStyle(kDashed);
        ellipse->Draw("L SAME");

        // Подпись
        TLatex *ellipseLabel = new TLatex(ellipseCx + 8, ellipseCy - 12, "Ellipse cut");
        ellipseLabel->SetTextColor(kGreen);
        ellipseLabel->SetTextSize(0.03);
        ellipseLabel->Draw();

        // Создаём рамку с параметрами
        double textX = 0.48, textY = 0.7;
        double boxWidth = 0.25, boxHeight = 0.18;

        // Фон под текстом
        auto *paramBox = new TPaveText(textX, textY, textX + boxWidth, textY + boxHeight, "NDC NB");
        paramBox->SetFillColor(kWhite);
        paramBox->SetFillStyle(1001);
        paramBox->SetLineColor(kGray + 2);
        paramBox->SetLineWidth(1);
        paramBox->SetTextAlign(12);
        paramBox->SetTextSize(0.028);

        // Формируем строки с параметрами
        double thetaDeg = ellipseTheta * 180.0 / M_PI;
        paramBox->AddText("Ellipse cut:");
        paramBox->AddText(Form("Center: (%.1f, %.1f) GeV", ellipseCx, ellipseCy));
        paramBox->AddText(Form("Semi-axes: a=%.2f, b=%.2f GeV", ellipseA, ellipseB));
        paramBox->AddText(Form("Rotation: %.1f#circ", thetaDeg));

        paramBox->Draw();
    }

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
    const std::string OUTPUT_PHOTON_E_VS_RECOIL = makeOutputPath("photonE_vs_recoil_2d");
    const std::string OUTPUT_COS_THETA_JET = makeOutputPath("cosTheta_jets");
    const std::string OUTPUT_MET_PFO = makeOutputPath("MET_pfo");
    const std::string OUTPUT_MET_JET = makeOutputPath("MET_jets");

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
    std::vector<double> *pfoPx = nullptr;
    std::vector<double> *pfoPy = nullptr;
    std::vector<double> *pfoPz = nullptr;

    if (APPLY_LEPTON_VETO) {
        tree->SetBranchAddress("isIsolatedLeptonFlag", &isIsolatedLeptonFlag);
    }
    if (DEBUG_LEPTON_ISOLATION || APPLY_HIGH_E_PHOTON_VETO || APPLY_ISOLATED_PHOTON_VETO) {
        tree->SetBranchAddress("particleType", &particleType);
        tree->SetBranchAddress("pfoE", &pfoE);
        tree->SetBranchAddress("pfoPx", &pfoPx);
        tree->SetBranchAddress("pfoPy", &pfoPy);
        tree->SetBranchAddress("pfoPz", &pfoPz);
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

    TH2F *hPhotonE_vs_Recoil =
        new TH2F("hPhotonE_vs_Recoil",
                 "Photon Energy vs Recoil Mass (pre-veto);E_{#gamma} [GeV];M_{recoil} [GeV]",
                 PHOTON_E_BINS, PHOTON_E_MIN_GEV, PHOTON_E_MAX_GEV, RECOIL_MASS_2D_BINS,
                 RECOIL_MASS_2D_MIN_GEV, RECOIL_MASS_2D_MAX_GEV);

    TH1F *hCosThetaJet = new TH1F("hCosThetaJet", "cos#theta of Jets;cos#theta;Events",
                                  COS_THETA_JET_BINS, COS_THETA_JET_MIN, COS_THETA_JET_MAX);

    TH1F *hMETpfo = new TH1F("hMETpfo", "MET from all PFOs;MET_{PFO} [GeV];Events", MET_PFO_BINS,
                             MET_PFO_MIN, MET_PFO_MAX);
    TH1F *hMETjet = new TH1F("hMETjet", "MET from Two Jets;MET_{jet} [GeV];Events", MET_JET_BINS,
                             MET_JET_MIN, MET_JET_MAX);

    // Статистика прохождения катов
    CutStatistics stats;

    // Основной цикл по событиям
    Long64_t nEntries = tree->GetEntries();
    std::cout << "\nНачало обработки событий..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        logProgress(i + 1, nEntries, "Processing");
        stats.totalEvents++;

        // Cut 1: Veto на изолированные лептоны (сравнение Gaudi флага и ROOT-расчета)
        if (APPLY_LEPTON_VETO) {
            bool gaudiHasIsoLepton = hasIsolatedLepton(isIsolatedLeptonFlag);
            bool rootHasIsoLepton = false;

            if (DEBUG_LEPTON_ISOLATION) {
                rootHasIsoLepton = hasIsolatedLeptonROOT(particleType, pfoE, pfoPx, pfoPy, pfoPz);

                // Печатаем таблицу при срабатывании вето или при рассинхроне флагов (лимит 20
                // событий)
                // static int debugPrintCount = 0;
                // if (rootHasIsoLepton || gaudiHasIsoLepton != rootHasIsoLepton) { // &&
                // debugPrintCount < 20) {
                // debugPrintCount++;
                printLeptonIsolationDebugTable(i + 1, particleType, pfoE, pfoPx, pfoPy, pfoPz,
                                               isIsolatedLeptonFlag);
                // }

                // В отладке вето срабатывает по независимому ROOT-расчёту
                if (rootHasIsoLepton)
                    continue;
            } else {
                // Стандартное поведение: veto по флагу из процессора
                if (gaudiHasIsoLepton)
                    continue;
            }
        }
        stats.afterLeptonVeto++;

        // Cut 2: Требование ровно 2 инклюзивных джета
        if (REQUIRE_EXACTLY_TWO_INCLUSIVE_JETS && inclJetE->size() != 2) {
            continue;
        }
        stats.afterJetCount++;

        // Cut 3: Требование минимального числа конституентов в каждом из двух джетов
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

        // Расчёт косинусов углов для джетов
        double cosTheta1 = (jet1.P() > 1e-9) ? jet1.Pz() / jet1.P() : 0.0;
        double cosTheta2 = (jet2.P() > 1e-9) ? jet2.Pz() / jet2.P() : 0.0;

        // Заполнение 2D гистограммы E_photon vs M_recoil (до photon veto)
        if (particleType && pfoE && particleType->size() == pfoE->size()) {
            for (size_t ip = 0; ip < particleType->size(); ++ip) {
                // Проверяем, что это фотон и его энергия выше порога
                if (std::abs(particleType->at(ip)) == PDG_PHOTON &&
                    pfoE->at(ip) > PHOTON_ENERGY_CUT_GEV) {
                    hPhotonE_vs_Recoil->Fill(pfoE->at(ip), recoilMass);
                }
            }
        }

        // Расчет Missing Transverse Energy
        double met_jet = dijet.Pt();
        double met_pfo = 0.0;
        if (pfoPx && pfoPy) {
            double sumPx = 0.0, sumPy = 0.0;
            for (size_t i = 0; i < pfoPx->size(); ++i) {
                sumPx += pfoPx->at(i);
                sumPy += pfoPy->at(i);
            }
            met_pfo = std::sqrt(sumPx * sumPx + sumPy * sumPy);
        }

        // Cut 4: Veto на высокоэнергетические фотоны
        if (APPLY_HIGH_E_PHOTON_VETO && hasHighEnergyPhoton(particleType, pfoE, photonEnergyCut)) {
            continue;
        }
        stats.afterHighEPhotonVeto++;

        // Cut 5: Veto на изолированные фотоны (критерии как для лептонов)
        if (APPLY_ISOLATED_PHOTON_VETO &&
            hasIsolatedPhoton(particleType, pfoE, pfoPx, pfoPy, pfoPz, PHOTON_ISO_MIN_ENERGY_GEV,
                              PHOTON_ISO_COS_CONE_ANGLE, PHOTON_ISO_MAX_CONE_ENERGY_GEV)) {
            continue;
        }
        stats.afterIsoPhotonVeto++;

        // Cut 6: Окно инвариантной массы диджета
        if (APPLY_DIJET_MASS_WINDOW) {
            if (invMass < DIJET_MASS_WINDOW_MIN_GEV || invMass > DIJET_MASS_WINDOW_MAX_GEV) {
                continue;
            }
        }
        stats.afterDijetMassWindow++;

        // Cut 7: Кат на |cos(theta_Z)| < 0.98 (угловое распределение Z-бозона)
        if (APPLY_COS_THETA_Z_CUT) {
            double cosThetaZ = std::cos(thetaZ);
            if (std::abs(cosThetaZ) >= COS_THETA_Z_CUT) {
                continue;
            }
        }
        stats.afterCosThetaZCut++;

        // Cut 8: Окно массы отдачи
        if (APPLY_RECOIL_MASS_WINDOW) {
            if (recoilMass < RECOIL_MASS_WINDOW_MIN_GEV ||
                recoilMass > RECOIL_MASS_WINDOW_MAX_GEV) {
                continue;
            }
        }
        stats.afterRecoilMassWindow++;

        // Cut 9: Эллиптический cut на M_jj vs M_recoil
        if (APPLY_ELLIPSE_CUT) {
            if (!isInsideEllipse(invMass, recoilMass, ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV,
                                 ELLIPSE_B_GEV, ELLIPSE_THETA)) {
                continue;
            }
        }
        stats.afterEllipseCut++;

        // Заполнение гистограмм
        hInvMass->Fill(invMass);
        hRecoilMass->Fill(recoilMass);
        h2D_Correlation->Fill(invMass, recoilMass);
        hThetaZ->Fill(thetaZ);
        hDeltaR->Fill(deltaR);
        hCosThetaJet->Fill(cosTheta1);
        hCosThetaJet->Fill(cosTheta2);
        hMETpfo->Fill(met_pfo);
        hMETjet->Fill(met_jet);

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
                    OUTPUT_2D_CORR, MZ_GEV, MH_GEV, "M_{Z}", "M_{H}",
#if APPLY_ELLIPSE_CUT
                    ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV, ELLIPSE_B_GEV, ELLIPSE_THETA,
                    true
#else
                    -1, -1, -1, -1, 0, false
#endif
    );

    drawHistogram1D(hThetaZ, "cThetaZ", "#theta_{Z} [rad]", OUTPUT_THETA_Z, -1, "", kGreen + 2, 2);

    drawHistogram1D(hDeltaR, "cDeltaR", "#Delta R", OUTPUT_DELTA_R, -1, "", kMagenta, 2);

    drawHistogram2D(hPhotonE_vs_Recoil, "cPhotonE_vs_Recoil", "E_{#gamma} [GeV] (pre-veto)",
                    "M_{recoil} [GeV]", OUTPUT_PHOTON_E_VS_RECOIL, -1, -1, "", "");

    drawHistogram1D(hCosThetaJet, "cCosThetaJet", "cos#theta", OUTPUT_COS_THETA_JET, -1, "", kCyan,
                    2);

    drawHistogram1D(hMETpfo, "cMETpfo", "MET_{PFO} [GeV]", OUTPUT_MET_PFO, -1, "", kOrange + 1, 2);
    drawHistogram1D(hMETjet, "cMETjet", "MET_{jet} [GeV]", OUTPUT_MET_JET, -1, "", kViolet, 2);

    // Очистка памяти
    delete hInvMass;
    delete hRecoilMass;
    delete h2D_Correlation;
    delete hThetaZ;
    delete hDeltaR;
    delete hPhotonE_vs_Recoil;
    delete hCosThetaJet;
    delete hMETpfo;
    delete hMETjet;
    inputFile->Close();
    delete inputFile;

    std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
              << std::endl;
    return 0;
}
