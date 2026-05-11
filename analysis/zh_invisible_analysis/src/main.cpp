// Анализ для выделения процесса ee -> ZH -> (Z -> qq; H -> invisible)
//
// Поддерживает три режима работы:
// 1. NORMAL (по умолчанию): применяет предотборы и основные отборы, строит гистограммы
// 2. TRAINING: выполняет только предотборы, сохраняет фичи в ROOT-файл для обучения ML
// 3. INFERENCE: пропускает предотборы и основные отборы, применяет ML-классификатор,
//              строит гистограммы для отобранных событий
//
// Кинематические переменные:
// - Инвариантная масса двух джетов (M_jj)
// - Масса отдачи (M_recoil)
// - Полярный угол системы двух джетов (cosθ_Z)
// - Расстояние ΔR между джетами
// - Поперечный импульс диджетной системы (MET_jet ≈ p_T^miss)
// - Недостающий 3-импульс и его направление
//
// Отрисовка: 1D и 2D гистограммы с опциональными линиями отборов и эллипсом

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

// Вычисление энергии в конусе, исключая фотоны (для подавления FSR-эффекта)
double calculateConeEnergyExclPhotons(size_t centerIdx, const std::vector<int> *types,
                                      const std::vector<double> *pfoE,
                                      const std::vector<double> *px, const std::vector<double> *py,
                                      const std::vector<double> *pz, double cosConeCut) {
    if (!types || !pfoE || !px || !py || !pz || centerIdx >= pfoE->size())
        return 0.0;
    double coneE = 0.0;
    double px1 = px->at(centerIdx), py1 = py->at(centerIdx), pz1 = pz->at(centerIdx);
    double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
    if (p1 < 1e-9)
        return 0.0;

    for (size_t i = 0; i < pfoE->size(); ++i) {
        if (i == centerIdx)
            continue;
        // Игнорируем фотоны
        if (std::abs(types->at(i)) == PDG_PHOTON)
            continue;

        double px2 = px->at(i), py2 = py->at(i), pz2 = pz->at(i);
        double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
        if (p2 < 1e-9)
            continue;

        double cosTheta = (px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
        if (cosTheta >= cosConeCut)
            coneE += pfoE->at(i);
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

// Проверка изоляции одного лептона с игнорированием фотонов в конусе
bool isLeptonIsolatedROOT_FSR(size_t idx, const std::vector<int> *types,
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

    // Энергия в конусе без учёта фотонов
    double coneE =
        calculateConeEnergyExclPhotons(idx, types, energies, px, py, pz, LEPTON_ISO_COS_CONE_ANGLE);
    if (coneE < LEPTON_ISO_MIN_CONE_E_GEV || coneE > LEPTON_ISO_MAX_CONE_E_GEV)
        return false;

    return true;
}

// Проверка наличия изолированного лептона в событии
bool hasIsolatedLeptonROOT_FSR(const std::vector<int> *types, const std::vector<double> *energies,
                               const std::vector<double> *px, const std::vector<double> *py,
                               const std::vector<double> *pz) {
    if (!types || !energies)
        return false;
    for (size_t i = 0; i < types->size(); ++i) {
        if (isLeptonIsolatedROOT_FSR(i, types, energies, px, py, pz))
            return true;
    }
    return false;
}

// Отрисовка 1D гистограммы с линиями маркерами
void drawHistogram1D(TH1F *hist, const std::string &canvasTitle, const std::string &xTitle,
                     const std::string &outputFile,
                     const std::vector<std::pair<double, std::string>> &markLines = {},
                     Color_t markColor = kRed, int lineWidth = 2) {
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

    // Вызываем Update() чтобы ROOT финализировал оси, затем берём
    // реальные границы из gPad, т.к. они учитывают автомасштаб
    c->Update();
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();
    double xmin = gPad->GetUxmin();
    double xmax = gPad->GetUxmax();
    double xRange = xmax - xmin;

    for (const auto &mark : markLines) {
        TLine *line = new TLine(mark.first, ymin, mark.first, ymax);
        line->SetLineColor(markColor);
        line->SetLineWidth(lineWidth);
        line->SetLineStyle(kDashed);
        line->Draw();

        // Смещение подписи в единицах оси X
        TLatex *label = new TLatex(mark.first + xRange * 0.015, ymax * 0.88, mark.second.c_str());
        label->SetTextColor(markColor);
        label->SetTextSize(0.035);
        label->SetTextAlign(12);
        label->Draw();
    }

    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;
    delete c;
}

// Отрисовка 2D гистограммы с опциональными линиями и эллипсом
void drawHistogram2D(TH2F *hist, const std::string &canvasTitle, const std::string &xTitle,
                     const std::string &yTitle, const std::string &outputFile, double markX = -1,
                     double markY = -1, const std::string &markLabelX = "",
                     const std::string &markLabelY = "", double ellipseCx = -1,
                     double ellipseCy = -1, double ellipseA = -1, double ellipseB = -1,
                     double thetaDeg = 0, bool drawEllipse = false) {
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

    // Берём границы из самой гистограммы (для 2D они фиксированы при создании)
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();
    double ymin = hist->GetYaxis()->GetXmin();
    double ymax = hist->GetYaxis()->GetXmax();
    double xRange = xmax - xmin;
    double yRange = ymax - ymin;

    if (markX > 0) {
        TLine *lineX = new TLine(markX, ymin, markX, ymax);
        lineX->SetLineColor(kRed);
        lineX->SetLineWidth(2);
        lineX->SetLineStyle(kDashed);
        lineX->Draw();

        TLatex *labelX =
            new TLatex(markX + xRange * 0.02, ymax - yRange * 0.06, markLabelX.c_str());
        labelX->SetTextColor(kRed);
        labelX->SetTextSize(0.03);
        labelX->SetTextAlign(12);
        labelX->Draw();
    }

    if (markY > 0) {
        TLine *lineY = new TLine(xmin, markY, xmax, markY);
        lineY->SetLineColor(kBlue);
        lineY->SetLineWidth(2);
        lineY->SetLineStyle(kDashed);
        lineY->Draw();

        TLatex *labelY =
            new TLatex(xmax - xRange * 0.25, markY + yRange * 0.02, markLabelY.c_str());
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
            double xl = ellipseA * std::cos(phi);
            double yl = ellipseB * std::sin(phi);
            double ellipseTheta = (thetaDeg * M_PI) / 180.0;
            ex[i] = ellipseCx + xl * std::cos(ellipseTheta) - yl * std::sin(ellipseTheta);
            ey[i] = ellipseCy + xl * std::sin(ellipseTheta) + yl * std::cos(ellipseTheta);
        }
        TPolyLine *ellipse = new TPolyLine(nPts + 1, ex.data(), ey.data());
        ellipse->SetLineColor(kGreen + 1);
        ellipse->SetLineWidth(3);
        ellipse->SetLineStyle(kDashed);
        ellipse->Draw("L SAME");

        TLatex *ellipseLabel =
            new TLatex(ellipseCx + xRange * 0.04, ellipseCy - yRange * 0.05, "Ellipse cut");
        ellipseLabel->SetTextColor(kGreen + 1);
        ellipseLabel->SetTextSize(0.03);
        ellipseLabel->Draw();

        double textX = 0.48, textY = 0.7;
        double boxWidth = 0.25, boxHeight = 0.18;
        auto *paramBox = new TPaveText(textX, textY, textX + boxWidth, textY + boxHeight, "NDC NB");
        paramBox->SetFillColor(kWhite);
        paramBox->SetFillStyle(1001);
        paramBox->SetLineColor(kGray + 2);
        paramBox->SetLineWidth(1);
        paramBox->SetTextAlign(12);
        paramBox->SetTextSize(0.028);
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

// Заглушка для ML-инференса
#ifdef USE_XGBOOST_INFERENCE
#include "../ml/xgboost_model.h"
bool predictML(const std::vector<double> &features) {
    // xgboost_model.h экспортирует функцию predict_xgboost(features.data())
    return predict_xgboost(features.data()) >= 0.5;
}
#else
bool predictML(const std::vector<double> &) { return true; } // Заглушка для компиляции без модели
#endif

// Вывод справки по использованию
void printUsage(const char *progName) {
    std::cout << "Анализ для выделения ZH -> qq + invisible\n\n"
              << "Использование:\n"
              << "  " << progName << " <input_root_file> [options]\n\n"
              << "Обязательные аргументы:\n"
              << "  input_root_file    Путь к входному ROOT-файлу с деревом анализа\n\n"
              << "Опции:\n"
              << "  -h, --help         Показать эту справку\n"
              << "  -m, --mode         Режим работы: normal (по умолчанию), train, infer\n"
              << "  -l, --label        Метка класса для TRAINING: 0=фон, 1=сигнал\n"
              << "\nПримеры:\n"
              << "  # Обычный режим (предотбор + основные отборы):\n"
              << "  " << progName << " merged_E240_qqHX.root\n\n"
              << "  # Генерация данных для обучения (сигнал):\n"
              << "  " << progName << " merged_E240_qqHinvi.root -m train -l 1\n\n"
              << "  # Генерация данных для обучения (фон):\n"
              << "  " << progName << " merged_E240_qqHX.root -m train -l 0\n\n"
              << "  # Инференс с моделью:\n"
              << "  " << progName << " merged_E240_qqHX.root -m infer\n";
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

    AnalysisMode mode = AnalysisMode::NORMAL;
    std::string inputRootFile;
    std::string labelStr = "0"; // 0=фон, 1=сигнал (для TRAINING)
    std::string modelPath = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-m" || arg == "--mode") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -m/--mode требуется режим (normal/train/infer)"
                          << std::endl;
                return 1;
            }
            std::string val = argv[++i];
            if (val == "normal")
                mode = AnalysisMode::NORMAL;
            else if (val == "train")
                mode = AnalysisMode::TRAINING;
            else if (val == "infer")
                mode = AnalysisMode::INFERENCE;
            else {
                std::cerr << "Ошибка: режим должен быть normal, train или infer" << std::endl;
                return 1;
            }
        } else if (arg == "-l" || arg == "--label") {
            if (i + 1 >= argc) {
                std::cerr << "Ошибка: после -l/--label требуется значение (0 или 1)" << std::endl;
                return 1;
            }
            labelStr = argv[++i];
            if (labelStr != "0" && labelStr != "1") {
                std::cerr << "Ошибка: метка должна быть 0 (фон) или 1 (сигнал)" << std::endl;
                return 1;
            }
        } else if (arg[0] != '-') {
            // Позиционный аргумент — входной файл
            if (!inputRootFile.empty()) {
                std::cerr << "Ошибка: указан более одного входного файла" << std::endl;
                return 1;
            }
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
    fs::path mlOutputDir;
    if (mode == AnalysisMode::TRAINING) {
        mlOutputDir = fs::path(OUTPUT_ML_BASE_DIR) / "ml_training";
        fs::create_directories(mlOutputDir);
        std::cout << "Режим TRAINING: данные будут в " << fs::absolute(mlOutputDir) << std::endl;
    }

    fs::path processOutputDir = fs::path(OUTPUT_PDF_BASE_DIR) / processName;
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
    const std::string OUTPUT_COS_THETA_Z = makeOutputPath("cos_theta_Z_polar_angle");
    const std::string OUTPUT_DELTA_R = makeOutputPath("deltaR_jet1_jet2");
    const std::string OUTPUT_PHOTON_E_VS_RECOIL = makeOutputPath("photonE_vs_recoil_2d");
    const std::string OUTPUT_COS_THETA_JET = makeOutputPath("cosTheta_jets");
    const std::string OUTPUT_MET_PFO = makeOutputPath("MET_pfo");
    const std::string OUTPUT_MET_JET = makeOutputPath("MET_jets");
    const std::string OUTPUT_PMISS_MAG = makeOutputPath("Pmiss_magnitude");
    const std::string OUTPUT_COS_THETA_PMISS = makeOutputPath("cosTheta_Pmiss");

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
    std::vector<int> *particleType = nullptr;
    std::vector<double> *pfoE = nullptr;
    std::vector<double> *pfoPx = nullptr;
    std::vector<double> *pfoPy = nullptr;
    std::vector<double> *pfoPz = nullptr;

    if (APPLY_PRE_LEPTON_VETO || APPLY_PRE_HIGH_E_PHOTON_VETO || APPLY_PRE_ISOLATED_PHOTON_VETO) {
        tree->SetBranchAddress("particleType", &particleType);
        tree->SetBranchAddress("pfoE", &pfoE);
        tree->SetBranchAddress("pfoPx", &pfoPx);
        tree->SetBranchAddress("pfoPy", &pfoPy);
        tree->SetBranchAddress("pfoPz", &pfoPz);
    }

    // Гистограммы создаём только если не в режиме TRAINING
    TH1F *hInvMass = nullptr, *hRecoilMass = nullptr;
    TH1F *hCosThetaZ = nullptr, *hDeltaR = nullptr, *hCosThetaJet = nullptr;
    TH1F *hMETpfo = nullptr, *hMETjet = nullptr, *hPmissMag = nullptr, *hCosThetaPmiss = nullptr;
    TH2F *h2D_Mrecoil_vs_MET = nullptr, *h2D_Correlation = nullptr, *h2D_Mrecoil_vs_Pmiss = nullptr;
    TH2F *h2D_MET_vs_Pmiss = nullptr, *h2D_Mjj_vs_MET = nullptr, *h2D_Mjj_vs_Pmiss = nullptr;
    TH2F *h2D_CosThetaZ_vs_CosThetaPmiss = nullptr;

    // Создание гистограмм
    if (mode != AnalysisMode::TRAINING) {
        hInvMass = new TH1F("hInvMass", "Invariant Mass of Two Jets;M_{jj} [GeV];Events", MASS_BINS,
                            MASS_MIN_GEV, MASS_MAX_GEV);

        hRecoilMass =
            new TH1F("hRecoilMass", "Recoil Mass Against Two Jets;M_{recoil} [GeV];Events",
                     RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);

        h2D_Correlation = new TH2F(
            "h2D_Correlation", "Invariant Mass vs Recoil Mass;M_{jj} [GeV];M_{recoil} [GeV]",
            MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);

        hCosThetaZ =
            new TH1F("hCosThetaZ ", "cos#theta of Z Boson (Two-Jet System);cos#theta_{Z};Events ",
                     COS_THETA_Z_BINS, COS_THETA_Z_MIN, COS_THETA_Z_MAX);

        hDeltaR = new TH1F("hDeltaR",
                           "Distance #Delta R Between Two Jets;#Delta R = #sqrt{#Delta#eta^{2} + "
                           "#Delta#phi^{2}};Events",
                           DELTA_R_BINS, DELTA_R_MIN, DELTA_R_MAX);

        hCosThetaJet = new TH1F("hCosThetaJet", "cos#theta of Jets;cos#theta;Events",
                                COS_THETA_JET_BINS, COS_THETA_JET_MIN, COS_THETA_JET_MAX);

        hMETpfo = new TH1F("hMETpfo", "MET from all PFOs;MET_{PFO} [GeV];Events", MET_PFO_BINS,
                           MET_PFO_MIN, MET_PFO_MAX);
        hMETjet = new TH1F("hMETjet", "MET from Two Jets;MET_{jet} [GeV];Events", MET_JET_BINS,
                           MET_JET_MIN, MET_JET_MAX);
        hPmissMag = new TH1F("hPmissMag", "Magnitude of Missing 3-Momentum;|P_{miss}| [GeV];Events",
                             PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV);

        hCosThetaPmiss =
            new TH1F("hCosThetaPmiss", "cos#theta of Missing 3-Momentum;cos#theta_{miss};Events",
                     COS_THETA_PMISS_BINS, COS_THETA_PMISS_MIN, COS_THETA_PMISS_MAX);

        h2D_Mrecoil_vs_MET = new TH2F(
            "h2D_Mrecoil_vs_MET", "M_{recoil} vs MET_{jet};MET_{jet} [GeV];M_{recoil} [GeV]",
            MET_JET_BINS, MET_JET_MIN, MET_JET_MAX, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);
        h2D_Mrecoil_vs_Pmiss = new TH2F(
            "h2D_Mrecoil_vs_Pmiss", "M_{recoil} vs |P_{miss}|;|P_{miss}| [GeV];M_{recoil} [GeV]",
            PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);
        h2D_MET_vs_Pmiss = new TH2F(
            "h2D_MET_vs_Pmiss", "MET_{jet} vs |P_{miss}|;|P_{miss}| [GeV];MET_{jet} [GeV]",
            PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV, MET_JET_BINS, MET_JET_MIN, MET_JET_MAX);
        h2D_Mjj_vs_MET =
            new TH2F("h2D_Mjj_vs_MET", "M_{jj} vs MET_{jet};MET_{jet} [GeV];M_{jj} [GeV]",
                     MET_JET_BINS, MET_JET_MIN, MET_JET_MAX, MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV);
        h2D_Mjj_vs_Pmiss = new TH2F(
            "h2D_Mjj_vs_Pmiss", "M_{jj} vs |P_{miss}|;|P_{miss}| [GeV];M_{jj} [GeV]", PMISS_BINS,
            PMISS_MIN_GEV, PMISS_MAX_GEV, MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV);
        h2D_CosThetaZ_vs_CosThetaPmiss =
            new TH2F("h2D_CosThetaZ_vs_CosThetaPmiss",
                     "cos#theta_{Z} vs cos#theta_{miss};cos#theta_{miss};cos#theta_{Z}",
                     COS_THETA_PMISS_BINS, COS_THETA_PMISS_MIN, COS_THETA_PMISS_MAX,
                     COS_THETA_Z_BINS, COS_THETA_Z_MIN, COS_THETA_Z_MAX);
    }

    // Статистики
    CutStatistics stats;
    IsoElectronStats elecStats;

    // Структура фич
    MLFeatures trFeat;

    // Создаём TTree для фич
    TFile *trainFile = nullptr;
    TTree *trainTree = nullptr;
    if (mode == AnalysisMode::TRAINING) {
        std::string trainPath = (mlOutputDir / (processName + "_train.root")).string();
        trainFile = TFile::Open(trainPath.c_str(), "RECREATE");
        trainTree = new TTree("ml_features", "ML Training Features");
        trainTree->Branch("invMass", &trFeat.invMass, "invMass/D");
        trainTree->Branch("recoilMass", &trFeat.recoilMass, "recoilMass/D");
        trainTree->Branch("cosThetaZ", &trFeat.cosThetaZ, "cosThetaZ/D");
        trainTree->Branch("deltaR", &trFeat.deltaR, "deltaR/D");
        trainTree->Branch("ej1", &trFeat.ej1, "ej1/D");
        trainTree->Branch("ej2", &trFeat.ej2, "ej2/D");
        trainTree->Branch("eta_j1", &trFeat.eta_j1, "eta_j1/D");
        trainTree->Branch("eta_j2", &trFeat.eta_j2, "eta_j2/D");
        trainTree->Branch("pt_jj", &trFeat.pt_jj, "pt_jj/D");
        trainTree->Branch("met_pfo", &trFeat.met_pfo, "met_pfo/D");
        trainTree->Branch("pmag_miss", &trFeat.pmag_miss, "pmag_miss/D");
        trainTree->Branch("costh_miss", &trFeat.costh_miss, "costh_miss/D");
        trainTree->Branch("energy_asym", &trFeat.energy_asym, "energy_asym/D");
        trainTree->Branch("nconst_j1", &trFeat.nconst_j1, "nconst_j1/I");
        trainTree->Branch("nconst_j2", &trFeat.nconst_j2, "nconst_j2/I");
        trainTree->Branch("label", &trFeat.label, "label/I");
        trFeat.label = std::stoi(labelStr);
    }

    // Основной цикл по событиям
    Long64_t nEntries = tree->GetEntries();
    std::cout << "\nНачало обработки событий..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        logProgress(i + 1, nEntries, "Processing");
        stats.totalEvents++;

        // Сбор статистики по лептонам (до всех катов)
        if (particleType && pfoE && pfoPx && pfoPy && pfoPz) {
            for (size_t k = 0; k < particleType->size(); ++k) {

                // Отбираем электроны
                if (std::abs(particleType->at(k)) == PDG_ELECTRON) {
                    if (isLeptonIsolatedROOT_FSR(k, particleType, pfoE, pfoPx, pfoPy, pfoPz)) {
                        double px = pfoPx->at(k);
                        double py = pfoPy->at(k);
                        double pz = pfoPz->at(k);
                        double p = std::sqrt(px * px + py * py + pz * pz);
                        double cosTheta = (p > 1e-9) ? pz / p : 0.0;

                        elecStats.total++;
                        if (std::abs(cosTheta) < 0.7)
                            elecStats.barrel++;
                        else
                            elecStats.endcap++;
                    }
                }
            }
        }

        // === ПРЕДОТБОРЫ ===
        if (APPLY_PRE_LEPTON_VETO &&
            hasIsolatedLeptonROOT_FSR(particleType, pfoE, pfoPx, pfoPy, pfoPz))
            continue;
        stats.afterPreLeptonVeto++;

        if (APPLY_PRE_HIGH_E_PHOTON_VETO &&
            hasHighEnergyPhoton(particleType, pfoE, PHOTON_ENERGY_CUT_GEV))
            continue;
        stats.afterPreHighEPhotonVeto++;

        if (APPLY_PRE_ISOLATED_PHOTON_VETO &&
            hasIsolatedPhoton(particleType, pfoE, pfoPx, pfoPy, pfoPz, PHOTON_ISO_MIN_ENERGY_GEV,
                              PHOTON_ISO_COS_CONE_ANGLE, PHOTON_ISO_MAX_CONE_ENERGY_GEV))
            continue;
        stats.afterPreIsoPhotonVeto++;

        if (APPLY_PRE_TWO_JETS_REQUIREMENT && inclJetE->size() != 2)
            continue;
        stats.afterPreJetCount++;

        if (APPLY_PRE_CONSTITUENTS_REQUIREMENT) {
            int n1 = static_cast<int>(inclJetSize->at(0));
            int n2 = static_cast<int>(inclJetSize->at(1));
            if (n1 < MIN_CONSTITUENTS_PER_JET || n2 < MIN_CONSTITUENTS_PER_JET)
                continue;
        }
        stats.afterPreConstituents++;

        // ==================== КИНЕМАТИКА ====================
        TLorentzVector jet1(inclJetPx->at(0), inclJetPy->at(0), inclJetPz->at(0), inclJetE->at(0));
        TLorentzVector jet2(inclJetPx->at(1), inclJetPy->at(1), inclJetPz->at(1), inclJetE->at(1));
        TLorentzVector dijet = jet1 + jet2;

        double invMass = dijet.M();
        double recoilMass = calculateRecoilMass(dijet, SQRT_S_GEV);
        double cosThetaZ = std::cos(calculatePolarAngle(dijet));
        double deltaR = calculateDeltaR(jet1, jet2);
        double cosTheta1 = (jet1.P() > 1e-9) ? jet1.Pz() / jet1.P() : 0.0;
        double cosTheta2 = (jet2.P() > 1e-9) ? jet2.Pz() / jet2.P() : 0.0;
        double met_jet = dijet.Pt();

        double met_pfo = 0.0, pmiss_x = 0.0, pmiss_y = 0.0, pmiss_z = 0.0;
        if (pfoPx && pfoPy && pfoPz) {
            double sPx = 0, sPy = 0, sPz = 0;
            for (size_t k = 0; k < pfoPx->size(); ++k) {
                sPx += pfoPx->at(k);
                sPy += pfoPy->at(k);
                sPz += pfoPz->at(k);
            }
            met_pfo = std::sqrt(sPx * sPx + sPy * sPy);
            pmiss_x = -sPx;
            pmiss_y = -sPy;
            pmiss_z = -sPz;
        }
        double pmiss_mag = std::sqrt(pmiss_x * pmiss_x + pmiss_y * pmiss_y + pmiss_z * pmiss_z);
        double cosThetaPmiss =
            (pmiss_mag > 1e-9) ? std::max(-1.0, std::min(1.0, pmiss_z / pmiss_mag)) : 0.0;

        // === РЕЖИМ TRAINING: сохраняем фичи и переходим к следующему событию ===
        if (mode == AnalysisMode::TRAINING) {
            trFeat.invMass = invMass;
            trFeat.recoilMass = recoilMass;
            trFeat.cosThetaZ = cosThetaZ;
            trFeat.deltaR = deltaR;
            trFeat.ej1 = jet1.E();
            trFeat.ej2 = jet2.E();
            trFeat.eta_j1 = calculatePseudorapidity(jet1);
            trFeat.eta_j2 = calculatePseudorapidity(jet2);
            trFeat.pt_jj = met_jet;
            trFeat.met_pfo = met_pfo;
            trFeat.pmag_miss = pmiss_mag;
            trFeat.costh_miss = cosThetaPmiss;
            trFeat.energy_asym = (jet1.E() + jet2.E() > 1e-9)
                                     ? std::abs(jet1.E() - jet2.E()) / (jet1.E() + jet2.E())
                                     : 0.0;
            trFeat.nconst_j1 = static_cast<int>(inclJetSize->at(0));
            trFeat.nconst_j2 = static_cast<int>(inclJetSize->at(1));
            trainTree->Fill();
            continue; // не заполняем гистограммы в режиме TRAINING
        }

        // === РЕЖИМ INFERENCE: применяем ML вместо основных отборов ===
        if (mode == AnalysisMode::INFERENCE) {
            // Заглушка: замените на реальный вызов XGBoost
            // bool mlPass = predict_xgboost({invMass, recoilMass, ...}) >= 0.5;
            bool mlPass = true; // TODO: подключить модель
            if (!mlPass)
                continue;
            stats.finalSelected++;
        }

        // === ЗАПОЛНЕНИЕ ГИСТОГРАММ (только NORMAL и INFERENCE) ===
        if (mode != AnalysisMode::TRAINING) {
            hInvMass->Fill(invMass);
            hRecoilMass->Fill(recoilMass);
            h2D_Correlation->Fill(invMass, recoilMass);
            hCosThetaZ->Fill(cosThetaZ);
            hDeltaR->Fill(deltaR);
            hCosThetaJet->Fill(cosTheta1);
            hCosThetaJet->Fill(cosTheta2);
            hMETpfo->Fill(met_pfo);
            hMETjet->Fill(met_jet);
            hPmissMag->Fill(pmiss_mag);
            hCosThetaPmiss->Fill(cosThetaPmiss);

            h2D_Mrecoil_vs_MET->Fill(met_jet, recoilMass);
            h2D_Mrecoil_vs_Pmiss->Fill(pmiss_mag, recoilMass);
            h2D_MET_vs_Pmiss->Fill(pmiss_mag, met_jet);
            h2D_Mjj_vs_MET->Fill(met_jet, invMass);
            h2D_Mjj_vs_Pmiss->Fill(pmiss_mag, invMass);
            h2D_CosThetaZ_vs_CosThetaPmiss->Fill(cosThetaPmiss, cosThetaZ);
        }

        // === ОСНОВНЫЕ ОТБОРЫ (только NORMAL режим) ===
        if (mode == AnalysisMode::NORMAL) {
            if (APPLY_MAIN_MET_CUT && met_jet < MET_CUT_MIN_GEV)
                continue;
            stats.afterMetCut++;

            if (APPLY_MAIN_DIJET_MASS_WINDOW &&
                (invMass < DIJET_MASS_WINDOW_MIN_GEV || invMass > DIJET_MASS_WINDOW_MAX_GEV))
                continue;
            stats.afterDijetMassWindow++;

            if (APPLY_MAIN_COS_THETA_Z_CUT && std::abs(cosThetaZ) >= COS_THETA_Z_CUT)
                continue;
            stats.afterCosThetaZCut++;

            if (APPLY_MAIN_RECOIL_MASS_WINDOW && (recoilMass < RECOIL_MASS_WINDOW_MIN_GEV ||
                                                  recoilMass > RECOIL_MASS_WINDOW_MAX_GEV))
                continue;
            stats.afterRecoilMassWindow++;

            if (APPLY_MAIN_ELLIPSE_CUT &&
                !isInsideEllipse(invMass, recoilMass, ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV,
                                 ELLIPSE_B_GEV, ELLIPSE_THETA))
                continue;
            stats.afterEllipseCut++;
            stats.finalSelected++;
        }
    }

    // TRAINING: сохраняем TTree и выходим
    if (mode == AnalysisMode::TRAINING) {
        trainTree->Write();
        trainFile->Close();
        std::cout << "\nДанные для обучения сохранены: "
                  << (mlOutputDir / (processName + "_train.root")).string() << std::endl;
    }

    // NORMAL или INFERENCE
    else {
        // Итоговая статистика
        auto endTime = std::chrono::high_resolution_clock::now();
        auto totalSec =
            std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

        std::cout << "\n=======================================================" << std::endl;
        std::cout << "Обработка завершена!" << std::endl;
        std::cout << "Всего событий: " << nEntries << std::endl;
        std::cout << "Прошло времени: " << totalSec << " с (" << totalSec / 60.0 << " мин)"
                  << std::endl;

        stats.print(processName);
        elecStats.print();

        // Отрисовка гистограмм с условным отображением основных отборов
        std::vector<std::pair<double, std::string>> invMassMarks = {{MZ_GEV, "M_{Z}"}};
#if APPLY_MAIN_DIJET_MASS_WINDOW
        invMassMarks.emplace_back(DIJET_MASS_WINDOW_MIN_GEV, "M_{jj}^{min}");
        invMassMarks.emplace_back(DIJET_MASS_WINDOW_MAX_GEV, "M_{jj}^{max}");
#endif
        drawHistogram1D(hInvMass, "cInvMass", "M_{jj} [GeV]", OUTPUT_INV_MASS, invMassMarks, kRed,
                        2);

        std::vector<std::pair<double, std::string>> recoilMarks = {{MH_GEV, "M_{H}"}};
#if APPLY_MAIN_RECOIL_MASS_WINDOW
        recoilMarks.emplace_back(RECOIL_MASS_WINDOW_MIN_GEV, "M_{rec}^{min}");
        recoilMarks.emplace_back(RECOIL_MASS_WINDOW_MAX_GEV, "M_{rec}^{max}");
#endif
        drawHistogram1D(hRecoilMass, "cRecoilMass", "M_{recoil} [GeV]", OUTPUT_RECOIL_MASS,
                        recoilMarks, kBlue, 2);

        drawHistogram2D(h2D_Correlation, "c2D_Correlation", "M_{jj} [GeV]", "M_{recoil} [GeV]",
                        OUTPUT_2D_CORR, MZ_GEV, MH_GEV, "M_{Z}", "M_{H}",
#if APPLY_MAIN_ELLIPSE_CUT
                        ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV, ELLIPSE_B_GEV, ELLIPSE_THETA,
                        true
#else
                        -1, -1, -1, -1, 0, false
#endif
        );

        std::vector<std::pair<double, std::string>> cosThetaMarks;
#if APPLY_MAIN_COS_THETA_Z_CUT
        cosThetaMarks.emplace_back(COS_THETA_Z_CUT, "|cos#theta|^{cut}");
        cosThetaMarks.emplace_back(-COS_THETA_Z_CUT, "-|cos#theta|^{cut}");
#endif
        drawHistogram1D(hCosThetaZ, "cCosThetaZ", "cos#theta_{Z}", OUTPUT_COS_THETA_Z,
                        cosThetaMarks, kRed, 2);

        drawHistogram1D(hDeltaR, "cDeltaR", "#Delta R", OUTPUT_DELTA_R, {}, kMagenta, 2);

        drawHistogram1D(hCosThetaJet, "cCosThetaJet", "cos#theta", OUTPUT_COS_THETA_JET, {}, kCyan,
                        2);

        drawHistogram1D(hMETpfo, "cMETpfo", "MET_{PFO} [GeV]", OUTPUT_MET_PFO, {}, kOrange + 1, 2);

        std::vector<std::pair<double, std::string>> metMarks;
#if APPLY_MAIN_MET_CUT
        metMarks.emplace_back(MET_CUT_MIN_GEV, "MET_{min}");
#endif
        drawHistogram1D(hMETjet, "cMETjet", "MET_{jet} [GeV]", OUTPUT_MET_JET, metMarks, kViolet,
                        2);

        // h2D_Mrecoil_vs_MET: вертикальная линия MET > 20 GeV
        drawHistogram2D(h2D_Mrecoil_vs_MET, "c2D_Mrecoil_vs_MET", "MET_{jet} [GeV]",
                        "M_{recoil} [GeV]", makeOutputPath("2D_Mrecoil_vs_MET"),
#if APPLY_MAIN_MET_CUT
                        MET_CUT_MIN_GEV, -1, "MET_{min}", "");
#else
                        -1, -1, "", "");
#endif

        // h2D_Mrecoil_vs_Pmiss: без активных отборов по этим осям
        drawHistogram2D(h2D_Mrecoil_vs_Pmiss, "c2D_Mrecoil_vs_Pmiss", "|P_{miss}| [GeV]",
                        "M_{recoil} [GeV]", makeOutputPath("2D_Mrecoil_vs_Pmiss"));

        // h2D_MET_vs_Pmiss: вертикальная линия MET > 20 GeV (ось Y здесь MET)
        drawHistogram2D(h2D_MET_vs_Pmiss, "c2D_MET_vs_Pmiss", "|P_{miss}| [GeV]", "MET_{jet} [GeV]",
                        makeOutputPath("2D_MET_vs_Pmiss"),
#if APPLY_MAIN_MET_CUT
                        -1, MET_CUT_MIN_GEV, "", "MET_{min}");
#else
                        -1, -1, "", "");
#endif

        // h2D_Mjj_vs_MET: вертикальная линия MET > 20 GeV
        drawHistogram2D(h2D_Mjj_vs_MET, "c2D_Mjj_vs_MET", "MET_{jet} [GeV]", "M_{jj} [GeV]",
                        makeOutputPath("2D_Mjj_vs_MET"),
#if APPLY_MAIN_MET_CUT
                        MET_CUT_MIN_GEV, -1, "MET_{min}", "");
#else
                        -1, -1, "", "");
#endif

        // h2D_Mjj_vs_Pmiss: без активных отборов по этим осям
        drawHistogram2D(h2D_Mjj_vs_Pmiss, "c2D_Mjj_vs_Pmiss", "|P_{miss}| [GeV]", "M_{jj} [GeV]",
                        makeOutputPath("2D_Mjj_vs_Pmiss"));

        // h2D_CosThetaZ_vs_CosThetaPmiss: горизонтальные линии |cosθ_Z| < 0.98
        drawHistogram2D(h2D_CosThetaZ_vs_CosThetaPmiss, "c2D_CosThetaZ_vs_CosThetaPmiss",
                        "cos#theta_{miss}", "cos#theta_{Z}",
                        makeOutputPath("2D_CosThetaZ_vs_CosThetaPmiss"),
#if APPLY_MAIN_COS_THETA_Z_CUT
                        -1, COS_THETA_Z_CUT, "", "|cos#theta|^{cut}");
#else
                        -1, -1, "", "");
#endif

        // Остальные гистограммы без линий отборов
        drawHistogram1D(hPmissMag, "cPmissMag", "|P_{miss}| [GeV]", OUTPUT_PMISS_MAG, {}, kOrange,
                        2);
        drawHistogram1D(hCosThetaPmiss, "cCosThetaPmiss", "cos#theta_{miss}",
                        OUTPUT_COS_THETA_PMISS, {}, kMagenta, 2);

        // Очистка памяти
        delete hInvMass;
        delete hRecoilMass;
        delete h2D_Correlation;
        delete hCosThetaZ;
        delete hDeltaR;
        delete hCosThetaJet;
        delete hMETpfo;
        delete hMETjet;
        delete hPmissMag;
        delete hCosThetaPmiss;
        delete h2D_Mrecoil_vs_MET;
        delete h2D_Mrecoil_vs_Pmiss;
        delete h2D_MET_vs_Pmiss;
        delete h2D_Mjj_vs_MET;
        delete h2D_Mjj_vs_Pmiss;
        delete h2D_CosThetaZ_vs_CosThetaPmiss;
        inputFile->Close();
        delete inputFile;

        std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
                  << std::endl;
    }
    return 0;
}
