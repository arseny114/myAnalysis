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

#include <RooAddPdf.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooKeysPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace RooFit;

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

// Отрисовка стек гистограммы
void drawRecoilStack(const std::map<std::string, std::pair<TH1F *, ProcessInfo>> &processes,
                     const std::vector<std::string> &order, const std::string &outputFile) {

    TCanvas *c = new TCanvas("cRecoilStack", "Recoil Mass Stack", 1000, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetLogy(RECOIL_STACK_LOG_Y);

    THStack *stack =
        new THStack("recoilStack", "Recoil Mass Distribution;M_{recoil} [GeV];Expected events");
    TLegend *leg = new TLegend(0.75, 0.7, 0.95, 0.88);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);

    // Для отрисовки сигнала от низа стековой гистограммы
    TH1F *signalHistBottom = nullptr;

    // Вектор из временных копий гистограмм без весов для построения стека
    std::vector<TH1F *> tempHists;

    // Добавляем в стек строго в заданном порядке, но только те процессы, которые реально есть во
    // входных данных
    for (const auto &procName : order) {
        auto it = processes.find(procName);
        if (it != processes.end()) {
            TH1F *hist = it->second.first;
            const ProcessInfo &info = it->second.second;

            // Делаем так чтобы ошибки считались правильно и перевзвешиваем
            hist->Sumw2();
            if (info.legendName.find("signal") != std::string::npos)
                hist->Scale(info.weight * RECOIL_STACK_SIGNAL_MULTIPLIER);
            else
                hist->Scale(info.weight);

            // Создаем временную гистограмму для добавления в стек. Нам приходится делать такой
            // костыль потому что в стековой гистограмме ломается заливка, если исходная гистограмма
            // заполнена с весами. Специально не переносим ошибки из исходной гистограммы, потому
            // что это ломает заливку
            TH1F *hForStack = new TH1F(Form("hStack_%s", info.legendName.c_str()), hist->GetTitle(),
                                       hist->GetNbinsX(), hist->GetXaxis()->GetXmin(),
                                       hist->GetXaxis()->GetXmax());
            // Заполняем временную гистограмму
            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                hForStack->SetBinContent(bin, hist->GetBinContent(bin));
            }
            tempHists.push_back(hForStack); // Сохраняем чтобы потом почистить

            // Применяем стиль и добавляем в стек, сигнал не добавляем
            if (info.legendName.find("signal") == std::string::npos) {
                hForStack->SetFillColor(info.color);
                hForStack->SetFillStyle(info.fillStyle);
                hForStack->SetMarkerStyle(21);
                hForStack->SetMarkerColor(info.color);
                hForStack->SetLineWidth(1);
                stack->Add(hForStack);
                leg->AddEntry(hForStack, info.legendName.c_str(), "f");
            }

            // Клонируем сигнал чтобы отрисовать его от низа стековой гистограммы и настраиваем
            // параметры отображения
            if (info.legendName.find("signal") != std::string::npos) {
                signalHistBottom = (TH1F *)hForStack->Clone("signal_bottom");
                signalHistBottom->SetFillStyle(0);
                signalHistBottom->SetLineStyle(7);
                signalHistBottom->SetLineWidth(3);
                signalHistBottom->SetLineColor(info.color);
                leg->AddEntry(signalHistBottom, "Signal", "L");
            }
        }
    }

    // В случае логарифмической шкалы сами настраиваем границы по Y
    if (RECOIL_STACK_LOG_Y) {
        stack->SetMinimum(RECOIL_STACK_MIN_Y);
        stack->SetMaximum(RECOIL_STACK_MAX_Y);
    }

    stack->Draw();
    stack->GetXaxis()->SetTitle("M_{recoil} [GeV]");
    stack->GetYaxis()->SetTitle("Expected events after selection");
    stack->GetXaxis()->SetRangeUser(RECOIL_STACK_MIN_GEV, RECOIL_STACK_MAX_GEV);

    // Рисуем оверлей сигнала, если он есть во входных данных
    if (signalHistBottom) {
        signalHistBottom->Draw("HIST SAME");
    }

    leg->Draw();
    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;

    // Очистка
    delete leg;
    for (auto &h : tempHists)
        delete h;
    if (signalHistBottom)
        delete signalHistBottom;
    delete stack;
    delete c;
}

// Отрисовка сравнительной гистограммы массы отдачи (qqHX и qqHinvi)
// Строит график только если оба процесса присутствуют во входных данных
void drawRecoilComparison(const std::map<std::string, std::pair<TH1F *, ProcessInfo>> &processes,
                          const std::string &outputFile) {
    // Ищем процессы по именам, которые возвращает extractProcessName (без merged_ и .root)
    auto it_qqHX = processes.find("E240_qqHX");
    auto it_signal = processes.find("E240_qqHinvi");

    if (it_qqHX == processes.end() || it_signal == processes.end())
        return;

    // Клонируем гистограммы, чтобы не повлиять на стек или другие отрисовки
    TH1F *hBg = (TH1F *)it_qqHX->second.first->Clone("hCompBg");
    TH1F *hSig = (TH1F *)it_signal->second.first->Clone("hCompSig");

    // Применяем кастомные веса
    hBg->Sumw2();
    hBg->Scale(RECOIL_COMP_W_QQHX);

    hSig->Sumw2();
    hSig->Scale(RECOIL_COMP_W_SIGNAL);

    TCanvas *c = new TCanvas("cRecoilComp", "Recoil Mass Comparison", 900, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetLogy(RECOIL_COMP_LOG_Y);

    // Настройки осей
    if (RECOIL_COMP_LOG_Y) {
        hBg->SetMinimum(RECOIL_STACK_MIN_Y);
        hBg->SetMaximum(RECOIL_STACK_MAX_Y);
    }
    hBg->GetXaxis()->SetTitle("M_{recoil} [GeV]");
    hBg->GetYaxis()->SetTitle("Expected events after selection");
    hBg->GetXaxis()->SetTitleSize(0.045);
    hBg->GetXaxis()->SetTitleOffset(1.1);
    hBg->GetYaxis()->SetTitleSize(0.045);
    hBg->GetYaxis()->SetTitleOffset(1.1);
    hBg->SetStats(0); // Отключаем стандартное окно статистики

    // Стиль фона (qqHX)
    hBg->SetFillColor(kYellow);
    hBg->SetLineColor(kBlack);
    hBg->SetLineWidth(2);
    hBg->Draw("HIST E");

    // Стиль сигнала (qqHinvi)
    hSig->SetLineColor(kRed + 1);
    hSig->SetLineWidth(3);
    hSig->SetLineStyle(kSolid);
    hSig->Draw("HIST E SAME");

    // Легенда
    TLegend *leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->AddEntry(hBg, "qqHX (Signal + Background)", "F");
    leg->AddEntry(hSig, "qqHinvi (Signal)", "L");
    leg->Draw();

    c->SaveAs(outputFile.c_str());

    // Очистка
    delete c;
    delete leg;
    delete hBg;
    delete hSig;
}

// Отрисовка гистограмм предотборов (линии разных цветов, режим same)
void drawPreselectionHistograms(
    const std::map<std::string, std::pair<TH1F *, ProcessInfo>> &processes,
    const std::string &histName, const std::string &title, const std::string &xTitle,
    const std::string &outputFile, double markValue = -1) {
    TCanvas *c = new TCanvas("cPreselection", title.c_str(), 900, 700);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetLogy(true);

    // Отключаем статистику в углу
    gStyle->SetOptStat(0);

    TLegend *leg = new TLegend(0.74, 0.75, 0.98, 0.98);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.025);

    TH1F *firstHist = nullptr;

    // Находим первую гистограмму для настройки осей
    for (const auto &proc : processes) {
        if (proc.second.first->GetEntries() > 0) {
            firstHist = proc.second.first;
            break;
        }
    }

    if (!firstHist) {
        std::cerr << "[drawPreselection] Нет данных для отрисовки!\n";
        delete c;
        delete leg;
        return;
    }

    // Рисуем первую гистограмму
    bool first = true;
    for (const auto &procEntry : processes) {
        TH1F *hist = procEntry.second.first;
        const ProcessInfo &info = procEntry.second.second;

        if (hist->GetEntries() == 0)
            continue;

        // Нормализуем на число событий
        TH1F *hNorm = (TH1F *)hist->Clone(Form("hNorm_%s", procEntry.first.c_str()));
        double scale = (hist->Integral() > 0) ? 1.0 / hist->Integral() : 1.0;
        hNorm->Scale(scale);

        hNorm->SetLineColor(info.color);
        hNorm->SetLineWidth(3);
        hNorm->SetFillColor(0);
        hNorm->SetMarkerStyle(20);
        hNorm->SetMarkerColor(info.color);
        hNorm->SetMarkerSize(0.8);

        if (first) {
            hNorm->Draw("HIST");
            hNorm->GetXaxis()->SetTitle(xTitle.c_str());
            hNorm->GetYaxis()->SetTitle("Normalized Events");
            hNorm->GetYaxis()->SetTitleSize(0.045);
            hNorm->GetXaxis()->SetTitleSize(0.045);
            first = false;
        } else {
            hNorm->Draw("HIST SAME");
        }

        leg->AddEntry(hNorm, info.legendName.c_str(), "L");
    }

    // Рисуем линию отбора
    if (markValue > 0) {
        c->Update();
        // В лог-режиме GetUymin/GetUymax возвращают log10 значений
        double yminLog = gPad->GetUymin();
        double ymaxLog = gPad->GetUymax();
        double yReal_min = std::pow(10.0, yminLog);
        double yReal_max = std::pow(10.0, ymaxLog);

        // Линия на весь диапазон по Y
        TLine *line = new TLine(markValue, yReal_min, markValue, yReal_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(6);
        line->SetLineStyle(kDashed);
        line->Draw();

        // Стрелка: начинается на 70% высоты, заканчивается на 15%
        // В log-шкале позиции вычисляем через экспоненту
        double arrowY2 = std::pow(10.0, yminLog + (ymaxLog - yminLog) * 0.15);
        double arrowY1 = std::pow(10.0, yminLog + (ymaxLog - yminLog) * 0.70);
        TArrow *arrow = new TArrow(markValue, arrowY1, markValue, arrowY2, 0.030, "|>");
        arrow->SetLineColor(kRed);
        arrow->SetLineWidth(3);
        arrow->SetFillColor(kRed);
        arrow->Draw();
    }

    leg->Draw();
    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;

    delete c;
    delete leg;
}

// Отрисовка суммарных гистограмм основных отборов
// Поддерживает до 2 стрелок (для окон типа [min, max])
void drawMainSelectionHistograms(
    const std::map<std::string, std::pair<TH1F *, ProcessInfo>> &processes,
    const std::string &title, const std::string &xTitle, const std::string &outputFile,
    std::vector<double> markValues = {}) {
    TCanvas *c = new TCanvas("cMainSel", title.c_str(), 900, 700);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetLogy(true);
    gStyle->SetOptStat(0);

    TLegend *leg = new TLegend(0.74, 0.75, 0.98, 0.98);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.025);

    bool first = true;
    for (const auto &procEntry : processes) {
        TH1F *hist = procEntry.second.first;
        const ProcessInfo &info = procEntry.second.second;
        if (hist->GetEntries() == 0)
            continue;

        TH1F *hNorm = (TH1F *)hist->Clone(Form("hNormMain_%s", procEntry.first.c_str()));
        double scale = (hist->Integral() > 0) ? 1.0 / hist->Integral() : 1.0;
        hNorm->Scale(scale);
        hNorm->SetLineColor(info.color);
        hNorm->SetLineWidth(3);
        hNorm->SetFillColor(0);
        hNorm->SetMarkerSize(0);

        if (first) {
            hNorm->Draw("HIST");
            hNorm->GetXaxis()->SetTitle(xTitle.c_str());
            hNorm->GetYaxis()->SetTitle("Normalized Events");
            hNorm->GetYaxis()->SetTitleSize(0.045);
            hNorm->GetXaxis()->SetTitleSize(0.045);
            first = false;
        } else {
            hNorm->Draw("HIST SAME");
        }
        leg->AddEntry(hNorm, info.legendName.c_str(), "L");
    }

    if (!first && !markValues.empty()) {
        c->Update();
        double yminLog = gPad->GetUymin();
        double ymaxLog = gPad->GetUymax();

        for (double val : markValues) {
            double yReal_min = std::pow(10.0, yminLog);
            double yReal_max = std::pow(10.0, ymaxLog);
            TLine *line = new TLine(val, yReal_min, val, yReal_max);
            line->SetLineColor(kRed);
            line->SetLineWidth(6);
            line->SetLineStyle(kDashed);
            line->Draw();

            double arrowY2 = std::pow(10.0, yminLog + (ymaxLog - yminLog) * 0.15);
            double arrowY1 = std::pow(10.0, yminLog + (ymaxLog - yminLog) * 0.70);
            TArrow *arrow = new TArrow(val, arrowY1, val, arrowY2, 0.030, "|>");
            arrow->SetLineColor(kRed);
            arrow->SetLineWidth(3);
            arrow->SetFillColor(kRed);
            arrow->Draw();
        }
    }

    leg->Draw();
    c->SaveAs(outputFile.c_str());
    std::cout << "Сохранено: " << outputFile << std::endl;
    delete c;
    delete leg;
}

// =============================================================================
// ФУНКЦИЯ runMrecoilTemplateFit
// =============================================================================
//
// Выполняет шаблонный фит распределения M_recoil методом расширенного
// максимума правдоподобия (Extended Maximum Likelihood, EML).
//
// ВХОДНЫЕ ДАННЫЕ:
//   vSignal    - MC-события процесса qqHinvi , пары (M_recoil, вес).
//   vBkg       - MC-события всех фоновых процессов, пары (M_recoil, вес).
//   outputPath - путь для сохранения PDF с графиком.
//
// АЛГОРИТМ:
//   1. Строим два шаблона (PDF) через ядерную оценку плотности RooKeysPdf:
//        fs(x) - форма распределения сигнала
//        fb(x) - форма распределения фона
//      Шаблоны нормированы на 1 (это просто формы, не нормировки).
//
//   2. Строим псевдоданные. Это смесь фона и сигнала с долей mu:
//        dsData = vBkg (с весами) + vSignal (с весами × mu)
//      Суммарный взвешенный интеграл dsData = sumW_B + mu * sumW_S.
//
//   3. Фит 1 (нулевая гипотеза H0): фитируем dsData только фоном.
//      Свободный параметр: nB (число событий фона).
//      nS зафиксирован = 0. Получаем NLL_b (значение -log(L)).
//
//   4. Фит 2 (альтернативная гипотеза H1): фитируем dsData суммой fs+fb.
//      Свободные параметры: nS и nB. Получаем NLL_sb.
//
//   5. Вычисляем статистику:
//      - Простая оценка: Z_approx = nS / sqrt(nS + nB)
//      - Отношение правдоподобий: q0 = -2*ln(L_b/L_sb) = 2*(NLL_b - NLL_sb)
//        По теореме Вилкса при большой статистике q0 ~ χ2(1),
//        значимость Z_LRT = sqrt(q0) в единицах σ.
//
// =============================================================================
void runMrecoilTemplateFit(const std::vector<std::pair<double, double>> &vSignal,
                           const std::vector<std::pair<double, double>> &vBkg,
                           const std::vector<std::pair<double, double>> &v_qqHX,
                           const std::string &outputPath) {

    if (vSignal.empty() || vBkg.empty()) {
        std::cerr << "[Fit] Ошибка: входные данные пусты!\n";
        return;
    }

    std::cout << "\n[Fit] =====================================================\n";
    std::cout << "[Fit] Запуск шаблонного фита M_recoil\n";
    std::cout << "[Fit] Сигнальных MC-событий (qqHinvi): " << vSignal.size() << "\n";
    std::cout << "[Fit] Фоновых MC-событий (остальные):  " << vBkg.size() << "\n";
    std::cout << "[Fit] qqHX MC-событий (сигнал+фон):    " << v_qqHX.size() << "\n";
    std::cout << "[Fit] mu (доля сигнала):               " << FIT_PSEUDO_MU << "\n";

    // =========================================================================
    // ШАГ 1: НАБЛЮДАЕМАЯ ПЕРЕМЕННАЯ
    // =========================================================================
    //
    // RooRealVar описывает физическую переменную, по которой делается фит.
    // Третий и четвёртый аргументы это диапазон допустимых значений.
    // События вне диапазона [fitMin, fitMax] будут игнорироваться.
    RooRealVar Mrecoil("Mrecoil", "M_{recoil} [GeV]", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);

    // Именованный диапазон fitRange нужен, чтобы явно передавать его в
    // fitTo() и plotOn().
    Mrecoil.setRange("fitRange", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);

    // =========================================================================
    // ШАГ 2: ШАБЛОННЫЕ ДАТАСЕТЫ (для построения PDF)
    // =========================================================================
    //
    // RooDataSet это контейнер событий для RooFit.
    // WeightVar(wVar) сообщает RooFit, что wVar хранит статистический вес
    // каждого события. Внутри RooFit это влияет на нормировку и ошибки.
    //
    // wVar должна входить в RooArgSet, переданный конструктору датасета,
    // иначе RooFit не будет знать о ней при чтении датасета.
    //
    // Диапазон wVar выбран с запасом.
    RooRealVar wVar("eventWeight", "Event weight", 1e-9, 1e9);
    RooArgSet argSet(Mrecoil, wVar);

    // --- Шаблон сигнала ---
    // dsSignal используется только для построения формы PDF сигнала fs(x).
    RooDataSet *dsSignal =
        new RooDataSet("dsSignal", "Signal template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vSignal) {
        // entry.first  - значение M_recoil для данного события
        // entry.second - вес события
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsSignal->add(argSet, entry.second);
    }

    // --- Шаблон фона ---
    // dsBkg используется только для построения формы PDF фона fb(x).
    RooDataSet *dsBkg =
        new RooDataSet("dsBkg", "Background template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        // entry.first  - значение M_recoil для данного события
        // entry.second - вес события
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second);
    }

    // Добавляем qqHX с положительным весом (содержит сигнал + фон от инклюзивных распадов H)
    for (const auto &entry : v_qqHX) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second); // добавляем с +весом
    }

    // Вычитаем чистый сигнал qqHinvi с отрицательным весом (чтобы убрать сигнал из qqHX)
    for (const auto &entry : vSignal) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, -entry.second); // вычитаем с -весом
    }

    // Суммарные взвешенные числа событий в шаблонах.
    // sumEntries() с весами возвращает сумму весов всех событий.
    // Это ожидаемое число событий при данной светимости.
    double sumW_S = dsSignal->sumEntries();
    double sumW_B = dsBkg->sumEntries();

    std::cout << "[Fit] sumW_S (ожидаемый сигнал, полный):  " << sumW_S << "\n";
    std::cout << "[Fit] sumW_B (ожидаемый фон с учётом qqHX - qqHinvi):  " << sumW_B << "\n";
    std::cout << "[Fit] Инжектируем mu * sumW_S =           " << FIT_PSEUDO_MU * sumW_S
              << " событий сигнала в псевдоданные\n";

    if (sumW_S <= 0 || sumW_B <= 0) {
        std::cerr << "[Fit] Ошибка: суммарные веса равны нулю. Проверь входные данные.\n";
        delete dsSignal;
        delete dsBkg;
        return;
    }

    // =========================================================================
    // ШАГ 3: ПОСТРОЕНИЕ ШАБЛОНОВ (RooKeysPdf)
    // =========================================================================
    //
    // RooKeysPdf реализует ядерную оценку плотности (Kernel Density Estimation).
    // Идея: вместо гистограммы каждое событие заменяется гауссовым ядром,
    // а PDF это сумма всех ядер, нормированная на 1.
    //
    // Параметр adaptivity (ρ) управляет шириной ядра:
    //   h_i = h_0 * (f(x_i))^(-1/2) * rho
    // где f(x_i) - локальная плотность событий. В плотных областях ядро
    // уже (не размывает пики), в разреженных шире (не даёт шума).
    //
    // MirrorBoth/NoMirror нужен для обработки граничных эффектов:
    //   NoMirror: ядро может "вытечь" за край диапазона, что привидет к занижению у границ.
    //   MirrorBoth: добавляет зеркальные копии событий, что корректирует края.
    auto mirrorOpt = FIT_KEYSPDF_MIRROR ? RooKeysPdf::MirrorBoth : RooKeysPdf::NoMirror;

    RooKeysPdf pdfSignal("pdfSignal", "Signal PDF (fs)", Mrecoil, *dsSignal, mirrorOpt,
                         FIT_KEYSPDF_ADAPTIVITY_SIGNAL);
    RooKeysPdf pdfBkg("pdfBkg", "Background PDF (fb)", Mrecoil, *dsBkg, mirrorOpt,
                      FIT_KEYSPDF_ADAPTIVITY_BGD);

    // Фиксируем диапазон нормировки PDF явно.
    pdfSignal.setNormRange("fitRange");
    pdfBkg.setNormRange("fitRange");

    // =========================================================================
    // ШАГ 4: ПРОВЕРКА КАЧЕСТВА ШАБЛОНОВ (CHI2 TEST)
    // =========================================================================
    //
    // С помощью каждого шаблона генерируем toy MC выборку большого размера,
    // строим гистограммы и сравниваем с исходными данными через Chi2Test.
    // Это позволяет оценить, насколько хорошо RooKeysPdf описывает MC.
    //
    std::cout << "\n[Fit] =====================================================\n";
    std::cout << "[Fit] Проверка качества шаблонов (Chi2 test)\n";

    // Переменные для хранения chi2 значений
    double chi2B = 0.0;
    double chi2S = 0.0;

    // --- Проверка шаблона фона с построением графика ---
    {
        const int nBinsTempl = TEMPLATE_BACKGROUND_BINS;

        // Создаём график для визуализации
        RooPlot *frameBkg = Mrecoil.frame(RooFit::Title("Background Template Quality Check"));
        dsBkg->plotOn(frameBkg, RooFit::Binning(nBinsTempl), RooFit::MarkerStyle(21),
                      RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
        pdfBkg.plotOn(frameBkg, RooFit::LineColor(kBlue), RooFit::LineWidth(2));

        // Вычисляем Chi2/NDF
        chi2B = frameBkg->chiSquare();

        // Сохраняем график
        TCanvas *cBkg = new TCanvas("cBkg", "Background Template", 800, 600);
        frameBkg->Draw();

        // Добавляем текст с Chi2 на график
        TLatex latexB;
        latexB.SetNDC();
        latexB.SetTextFont(42);
        latexB.SetTextSize(0.04);
        latexB.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.3f", chi2B));

        cBkg->SaveAs((outputPath + "/template_background.pdf").c_str());
        std::cout << "[Fit] График сохранен: " << outputPath << "/template_background.pdf\n";

        delete cBkg;
        delete frameBkg;
    }

    // --- Проверка шаблона сигнала с построением графика ---
    {
        const int nBinsTempl = TEMPLATE_SIGNAL_BINS;

        // Создаём график для визуализации
        RooPlot *frameSig = Mrecoil.frame(RooFit::Title("Signal Template Quality Check"));
        dsSignal->plotOn(frameSig, RooFit::Binning(nBinsTempl), RooFit::MarkerStyle(22),
                         RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
        pdfSignal.plotOn(frameSig, RooFit::LineColor(kRed), RooFit::LineWidth(2));

        // Вычисляем Chi2/NDF
        chi2S = frameSig->chiSquare();

        // Сохраняем график
        TCanvas *cSig = new TCanvas("cSig", "Signal Template", 800, 600);
        frameSig->Draw();

        // Добавляем текст с Chi2 на график
        TLatex latexS;
        latexS.SetNDC();
        latexS.SetTextFont(42);
        latexS.SetTextSize(0.04);
        latexS.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.3f", chi2S));

        cSig->SaveAs((outputPath + "/template_signal.pdf").c_str());
        std::cout << "[Fit] График сохранен: " << outputPath << "/template_signal.pdf\n";

        delete cSig;
        delete frameSig;
    }

    // =========================================================================
    // ШАГ 5: ПАРАМЕТРЫ НОРМИРОВКИ (число событий)
    // =========================================================================
    //
    // Начальные значения = ожидаемые числа событий из MC.
    // Для границ не допускаем отрицательных значений. Верхний предел берем с
    // запасом, чтобы фит не упирался в границу.
    double expectedNS = FIT_PSEUDO_MU * sumW_S;
    double expectedNB = sumW_B;

    // Верхняя граница nS: берём максимум из (10 * ожидаемого) и абсолютного минимума 100,
    // чтобы фит не был зажат при малом ожидаемом сигнале.
    double nS_max = std::max(10.0 * expectedNS, 100.0);
    double nB_max = 3.0 * expectedNB;

    RooRealVar nS("nS", "Signal yield", expectedNS, 0.0, nS_max);
    RooRealVar nB("nB", "Background yield", expectedNB, 0.0, nB_max);

    // =========================================================================
    // ШАГ 6: КОМБИНИРОВАННАЯ МОДЕЛЬ
    // =========================================================================
    //
    // RooAddPdf с двумя аргументами-нормировками создаёт расширенный PDF:
    //   model(x) = nS * fs(x) + nB * fb(x)
    //
    RooAddPdf model("model", "nS*fs + nB*fb", RooArgList(pdfSignal, pdfBkg), RooArgList(nS, nB));
    model.setNormRange("fitRange");

    // =========================================================================
    // ШАГ 7: ПСЕВДОДАННЫЕ
    // =========================================================================
    //
    // Псевдоданные (toy MC) - это то, что в реальном эксперименте было бы
    // экспериментальными данными. Здесь мы формируем их из MC, инжектируя
    // известную долю сигнала mu.
    //
    // dsData = все фоновые MC события (с весами) + сигнальные MC (с весами × mu)
    //
    // Важно: dsData это взвешенный датасет (MC с весами != 1). Поэтому при фите
    // используем SumW2Error(kTRUE). Это корректирует ковариационную матрицу
    // по формуле:
    //   V_corr = H^{-1} * (sum w_i^2 * grad_i * grad_i^T) * H^{-1}
    // вместо стандартной V = H^{-1}, где H это матрица Гессе NLL.
    // Без SumW2Error ошибки параметров будут занижены.
    //
    RooDataSet *dsData =
        new RooDataSet("dsData", "Pseudo-data (Bkg + mu*Signal)", argSet, RooFit::WeightVar(wVar));

    // Добавляем фоновые события с их оригинальными весами
    for (const auto &entry : vBkg) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsData->add(argSet, entry.second);
    }

    // Добавляем сигнальные события с весами, умноженными на mu.
    if (FIT_PSEUDO_MU > 0.0) {
        for (const auto &entry : vSignal) {
            if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
                continue;
            Mrecoil.setVal(entry.first);
            dsData->add(argSet, FIT_PSEUDO_MU * entry.second);
        }
    }

    double sumW_Data = dsData->sumEntries();
    std::cout << "[Fit] Псевдоданные: sumW = " << sumW_Data << "  (ожидалось "
              << expectedNB + expectedNS << ")\n";

    // =========================================================================
    // ШАГ 8: ФИТ 1. НУЛЕВАЯ ГИПОТЕЗА H0 (только фон, nS = 0)
    // =========================================================================
    //
    // H0: сигнала нет. Единственный свободный параметр nB.
    // nS фиксируем в 0 через setConstant(true).
    //
    // fitTo() выполняет минимизацию -log(L) по свободным параметрам.
    // Опции:
    //   Extended(kTRUE)    — использовать расширенный NLL (с пуассоновским членом)
    //   Range("fitRange")  — фитировать только в указанном диапазоне
    //   SumW2Error(kTRUE)  — корректировать ошибки для взвешенных данных
    //   Save(kTRUE)        — сохранить результат в RooFitResult
    //   PrintLevel(-1)     — подавить вывод MINUIT (0 = минимальный, -1 = тихий)
    //
    nS.setVal(0.0);
    nS.setConstant(kTRUE); // заморозить nS = 0

    // Сбрасываем nB к ожидаемому значению перед каждым фитом,
    // чтобы минимизатор не стартовал из плохой точки.
    nB.setVal(expectedNB);

    RooFitResult *fitResB =
        model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                    RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1));

    if (!fitResB || fitResB->status() != 0) {
        // status() == 0 означает, что MINUIT успешно сошёлся.
        // Ненулевой статус: 1 = covariance forced positive definite,
        //                   2 = HESSE failed, 3 = MINOS failed.
        std::cerr << "[Fit] Предупреждение: фит H0 завершился со статусом "
                  << (fitResB ? fitResB->status() : -1) << "\n";
    }
    double nll_b = fitResB ? fitResB->minNll() : 0.0;
    double fit_nB_H0 = nB.getVal();

    std::cout << "[Fit] H0 (только фон): nB = " << fit_nB_H0 << ",  NLL_b = " << nll_b << "\n";

    // =========================================================================
    // ШАГ 9: ФИТ 2. АЛЬТЕРНАТИВНАЯ ГИПОТЕЗА H1 (сигнал + фон)
    // =========================================================================
    //
    // H1: сигнал присутствует. Оба параметра nS и nB свободны.
    // Получаем NLL_sb - минимум NLL при наилучшем описании данных.
    //
    nS.setConstant(kFALSE); // размораживаем nS
    nS.setVal(expectedNS);  // стартуем из физически разумной точки
    nB.setVal(expectedNB);

    RooFitResult *fitResSB =
        model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                    RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1));

    if (!fitResSB || fitResSB->status() != 0) {
        std::cerr << "[Fit] Предупреждение: фит H1 завершился со статусом "
                  << (fitResSB ? fitResSB->status() : -1) << "\n";
    }
    double nll_sb = fitResSB ? fitResSB->minNll() : 0.0;

    // Извлекаем итоговые значения параметров из результата фита H1.
    double fit_nS = nS.getVal();
    double fit_nS_err = nS.getError();
    double fit_nB = nB.getVal();
    double fit_nB_err = nB.getError();

    // =========================================================================
    // ШАГ 10: ОЦЕНКА ЗНАЧИМОСТИ
    // =========================================================================
    //
    // --- Простая оценка ---
    // Z_approx = nS / sqrt(nS + nB)
    // Это приближение работает при большой статистике и гауссовых ошибках.
    // Смысл: число сигнальных событий в единицах статистической флуктуации фона.
    //
    double significance_approx = 0.0;
    if (fit_nS > 0 && (fit_nS + fit_nB) > 0)
        significance_approx = fit_nS / std::sqrt(fit_nS + fit_nB);

    // --- Оценка через отношение правдоподобий (Likelihood Ratio Test) ---
    // Статистика теста: q0 = -2 * ln(L_b / L_sb) = 2 * (NLL_b - NLL_sb)
    //
    // По теореме Вилкса, при H0 верной и большой статистике q0 распределена
    // как χ2(1) (хи-квадрат с одной степенью свободы) (так как H1 имеет на
    // один параметр больше, чем H0).
    //
    // Значимость в σ: Z_LRT = sqrt(q0)
    //
    double q0 = 2.0 * (nll_b - nll_sb);
    double pValue = TMath::Prob(q0, 1);
    double significance_lrt = std::sqrt(std::max(0.0, q0));

    std::cout << "\n[Fit] =====================================================\n";
    std::cout << "[Fit] РЕЗУЛЬТАТЫ ФИТА (mu = " << FIT_PSEUDO_MU << ")\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   nS (fitted)  = " << fit_nS << " ± " << fit_nS_err << "\n";
    std::cout << "[Fit]   nB (fitted)  = " << fit_nB << " ± " << fit_nB_err << "\n";
    std::cout << "[Fit]   nS (expected)= " << expectedNS << "\n";
    std::cout << "[Fit]   nB (expected)= " << expectedNB << "\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   Фон (Bkg):  chi2 / NDF = " << chi2B << "\n";
    std::cout << "[Fit]   Сигнал:     chi2 / NDF = " << chi2S << "\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   NLL (H0, фон):      " << nll_b << "\n";
    std::cout << "[Fit]   NLL (H1, сигн+фон): " << nll_sb << "\n";
    std::cout << "[Fit]   q0 = 2*(NLL_b - NLL_sb) = " << q0 << "\n";
    std::cout << "[Fit]   p-value (H0): = " << pValue << "\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   Значимость (простая):  Z = " << significance_approx << " σ\n";
    std::cout << "[Fit]   Значимость (LRT):      Z = " << significance_lrt << " σ\n";
    std::cout << "[Fit] =====================================================\n\n";

    // =========================================================================
    // ШАГ 11: ВИЗУАЛИЗАЦИЯ
    // =========================================================================
    //
    // RooPlot это объект для отрисовки данных и PDF поверх одной оси.
    // Важно: порядок plotOn() имеет значение, каждый вызов добавляет слой.
    // Данные рисуем первыми, потом кривые, так данные не перекрываются.
    //
    // Нормировка кривых: при plotOn() RooFit автоматически масштабирует PDF
    // на (число событий в данных) × (ширина бина), если не задано иное.
    // В расширенном фите используется Normalization(sumW_Data, RooAbsReal::NumEvent).
    //
    TCanvas *cFit = new TCanvas("cMrecoilFit", "M_{recoil} Template Fit", 900, 700);
    cFit->SetLeftMargin(0.13);
    cFit->SetRightMargin(0.05);
    cFit->SetBottomMargin(0.12);

    RooPlot *frame = Mrecoil.frame(RooFit::Title("M_{recoil} Template Fit"));
    frame->GetXaxis()->SetTitle("M_{recoil} [GeV]");
    frame->GetYaxis()->SetTitle("Events");

    // Данные: DataError(RooAbsData::SumW2) рисует ошибки как sqrt(sum w_i^2),
    // что правильно для взвешенных MC (а не просто sqrt(N)).
    dsData->plotOn(frame, RooFit::Binning(FIT_PLOT_BINS), RooFit::MarkerStyle(20),
                   RooFit::MarkerSize(1.0), RooFit::DataError(RooAbsData::SumW2),
                   RooFit::Name("data"));

    // Суммарный фит (H1): нормируем на суммарный вес псевдоданных
    model.plotOn(frame, RooFit::Normalization(sumW_Data, RooAbsReal::NumEvent),
                 RooFit::LineColor(kGreen + 1), RooFit::LineWidth(2), RooFit::Name("total"));

    // Фоновая компонента: масштабируем на fit_nB (число событий фона из фита)
    model.plotOn(frame, RooFit::Components(pdfBkg),
                 RooFit::Normalization(fit_nB, RooAbsReal::NumEvent), RooFit::LineColor(kRed),
                 RooFit::LineStyle(kDashed), RooFit::LineWidth(2), RooFit::Name("bkg"));

    // Сигнальная компонента: масштабируем на fit_nS
    model.plotOn(frame, RooFit::Components(pdfSignal),
                 RooFit::Normalization(fit_nS, RooAbsReal::NumEvent), RooFit::LineColor(kBlue),
                 RooFit::LineStyle(kDashed), RooFit::LineWidth(2), RooFit::Name("sig"));

    // Задаём диапазон Y ДО Draw(), но с учётом режима шкалы. В лог-режиме минимум обязательно > 0.
    if (FIT_PLOT_LOG_Y) {
        frame->GetYaxis()->SetRangeUser(0.1, FIT_PLOT_YMAX);
    } else {
        frame->GetYaxis()->SetRangeUser(0.0, FIT_PLOT_YMAX);
    }

    // Вертикальные пунктирные линии по границам фитового диапазона.
    TLine *lineMin = new TLine(FIT_MRECOIL_MIN, 0, FIT_MRECOIL_MIN, FIT_PLOT_YMAX * 0.95);
    TLine *lineMax = new TLine(FIT_MRECOIL_MAX, 0, FIT_MRECOIL_MAX, FIT_PLOT_YMAX * 0.95);
    lineMin->SetLineStyle(kDotted);
    lineMin->SetLineColor(kGray + 1);
    lineMax->SetLineStyle(kDotted);
    lineMax->SetLineColor(kGray + 1);

    frame->Draw();
    if (FIT_PLOT_LOG_Y)
        cFit->SetLogy();
    lineMin->Draw("SAME");
    lineMax->Draw("SAME");

    // Легенда
    TLegend *leg = new TLegend(0.76, 0.85, 0.98, 0.98);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.023);
    leg->AddEntry(frame->findObject("data"), "Pseudo-Data", "p");
    leg->AddEntry(frame->findObject("total"), "Total Fit", "l");
    leg->AddEntry(frame->findObject("bkg"), "Background Template", "l");
    leg->AddEntry(frame->findObject("sig"), "Signal Template", "l");
    leg->Draw();

    // Информационный блок с результатами фита
    TPaveText *info = new TPaveText(0.76, 0.6, 0.98, 0.83, "NDC NB");
    info->SetFillColor(0);
    info->SetFillStyle(1001);
    info->SetBorderSize(1);
    info->SetTextAlign(12);
    info->SetTextSize(0.020);
    info->AddText(Form("mu = %.2f", FIT_PSEUDO_MU));
    info->AddText(Form("nS = %.1f #pm %.1f", fit_nS, fit_nS_err));
    info->AddText(Form("nB = %.1f #pm %.1f", fit_nB, fit_nB_err));
    info->AddText(Form("#chi^{2}_{bkg} / NDF = %.3f", chi2B));
    info->AddText(Form("#chi^{2}_{sig} / NDF = %.3f", chi2S));
    info->AddText(Form("nB = %.1f #pm %.1f", fit_nB, fit_nB_err));
    info->AddText(Form("Z (simple) = %.2f #sigma", significance_approx));
    info->AddText(Form("Z (LRT)    = %.2f #sigma", significance_lrt));
    info->AddText(Form("p-value    = %.3g", pValue));
    info->Draw();

    cFit->SaveAs((outputPath + "/template_fit_recoil.pdf").c_str());
    std::cout << "[Fit] График сохранен: " << outputPath << "/template_fit_recoil.pdf\n";

    // =========================================================================
    // ОЧИСТКА ПАМЯТИ
    // =========================================================================
    delete leg;
    delete info;
    delete lineMin;
    delete lineMax;
    delete frame;
    delete cFit;
    delete dsData;
    delete dsSignal;
    delete dsBkg;
    if (fitResB)
        delete fitResB;
    if (fitResSB)
        delete fitResSB;
}

// =============================================================================
// ФУНКЦИЯ runMrecoilAnalyticalFit (RooFit: RooKeysPdf + Chebyshev)
// =============================================================================
//
// Выполняет аналитический фит распределения M_recoil методом расширенного
// максимума правдоподобия (Extended Maximum Likelihood, EML):
//
//   - Сигнал: RooKeysPdf (непараметрическая ядерная оценка плотности из MC).
//   - Фон: полином Чебышёва 4-й степени (RooChebychev).
//
void runMrecoilAnalyticalFit(const std::vector<std::pair<double, double>> &vSignal,
                             const std::vector<std::pair<double, double>> &vBkg,
                             const std::vector<std::pair<double, double>> &v_qqHX,
                             const std::string &outputPath) {
    if (vSignal.empty() || vBkg.empty()) {
        std::cerr << "[AnalyticalFit] Ошибка: входные данные пусты!\n";
        return;
    }

    std::cout << "\n[AnalyticalFit] =====================================================\n";
    std::cout << "[AnalyticalFit] Запуск аналитического фита M_recoil (RooFit)\n";
    std::cout << "[AnalyticalFit] Сигнал: RooKeysPdf\n";
    std::cout << "[AnalyticalFit] Фон:    Chebyshev(4)\n";

    // =========================================================================
    // ШАГ 1: НАБЛЮДАЕМАЯ ПЕРЕМЕННАЯ
    // =========================================================================
    RooRealVar Mrecoil("Mrecoil", "M_{recoil} [GeV]", ANALYTICAL_FIT_MRECOIL_MIN,
                       ANALYTICAL_FIT_MRECOIL_MAX);
    Mrecoil.setRange("fitRange", ANALYTICAL_FIT_MRECOIL_MIN, ANALYTICAL_FIT_MRECOIL_MAX);

    // =========================================================================
    // ШАГ 2: СОЗДАНИЕ ДАТАСЕТОВ
    // =========================================================================
    RooRealVar wVar("eventWeight", "Event weight", 1e-9, 1e9);
    RooArgSet argSet(Mrecoil, wVar);

    // --- Шаблон сигнала ---
    RooDataSet *dsSignalTemplate =
        new RooDataSet("dsSignalTemplate", "Signal template", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vSignal) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsSignalTemplate->add(argSet, entry.second);
    }

    // --- Шаблон фона (vBkg + v_qqHX − vSignal) ---
    RooDataSet *dsBkgTemplate =
        new RooDataSet("dsBkgTemplate", "Background template", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, entry.second);
    }
    for (const auto &entry : v_qqHX) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, entry.second);
    }
    for (const auto &entry : vSignal) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, -entry.second);
    }

    // --- Псевдо-данные для общего фита (фон + mu * сигнал) ---
    RooDataSet *dsData =
        new RooDataSet("dsData", "Pseudo-data (Bkg + mu*Signal)", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsData->add(argSet, entry.second);
    }
    for (const auto &entry : v_qqHX) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsData->add(argSet, entry.second);
    }
    for (const auto &entry : vSignal) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        // Вычитаем сигнал из фона (он туда попал через vBkg/v_qqHX) и добавляем mu*сигнал
        dsData->add(argSet, -entry.second);
        dsData->add(argSet, ANALYTICAL_FIT_PSEUDO_MU * entry.second);
    }

    double sumW_Data = dsData->sumEntries();
    std::cout << "[AnalyticalFit] Суммарный вес данных: " << sumW_Data << "\n";
    std::cout << "[AnalyticalFit] Событий в сигнальном шаблоне: " << dsSignalTemplate->sumEntries()
              << "\n";
    std::cout << "[AnalyticalFit] Событий в фоновом шаблоне:   " << dsBkgTemplate->sumEntries()
              << "\n";

    // =========================================================================
    // ШАГ 3: ПОСТРОЕНИЕ ШАБЛОНА СИГНАЛА (RooKeysPdf)
    // =========================================================================
    std::cout << "\n[AnalyticalFit] Построение шаблона сигнала (RooKeysPdf)...\n";

    auto mirrorOpt = FIT_KEYSPDF_MIRROR ? RooKeysPdf::MirrorBoth : RooKeysPdf::NoMirror;
    RooKeysPdf pdfSigInit("pdfSigInit", "Signal PDF (RooKeysPdf)", Mrecoil, *dsSignalTemplate,
                          mirrorOpt, FIT_KEYSPDF_ADAPTIVITY_SIGNAL);
    pdfSigInit.setNormRange("fitRange");

    if (dsSignalTemplate->numEntries() > 0) {
        // --- График шаблона сигнала с χ²/ndf ---
        TCanvas *cSigInit = new TCanvas("cSigInit", "Signal Template Check", 800, 600);
        cSigInit->SetLeftMargin(0.13);
        cSigInit->SetRightMargin(0.05);
        cSigInit->SetBottomMargin(0.12);
        cSigInit->SetTopMargin(0.08);

        RooPlot *frameSigInit = Mrecoil.frame(RooFit::Title("Signal Template (RooKeysPdf)"));
        dsSignalTemplate->plotOn(frameSigInit, RooFit::Binning(ANALYTICAL_FIT_PLOT_BINS),
                                 RooFit::MarkerStyle(20), RooFit::MarkerSize(1.0),
                                 RooFit::DataError(RooAbsData::SumW2));
        pdfSigInit.plotOn(frameSigInit, RooFit::LineColor(kBlue), RooFit::LineWidth(2));
        frameSigInit->Draw();

        double chi2_sig = frameSigInit->chiSquare();

        TPaveText *ptSig = new TPaveText(0.2, 0.8, 0.45, 0.92, "NDC NB");
        ptSig->SetFillColor(0);
        ptSig->SetBorderSize(1);
        ptSig->SetTextSize(0.020);
        ptSig->AddText("Signal Template (RooKeysPdf)");
        ptSig->AddText(Form("Adaptivity = %.2f", FIT_KEYSPDF_ADAPTIVITY_SIGNAL));
        ptSig->AddText(Form("#chi^{2}/ndf = %.2f", chi2_sig));
        ptSig->Draw();

        cSigInit->SaveAs((outputPath + "/analytical_fit_signal_only.pdf").c_str());
        delete ptSig;
        delete frameSigInit;
        delete cSigInit;
    }

    // =========================================================================
    // ШАГ 4: ФИТ ЧИСТОГО ФОНА (Только Chebyshev)
    // =========================================================================
    std::cout << "\n[AnalyticalFit] Фит чистого фона (Chebyshev)...\n";

    RooRealVar c0_init("c0_init", "Chebyshev coeff 0", 0.0, -1.0, 1.0);
    RooRealVar c1_init("c1_init", "Chebyshev coeff 1", 0.0, -1.0, 1.0);
    RooRealVar c2_init("c2_init", "Chebyshev coeff 2", 0.0, -1.0, 1.0);
    RooRealVar c3_init("c3_init", "Chebyshev coeff 3", 0.0, -1.0, 1.0);

    RooChebychev pdfBkgInit("pdfBkgInit", "Chebyshev Background", Mrecoil,
                            RooArgList(c0_init, c1_init, c2_init, c3_init));

    RooFitResult *fitResBkgInit =
        pdfBkgInit.fitTo(*dsBkgTemplate, RooFit::Range("fitRange"), RooFit::SumW2Error(kTRUE),
                         RooFit::Strategy(1), RooFit::PrintLevel(-1), RooFit::Save(true));

    std::cout << "[AnalyticalFit] Параметры фона (чистый фит):\n";
    std::cout << "  c0 = " << c0_init.getVal() << " ± " << c0_init.getError() << "\n";
    std::cout << "  c1 = " << c1_init.getVal() << " ± " << c1_init.getError() << "\n";
    std::cout << "  c2 = " << c2_init.getVal() << " ± " << c2_init.getError() << "\n";
    std::cout << "  c3 = " << c3_init.getVal() << " ± " << c3_init.getError() << "\n";

    // --- График чистого фона с χ²/ndf ---
    TCanvas *cBkgInit = new TCanvas("cBkgInit", "Background Template Fit", 800, 600);
    cBkgInit->SetLeftMargin(0.13);
    cBkgInit->SetRightMargin(0.05);
    cBkgInit->SetBottomMargin(0.12);
    cBkgInit->SetTopMargin(0.08);

    RooPlot *frameBkgInit = Mrecoil.frame(RooFit::Title("Background Template Fit (Chebyshev)"));
    dsBkgTemplate->plotOn(frameBkgInit, RooFit::Binning(ANALYTICAL_FIT_PLOT_BINS),
                          RooFit::MarkerStyle(20), RooFit::MarkerSize(1.0),
                          RooFit::DataError(RooAbsData::SumW2), RooFit::Name("dsBkgTemplateData"));
    pdfBkgInit.plotOn(frameBkgInit, RooFit::LineColor(kRed), RooFit::LineWidth(2),
                      RooFit::Name("pdfBkgInitFit"));
    frameBkgInit->Draw();

    double chi2_bkg = frameBkgInit->chiSquare();

    TLegend *legBkg = new TLegend(0.76, 0.8, 0.98, 0.92);
    legBkg->SetBorderSize(1);
    legBkg->SetFillColor(0);
    legBkg->SetTextSize(0.023);
    legBkg->AddEntry(frameBkgInit->findObject("dsBkgTemplateData"), "Background Data", "lep");
    legBkg->AddEntry(frameBkgInit->findObject("pdfBkgInitFit"), "Chebyshev Fit", "l");
    legBkg->Draw();

    TPaveText *ptBkg = new TPaveText(0.2, 0.8, 0.45, 0.92, "NDC NB");
    ptBkg->SetFillColor(0);
    ptBkg->SetBorderSize(1);
    ptBkg->SetTextSize(0.020);
    ptBkg->AddText("Background Fit (Chebyshev)");
    ptBkg->AddText(Form("c_{0} = %.3f #pm %.3f", c0_init.getVal(), c0_init.getError()));
    ptBkg->AddText(Form("c_{1} = %.3f #pm %.3f", c1_init.getVal(), c1_init.getError()));
    ptBkg->AddText(Form("c_{2} = %.3f #pm %.3f", c2_init.getVal(), c2_init.getError()));
    ptBkg->AddText(Form("c_{3} = %.3f #pm %.3f", c3_init.getVal(), c3_init.getError()));
    ptBkg->AddText(Form("#chi^{2}/ndf = %.2f", chi2_bkg));
    ptBkg->Draw();

    cBkgInit->SaveAs((outputPath + "/analytical_fit_background_only.pdf").c_str());
    delete legBkg;
    delete ptBkg;
    delete frameBkgInit;
    delete cBkgInit;

    // =========================================================================
    // ШАГ 5: ПАРАМЕТРЫ НОРМИРОВКИ
    // =========================================================================

    // Вычисляем ожидаемые значения фона и сигнала
    double expectedNS = ANALYTICAL_FIT_PSEUDO_MU * dsSignalTemplate->sumEntries();
    double expectedNB = sumW_Data - expectedNS;

    // Создаем свободные параметры числа фона и сигнала (фит будет их варьировать)
    RooRealVar nS("nS", "Signal yield",
                  expectedNS,                                // начальное значение
                  0.0,                                       // minValue
                  10.0 * std::max(expectedNS, 1.0) + 100.0); // maxValue
    RooRealVar nB("nB", "Background yield", expectedNB, 0.0,
                  3.0 * std::max(expectedNB, 1.0) + 1000.0);

    // =========================================================================
    // ШАГ 6: МОДЕЛЬ СИГНАЛА (RooKeysPdf)
    // Форма сигнала уже зафиксирована самим фактом построения RooKeysPdf из MC шаблона.
    // =========================================================================
    RooKeysPdf pdfSig("pdfSig", "Signal PDF (RooKeysPdf)", Mrecoil, *dsSignalTemplate, mirrorOpt,
                      FIT_KEYSPDF_ADAPTIVITY_SIGNAL);
    pdfSig.setNormRange("fitRange");

    // =========================================================================
    // ШАГ 7: МОДЕЛЬ ФОНА (Только Chebyshev)
    // =========================================================================
    RooRealVar c0("c0", "Chebyshev coeff 0", c0_init.getVal(), -1.0, 1.0);
    RooRealVar c1("c1", "Chebyshev coeff 1", c1_init.getVal(), -1.0, 1.0);
    RooRealVar c2("c2", "Chebyshev coeff 2", c2_init.getVal(), -1.0, 1.0);
    RooRealVar c3("c3", "Chebyshev coeff 3", c3_init.getVal(), -1.0, 1.0);

    RooChebychev pdfBkg("pdfBkg", "Chebyshev Background", Mrecoil, RooArgList(c0, c1, c2, c3));

    // =========================================================================
    // ШАГ 8: КОМБИНИРОВАННАЯ МОДЕЛЬ (extended)
    // =========================================================================
    RooAddPdf model("model", "Signal + Background", RooArgList(pdfSig, pdfBkg), RooArgList(nS, nB));

    // =========================================================================
    // ШАГ 9: ФИТ H0 (только фон, nS = 0)
    // =========================================================================
    std::cout << "\n[AnalyticalFit] Фит H0: только фон...\n";

    // Фиксируем число сигнала равным 0, выставляем самое близкое к правильному начальное значение
    // для фона
    nS.setVal(0.0);
    nS.setConstant(kTRUE);
    nB.setVal(expectedNB);

    RooFitResult *fitResB = model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                                        RooFit::SumW2Error(kTRUE), RooFit::Strategy(1),
                                        RooFit::Minimizer("Minuit2", "migrad"),
                                        RooFit::PrintLevel(-1), RooFit::Save(kTRUE));

    // Сохраняем параметры фона из фита H0
    double c0_H0 = c0.getVal();
    double c1_H0 = c1.getVal();
    double c2_H0 = c2.getVal();
    double c3_H0 = c3.getVal();
    double nll_b = fitResB ? fitResB->minNll() : 0.0;
    double fit_nB_H0 = nB.getVal();

    if (!fitResB || fitResB->status() != 0)
        std::cerr << "[AnalyticalFit] Предупреждение: фит H0 статус "
                  << (fitResB ? fitResB->status() : -1) << "\n";

    std::cout << "[AnalyticalFit] H0: nB = " << fit_nB_H0 << ", NLL = " << nll_b << "\n";

    // =========================================================================
    // ШАГ 10: ФИТ H1 (сигнал + фон)
    // =========================================================================
    std::cout << "[AnalyticalFit] Фит H1: сигнал + фон...\n";

    // Стартуем из параметров H0 для обеспечения вложенности
    c0.setVal(c0_H0);
    c1.setVal(c1_H0);
    c2.setVal(c2_H0);
    c3.setVal(c3_H0);
    nB.setVal(fit_nB_H0);

    nS.setConstant(kFALSE);
    nS.setVal(expectedNS);
    nB.setVal(expectedNB);

    RooFitResult *fitResSB = model.fitTo(
        *dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"), RooFit::SumW2Error(kTRUE),
        RooFit::Strategy(1), RooFit::Minimizer("Minuit2", "migrad"), RooFit::PrintLevel(-1),
        RooFit::Save(kTRUE));

    // Сохраняем параметры фона из фита H1
    double c0_H1 = c0.getVal();
    double c1_H1 = c1.getVal();
    double c2_H1 = c2.getVal();
    double c3_H1 = c3.getVal();
    double nll_sb = fitResSB ? fitResSB->minNll() : 0.0;
    double fit_nS = nS.getVal(), fit_nS_err = nS.getError();
    double fit_nB = nB.getVal(), fit_nB_err = nB.getError();

    if (!fitResSB || fitResSB->status() != 0)
        std::cerr << "[AnalyticalFit] Предупреждение: фит H1 статус "
                  << (fitResSB ? fitResSB->status() : -1) << "\n";

    std::cout << "[AnalyticalFit] H1: nS = " << fit_nS << ", nB = " << fit_nB
              << ", NLL = " << nll_sb << "\n";

    // =========================================================================
    // ШАГ 10b: ВАЛИДАЦИЯ: ФИТ С ЗАФИКСИРОВАННОЙ ФОРМОЙ ФОНА
    // =========================================================================
    // Цель: Проверить, что если форма фона (Чебышев) зафиксирована на "истинных"
    // значениях из MC, а форма сигнала зафиксирована (RooKeysPdf), то свободные
    // параметры нормировки (nS, nB) сойдутся к ожидаемым инжектированным значениям.
    std::cout << "\n[AnalyticalFit] =====================================================\n";
    std::cout << "[AnalyticalFit] ВАЛИДАЦИЯ: Фит с зафиксированной формой фона\n";

    // 1. Фиксируем коэффициенты Чебышева к значениям, полученным из чистого фона (Шаг 4)
    // Это наши истинные параметры формы фона.
    c0.setVal(c0_init.getVal());
    c1.setVal(c1_init.getVal());
    c2.setVal(c2_init.getVal());
    c3.setVal(c3_init.getVal());

    c0.setConstant(kTRUE);
    c1.setConstant(kTRUE);
    c2.setConstant(kTRUE);
    c3.setConstant(kTRUE);

    // 2. Параметры числа событий (nS и nB) оставляем свободными,
    // но задаем им хорошие начальные значения (ожидаемые)
    nS.setVal(expectedNS);
    nB.setVal(expectedNB);
    nS.setConstant(kFALSE);
    nB.setConstant(kFALSE);

    // 3. Запускаем фит
    RooFitResult *fitResFixed = model.fitTo(
        *dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"), RooFit::SumW2Error(kTRUE),
        RooFit::Strategy(1), RooFit::Minimizer("Minuit2", "migrad"), RooFit::PrintLevel(1),
        RooFit::Save(kTRUE));

    if (!fitResFixed || fitResFixed->status() != 0) {
        std::cerr << "[AnalyticalFit] Предупреждение: валидационный фит завершился со статусом "
                  << (fitResFixed ? fitResFixed->status() : -1) << "\n";
    }

    double fit_nS_fixed = nS.getVal();
    double fit_nS_fixed_err = nS.getError();
    double fit_nB_fixed = nB.getVal();
    double fit_nB_fixed_err = nB.getError();
    double nll_fixed = fitResFixed ? fitResFixed->minNll() : 0.0;

    std::cout << "[AnalyticalFit] -----------------------------------------------------\n";
    std::cout << "[AnalyticalFit] Результаты фита с фиксированной формой:\n";
    std::cout << "[AnalyticalFit]   nS (fitted)  = " << std::fixed << std::setprecision(2)
              << fit_nS_fixed << " ± " << fit_nS_fixed_err << "  (ожидалось: " << expectedNS
              << ")\n";
    std::cout << "[AnalyticalFit]   nB (fitted)  = " << std::fixed << std::setprecision(2)
              << fit_nB_fixed << " ± " << fit_nB_fixed_err << "  (ожидалось: " << expectedNB
              << ")\n";
    std::cout << "[AnalyticalFit]   NLL          = " << nll_fixed << "\n";
    std::cout << "[AnalyticalFit] =====================================================\n\n";

    // 4. Освобождаем параметры фона, чтобы не сломать последующую визуализацию (Шаг 12),
    // которая ожидает, что параметры можно менять.
    c0.setConstant(kFALSE);
    c1.setConstant(kFALSE);
    c2.setConstant(kFALSE);
    c3.setConstant(kFALSE);

    if (fitResFixed)
        delete fitResFixed;

    // =========================================================================
    // ШАГ 11: ОЦЕНКА ЗНАЧИМОСТИ
    // =========================================================================
    double q0 = 2.0 * (nll_b - nll_sb);
    if (q0 < 0.0)
        q0 = 0.0;

    double pValue = TMath::Prob(q0, 1);
    double significance_lrt = std::sqrt(q0);
    double significance_simple =
        (fit_nS > 0 && (fit_nS + fit_nB) > 0) ? fit_nS / std::sqrt(fit_nS + fit_nB) : 0.0;

    std::cout << "\n[AnalyticalFit] =====================================================\n";
    std::cout << "[AnalyticalFit] РЕЗУЛЬТАТЫ ФИТА\n";
    std::cout << "[AnalyticalFit] -----------------------------------------------------\n";
    std::cout << "[AnalyticalFit] Сигнал (RooKeysPdf):\n";
    std::cout << "[AnalyticalFit]   nS      = " << fit_nS << " ± " << fit_nS_err << "\n";
    std::cout << "[AnalyticalFit]   expnS   = " << expectedNS << "\n";
    std::cout << "[AnalyticalFit] Фон (Chebyshev):\n";
    std::cout << "[AnalyticalFit]   c0    = " << c0.getVal() << " ± " << c0.getError() << "\n";
    std::cout << "[AnalyticalFit]   c1    = " << c1.getVal() << " ± " << c1.getError() << "\n";
    std::cout << "[AnalyticalFit]   c2    = " << c2.getVal() << " ± " << c2.getError() << "\n";
    std::cout << "[AnalyticalFit]   c3    = " << c3.getVal() << " ± " << c3.getError() << "\n";
    std::cout << "[AnalyticalFit]   nB    = " << fit_nB << " ± " << fit_nB_err << "\n";
    std::cout << "[AnalyticalFit]   expnB = " << expectedNB << "\n";
    std::cout << "[AnalyticalFit] -----------------------------------------------------\n";
    std::cout << "[AnalyticalFit] NLL (H0)  = " << nll_b << "\n";
    std::cout << "[AnalyticalFit] NLL (H1)  = " << nll_sb << "\n";
    std::cout << "[AnalyticalFit] q0        = " << q0 << "\n";
    std::cout << "[AnalyticalFit] p-value   = " << pValue << "\n";
    std::cout << "[AnalyticalFit] -----------------------------------------------------\n";
    std::cout << "[AnalyticalFit] Significance (simple): Z = " << significance_simple << " σ\n";
    std::cout << "[AnalyticalFit] Significance (LRT):    Z = " << significance_lrt << " σ\n";
    std::cout << "[AnalyticalFit] =====================================================\n\n";

    // =========================================================================
    // ШАГ 12: ВИЗУАЛИЗАЦИЯ ОБЩЕГО ФИТА (с гипотезой H0)
    // =========================================================================
    TCanvas *cFit = new TCanvas("cAnalyticalFit", "M_{recoil} Analytical Fit", 900, 700);
    cFit->SetLeftMargin(0.13);
    cFit->SetRightMargin(0.05);
    cFit->SetBottomMargin(0.12);
    cFit->SetTopMargin(0.08);

    RooPlot *frame = Mrecoil.frame(RooFit::Title("Total Data (Bkg + Signal)"));

    // 1. Данные
    dsData->plotOn(frame, RooFit::Binning(ANALYTICAL_FIT_PLOT_BINS), RooFit::MarkerStyle(20),
                   RooFit::MarkerSize(1.0), RooFit::DataError(RooAbsData::SumW2),
                   RooFit::Name("data"));

    // 2. H1: полный фит (сигнал + фон) - зелёная сплошная линия
    model.plotOn(frame, RooFit::LineColor(kGreen + 2), RooFit::LineWidth(2),
                 RooFit::Name("H1_total"));

    // 3. H0: только фон.
    //
    // H0: восстанавливаем параметры из фита H0
    c0.setVal(c0_H0);
    c1.setVal(c1_H0);
    c2.setVal(c2_H0);
    c3.setVal(c3_H0);
    nB.setVal(fit_nB_H0);

    // Не меняем nS и nB у основного model, а создаём отдельный PDF
    RooAddPdf modelH0("modelH0", "Background only (H0)", RooArgList(pdfBkg), RooArgList(nB));
    modelH0.plotOn(frame, RooFit::LineColor(kBlack), RooFit::LineStyle(kDotted),
                   RooFit::LineWidth(2), RooFit::Name("H0_only_bkg"));

    // 4. Компоненты H1 (сигнал и фон из полного фита)
    // Восстанавливаем параметры H1 для отрисовки компонент
    c0.setVal(c0_H1);
    c1.setVal(c1_H1);
    c2.setVal(c2_H1);
    c3.setVal(c3_H1);
    nB.setVal(fit_nB);

    model.plotOn(frame, RooFit::Components(pdfBkg), RooFit::LineColor(kRed),
                 RooFit::LineStyle(kDashed), RooFit::LineWidth(2), RooFit::Name("H1_bkg"));
    model.plotOn(frame, RooFit::Components(pdfSig), RooFit::LineColor(kBlue),
                 RooFit::LineStyle(kDashed), RooFit::LineWidth(2), RooFit::Name("H1_sig"));

    frame->Draw();

    double yMax = frame->GetYaxis()->GetXmax() * 0.95;
    TLine *lineMin = new TLine(ANALYTICAL_FIT_MRECOIL_MIN, 0, ANALYTICAL_FIT_MRECOIL_MIN, yMax);
    TLine *lineMax = new TLine(ANALYTICAL_FIT_MRECOIL_MAX, 0, ANALYTICAL_FIT_MRECOIL_MAX, yMax);
    lineMin->SetLineStyle(kDotted);
    lineMin->SetLineColor(kGray + 1);
    lineMax->SetLineStyle(kDotted);
    lineMax->SetLineColor(kGray + 1);
    lineMin->Draw("SAME");
    lineMax->Draw("SAME");

    TLegend *leg = new TLegend(0.7, 0.85, 0.98, 0.98);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.020);
    leg->AddEntry(frame->findObject("data"), "Pseudo-Data", "lep");
    leg->AddEntry(frame->findObject("H1_total"), "H1: Signal + Background", "l");
    leg->AddEntry(frame->findObject("H0_only_bkg"), "H0: Background only", "l");
    leg->AddEntry(frame->findObject("H1_sig"), "Signal component (H1)", "l");
    leg->AddEntry(frame->findObject("H1_bkg"), "Background component (H1)", "l");
    leg->Draw();

    TPaveText *info = new TPaveText(0.15, 0.68, 0.35, 0.92, "NDC NB");
    info->SetFillColor(0);
    info->SetFillStyle(1001);
    info->SetBorderSize(1);
    info->SetTextAlign(12);
    info->SetTextSize(0.020);
    info->AddText("Total Fit Results");
    info->AddText("Signal: RooKeysPdf");
    info->AddText(Form("nS = %.1f #pm %.1f", fit_nS, fit_nS_err));
    info->AddText(Form("nB = %.1f #pm %.1f", fit_nB, fit_nB_err));

    // Вычисляем chi2 явно для H0 и H1 и данных data
    info->AddText(Form("H0: #chi^{2}/ndf = %.2f", frame->chiSquare("H0_only_bkg", "data")));
    info->AddText(Form("H1: #chi^{2}/ndf = %.2f", frame->chiSquare("H1_total", "data")));
    info->AddText(Form("Z (LRT)    = %.2f #sigma", significance_lrt));
    info->AddText(Form("p-value    = %.3g", pValue));
    info->Draw();

    cFit->SaveAs((outputPath + "/analytical_fit_total.pdf").c_str());
    std::cout << "[AnalyticalFit] График сохранен: " << outputPath << "/analytical_fit_total.pdf\n";

    // =========================================================================
    // ОЧИСТКА ПАМЯТИ
    // =========================================================================
    delete leg;
    delete info;
    delete lineMin;
    delete lineMax;
    delete frame;
    delete cFit;
    delete dsData;
    delete dsSignalTemplate;
    delete dsBkgTemplate;
    if (fitResBkgInit)
        delete fitResBkgInit;
    if (fitResB)
        delete fitResB;
    if (fitResSB)
        delete fitResSB;
}

// =============================================================================
// ФУНКЦИЯ runMrecoilAnalyticalScanMu
// =============================================================================
//
// Выполняет сканирование по параметру mu для аналитического фита (RooKeysPdf + Chebyshev)
// с целью построения зависимости Z_LRT от mu и поиска mu для достижения 5 сигм.
void runMrecoilAnalyticalScanMu(const std::vector<std::pair<double, double>> &vSignal,
                                const std::vector<std::pair<double, double>> &vBkg,
                                const std::vector<std::pair<double, double>> &v_qqHX,
                                const std::string &outputPath) {
    std::cout << "\n[AnalyticalScan] ====================================================\n";
    std::cout << "[AnalyticalScan] Запуск сканирования по mu (Аналитический фит)\n";
    std::cout << "[AnalyticalScan] Диапазон mu: [" << ANALYTICAL_SCAN_MU_MIN << ", "
              << ANALYTICAL_SCAN_MU_MAX << "]\n";
    std::cout << "[AnalyticalScan] Шаг mu: " << ANALYTICAL_SCAN_MU_STEP << "\n";
    std::cout << "[AnalyticalScan] Целевая значимость: " << ANALYTICAL_SCAN_Z_TARGET << " σ\n";

    // 1. Подготовка наблюдаемой переменной и аргументов
    RooRealVar Mrecoil("Mrecoil", "M_{recoil} [GeV]", ANALYTICAL_FIT_MRECOIL_MIN,
                       ANALYTICAL_FIT_MRECOIL_MAX);
    Mrecoil.setRange("fitRange", ANALYTICAL_FIT_MRECOIL_MIN, ANALYTICAL_FIT_MRECOIL_MAX);
    RooRealVar wVar("eventWeight", "Event weight", 1e-9, 1e9);
    RooArgSet argSet(Mrecoil, wVar);

    // 2. Создание шаблонов
    RooDataSet *dsSignalTemplate =
        new RooDataSet("dsSignalTemplate", "Signal template", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vSignal) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsSignalTemplate->add(argSet, entry.second);
    }

    RooDataSet *dsBkgTemplate =
        new RooDataSet("dsBkgTemplate", "Background template", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, entry.second);
    }
    for (const auto &entry : v_qqHX) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, entry.second);
    }
    for (const auto &entry : vSignal) {
        if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN || entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkgTemplate->add(argSet, -entry.second); // Вычитаем сигнал из qqHX
    }

    // 3. Построение PDF
    auto mirrorOpt = FIT_KEYSPDF_MIRROR ? RooKeysPdf::MirrorBoth : RooKeysPdf::NoMirror;
    RooKeysPdf pdfSig("pdfSig_scan", "Signal PDF", Mrecoil, *dsSignalTemplate, mirrorOpt,
                      FIT_KEYSPDF_ADAPTIVITY_SIGNAL);
    pdfSig.setNormRange("fitRange");

    // Предварительный фит чистого фона для получения хороших стартовых значений Чебышева
    RooRealVar c0_init("c0_init", "c0", 0.0, -1.0, 1.0);
    RooRealVar c1_init("c1_init", "c1", 0.0, -1.0, 1.0);
    RooRealVar c2_init("c2_init", "c2", 0.0, -1.0, 1.0);
    RooRealVar c3_init("c3_init", "c3", 0.0, -1.0, 1.0);
    RooChebychev pdfBkgInit("pdfBkgInit", "Background", Mrecoil,
                            RooArgList(c0_init, c1_init, c2_init, c3_init));
    pdfBkgInit.fitTo(*dsBkgTemplate, RooFit::Range("fitRange"), RooFit::SumW2Error(kTRUE),
                     RooFit::PrintLevel(-1));

    // 4. Параметры модели
    RooRealVar c0("c0", "c0", c0_init.getVal(), -1.0, 1.0);
    RooRealVar c1("c1", "c1", c1_init.getVal(), -1.0, 1.0);
    RooRealVar c2("c2", "c2", c2_init.getVal(), -1.0, 1.0);
    RooRealVar c3("c3", "c3", c3_init.getVal(), -1.0, 1.0);
    RooChebychev pdfBkg("pdfBkg", "Chebyshev Background", Mrecoil, RooArgList(c0, c1, c2, c3));

    double sumW_S = dsSignalTemplate->sumEntries();
    double sumW_B = dsBkgTemplate->sumEntries();

    RooRealVar nS("nS", "Signal yield", 0.0, 0.0, 50.0 * sumW_S);
    RooRealVar nB("nB", "Background yield", sumW_B, 0.0, 10.0 * sumW_B);
    RooAddPdf model("model", "Signal + Background", RooArgList(pdfSig, pdfBkg), RooArgList(nS, nB));

    // Контейнеры для результатов
    std::vector<double> muValues, zValues;

    // 5. Цикл сканирования
    int nSteps = 1 + static_cast<int>((ANALYTICAL_SCAN_MU_MAX - ANALYTICAL_SCAN_MU_MIN) /
                                      ANALYTICAL_SCAN_MU_STEP);
    for (int i = 0; i < nSteps; ++i) {
        double mu = ANALYTICAL_SCAN_MU_MIN + i * ANALYTICAL_SCAN_MU_STEP;

        // Формируем псевдоданные для текущего mu
        RooDataSet *dsData =
            new RooDataSet("dsData_scan", "Pseudo-data", argSet, RooFit::WeightVar(wVar));
        for (const auto &entry : vBkg) {
            if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN ||
                entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
                continue;
            Mrecoil.setVal(entry.first);
            dsData->add(argSet, entry.second);
        }
        for (const auto &entry : v_qqHX) {
            if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN ||
                entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
                continue;
            Mrecoil.setVal(entry.first);
            dsData->add(argSet, entry.second);
        }
        for (const auto &entry : vSignal) {
            if (entry.first < ANALYTICAL_FIT_MRECOIL_MIN ||
                entry.first > ANALYTICAL_FIT_MRECOIL_MAX)
                continue;
            Mrecoil.setVal(entry.first);
            dsData->add(argSet, -entry.second);     // Вычитаем сигнал из qqHX
            dsData->add(argSet, mu * entry.second); // Добавляем mu * сигнал
        }

        // --- ФИТ H0 (только фон) ---
        nS.setVal(0.0);
        nS.setConstant(kTRUE);
        nB.setVal(sumW_B);
        c0.setVal(c0_init.getVal());
        c1.setVal(c1_init.getVal());
        c2.setVal(c2_init.getVal());
        c3.setVal(c3_init.getVal());

        RooFitResult *fitResB =
            model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                        RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::Save(kTRUE));
        double nll_b = fitResB ? fitResB->minNll() : 0.0;
        delete fitResB;

        // --- ФИТ H1 (сигнал + фон) ---
        nS.setConstant(kFALSE);
        nS.setVal(mu * sumW_S);
        nB.setVal(sumW_B);

        RooFitResult *fitResSB =
            model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                        RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::Save(kTRUE));
        double nll_sb = fitResSB ? fitResSB->minNll() : 0.0;
        delete fitResSB;

        // Вычисление значимости
        double q0 = 2.0 * (nll_b - nll_sb);
        if (q0 < 0.0)
            q0 = 0.0; // Защита от численных артефактов
        double z_lrt = std::sqrt(q0);

        muValues.push_back(mu);
        zValues.push_back(z_lrt);

        if (i % 5 == 0 || i == nSteps - 1) {
            std::cout << "[AnalyticalScan] mu = " << std::fixed << std::setprecision(1) << mu
                      << ", Z_LRT = " << std::setprecision(2) << z_lrt << " σ\n";
        }
        delete dsData;
    }

    // 6. ПОИСК mu для 5 сигм (линейная интерполяция)
    double mu_5sigma = -1.0;
    for (size_t i = 0; i < muValues.size() - 1; ++i) {
        if ((zValues[i] <= ANALYTICAL_SCAN_Z_TARGET &&
             zValues[i + 1] >= ANALYTICAL_SCAN_Z_TARGET) ||
            (zValues[i] >= ANALYTICAL_SCAN_Z_TARGET &&
             zValues[i + 1] <= ANALYTICAL_SCAN_Z_TARGET)) {
            double frac = (ANALYTICAL_SCAN_Z_TARGET - zValues[i]) / (zValues[i + 1] - zValues[i]);
            mu_5sigma = muValues[i] + frac * (muValues[i + 1] - muValues[i]);
            break;
        }
    }

    std::cout << "\n[AnalyticalScan] ====================================================\n";
    if (mu_5sigma > 0) {
        std::cout << "[AnalyticalScan] НАЙДЕНО: mu_5sigma = " << std::fixed << std::setprecision(2)
                  << mu_5sigma << " (для достижения 5 σ)\n";
    } else {
        std::cout << "[AnalyticalScan] Значимость 5 σ НЕ ДОСТИГНУТА в заданном диапазоне mu.\n";
        std::cout << "[AnalyticalScan] Попробуйте увеличить ANALYTICAL_SCAN_MU_MAX.\n";
    }
    std::cout << "[AnalyticalScan] ====================================================\n";

    // 7. Визуализация
    TCanvas *cScan = new TCanvas("cAnalyticalMuScan", "Analytical Fit - Z vs mu", 900, 700);
    cScan->SetLeftMargin(0.13);
    cScan->SetRightMargin(0.05);
    cScan->SetBottomMargin(0.12);
    cScan->SetTopMargin(0.08);

    TGraph *grZ = new TGraph(muValues.size(), muValues.data(), zValues.data());
    grZ->SetTitle(";#mu;Z_{LRT} [#sigma]");
    grZ->SetMarkerStyle(20);
    grZ->SetMarkerSize(1.0);
    grZ->SetLineColor(kBlue);
    grZ->SetLineWidth(2);
    grZ->GetYaxis()->SetRangeUser(
        0.0, std::max(6.0, *std::max_element(zValues.begin(), zValues.end()) * 1.2));
    grZ->Draw("APL");

    // Линия целевой значимости (5 сигм)
    TLine *lineZ5 = new TLine(ANALYTICAL_SCAN_MU_MIN, ANALYTICAL_SCAN_Z_TARGET,
                              ANALYTICAL_SCAN_MU_MAX, ANALYTICAL_SCAN_Z_TARGET);
    lineZ5->SetLineColor(kRed);
    lineZ5->SetLineStyle(kDashed);
    lineZ5->SetLineWidth(2);
    lineZ5->Draw("SAME");

    // Вертикальная линия на mu_5sigma
    TLine *lineMu5 = nullptr;
    if (mu_5sigma > 0 && mu_5sigma >= ANALYTICAL_SCAN_MU_MIN &&
        mu_5sigma <= ANALYTICAL_SCAN_MU_MAX) {
        lineMu5 = new TLine(mu_5sigma, 0.0, mu_5sigma, grZ->GetYaxis()->GetXmax());
        lineMu5->SetLineColor(kGreen + 1);
        lineMu5->SetLineStyle(kDashed);
        lineMu5->SetLineWidth(2);
        lineMu5->Draw("SAME");
    }

    TLegend *leg = new TLegend(0.6, 0.75, 0.98, 0.98);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(grZ, "Z_{LRT} vs #mu", "pl");
    leg->AddEntry(lineZ5, "Target: 5 #sigma", "l");
    if (lineMu5)
        leg->AddEntry(lineMu5, Form("#mu_{5#sigma} = %.2f", mu_5sigma), "l");
    leg->Draw();

    std::string outPdf = (fs::path(outputPath) / "analytical_scan_z_vs_mu.pdf").string();
    cScan->SaveAs(outPdf.c_str());
    std::cout << "[AnalyticalScan] График сохранен: " << outPdf << "\n";

    // Очистка
    delete leg;
    delete lineZ5;
    if (lineMu5)
        delete lineMu5;
    delete grZ;
    delete cScan;
    delete dsSignalTemplate;
    delete dsBkgTemplate;
}

// =============================================================================
// ФУНКЦИЯ runMrecoilScanMu
// =============================================================================
//
// Выполняет сканирование по параметру mu (доля сигнала) для построения
// зависимости p-value от mu и определения чувствительности на уровне 95% CL.
//
// ВХОДНЫЕ ДАННЫЕ:
//   vSignal    - MC-события процесса qqHinvi, пары (M_recoil, вес).
//   vBkg       - MC-события всех фоновых процессов, пары (M_recoil, вес).
//   v_qqHX     - MC-события процесса qqHX, пары (M_recoil, вес).
//   outputPath - путь для сохранения PDF с графиком p-value vs mu.
//
// АЛГОРИТМ:
//   1. Для каждого значения mu в диапазоне [SCAN_MU_MIN, SCAN_MU_MAX] с шагом SCAN_MU_STEP:
//      a. Формируем псевдоданные как смесь фона и сигнала с данным mu
//      b. Выполняем фит только фоном (H0: nS=0), получаем NLL_b
//      c. Выполняем фит сигнал+фон (H1: nS,nB свободны), получаем NLL_sb
//      d. Вычисляем тестовую статистику q0 = 2*(NLL_b - NLL_sb)
//      e. Вычисляем p-value = 1 - CDF_χ²(q0; ndf=1)
//   2. Строим график p-value(mu)
//   3. Находим mu_95, при котором p-value = 0.05 (95% CL)
//
// =============================================================================
void runMrecoilScanMu(const std::vector<std::pair<double, double>> &vSignal,
                      const std::vector<std::pair<double, double>> &vBkg,
                      const std::vector<std::pair<double, double>> &v_qqHX,
                      const std::string &outputPath) {

    std::cout << "\n[Scan] ====================================================\n";
    std::cout << "[Scan] Запуск сканирования по mu\n";
    std::cout << "[Scan] Диапазон mu: [" << SCAN_MU_MIN << ", " << SCAN_MU_MAX << "]\n";
    std::cout << "[Scan] Шаг mu: " << SCAN_MU_STEP << "\n";
    std::cout << "[Scan] Целевой CL: " << (1.0 - SCAN_CL_TARGET) * 100
              << "% (p-value = " << SCAN_CL_TARGET << ")\n";

    // =========================================================================
    // ПОДГООВКА: наблюдаемая переменная и шаблоны
    // =========================================================================
    RooRealVar Mrecoil("Mrecoil", "M_{recoil} [GeV]", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);
    Mrecoil.setRange("fitRange", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);

    RooRealVar wVar("eventWeight", "Event weight", 1e-9, 1e9);
    RooArgSet argSet(Mrecoil, wVar);

    // Шаблон сигнала
    RooDataSet *dsSignal =
        new RooDataSet("dsSignal", "Signal template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vSignal) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsSignal->add(argSet, entry.second);
    }

    // Шаблон фона (с учётом qqHX и вычитанием чистого сигнала)
    RooDataSet *dsBkg =
        new RooDataSet("dsBkg", "Background template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second);
    }
    for (const auto &entry : v_qqHX) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second);
    }
    for (const auto &entry : vSignal) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, -entry.second);
    }

    double sumW_S = dsSignal->sumEntries();
    double sumW_B = dsBkg->sumEntries();

    std::cout << "[Scan] sumW_S (сигнал): " << sumW_S << "\n";
    std::cout << "[Scan] sumW_B (фон): " << sumW_B << "\n";

    if (sumW_S <= 0 || sumW_B <= 0) {
        std::cerr << "[Scan] Ошибка: суммарные веса равны нулю.\n";
        delete dsSignal;
        delete dsBkg;
        return;
    }

    // Построение PDF
    auto mirrorOpt = FIT_KEYSPDF_MIRROR ? RooKeysPdf::MirrorBoth : RooKeysPdf::NoMirror;
    RooKeysPdf pdfSignal("pdfSignal_scan", "Signal PDF", Mrecoil, *dsSignal, mirrorOpt,
                         FIT_KEYSPDF_ADAPTIVITY_SIGNAL);
    RooKeysPdf pdfBkg("pdfBkg_scan", "Background PDF", Mrecoil, *dsBkg, mirrorOpt,
                      FIT_KEYSPDF_ADAPTIVITY_BGD);
    pdfSignal.setNormRange("fitRange");
    pdfBkg.setNormRange("fitRange");

    // Параметры нормировки
    RooRealVar nS("nS_scan", "Signal yield", 0.0, 0.0, 10.0 * sumW_S * SCAN_MU_MAX);
    RooRealVar nB("nB_scan", "Background yield", sumW_B, 0.0, 3.0 * sumW_B);

    // Комбинированная модель
    RooAddPdf model("model_scan", "nS*fs + nB*fb", RooArgList(pdfSignal, pdfBkg),
                    RooArgList(nS, nB));
    model.setNormRange("fitRange");

    // Контейнеры для результатов сканирования
    std::vector<double> muValues;
    std::vector<double> pValues;
    std::vector<double> q0Values;
    std::vector<double> significanceValues;

    // =========================================================================
    // СКАНИРОВАНИЕ ПО MU
    // =========================================================================
    int nSteps = static_cast<int>((SCAN_MU_MAX - SCAN_MU_MIN) / SCAN_MU_STEP) + 1;

    for (int i = 0; i < nSteps; ++i) {
        double mu = SCAN_MU_MIN + i * SCAN_MU_STEP;

        // Формируем псевдоданные для данного mu
        RooDataSet *dsData =
            new RooDataSet("dsData_scan", "Pseudo-data", argSet, RooFit::WeightVar(wVar));

        for (const auto &entry : vBkg) {
            if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
                continue;
            Mrecoil.setVal(entry.first);
            dsData->add(argSet, entry.second);
        }
        if (mu > 0.0) {
            for (const auto &entry : vSignal) {
                if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
                    continue;
                Mrecoil.setVal(entry.first);
                dsData->add(argSet, mu * entry.second);
            }
        }

        double sumW_Data = dsData->sumEntries();

        // Фит H0: только фон (nS = 0)
        nS.setVal(0.0);
        nS.setConstant(kTRUE);
        nB.setVal(sumW_B);

        RooFitResult *fitResB =
            model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                        RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1));

        double nll_b = fitResB ? fitResB->minNll() : 0.0;
        if (fitResB)
            delete fitResB;

        // Фит H1: сигнал + фон
        nS.setConstant(kFALSE);
        nS.setVal(mu * sumW_S);
        nB.setVal(sumW_B);

        RooFitResult *fitResSB =
            model.fitTo(*dsData, RooFit::Extended(kTRUE), RooFit::Range("fitRange"),
                        RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1));

        double nll_sb = fitResSB ? fitResSB->minNll() : 0.0;
        if (fitResSB)
            delete fitResSB;

        // Вычисляем q0 и p-value
        double q0 = 2.0 * (nll_b - nll_sb);
        if (q0 < 0.0)
            q0 = 0.0;

        // p-value из χ² распределения с 1 степенью свободы
        // P(χ² > q0) = 1 - CDF(q0) = 1 - (1 - erf(sqrt(q0/2))) = erf(sqrt(q0/2))
        double sqrt_q0 = std::sqrt(q0);
        double pValue = TMath::Prob(q0, 1);

        double significance = std::sqrt(q0);

        muValues.push_back(mu);
        pValues.push_back(pValue);
        q0Values.push_back(q0);
        significanceValues.push_back(significance);

        if (i % 10 == 0 || i == nSteps - 1) {
            std::cout << "[Scan] mu = " << mu << ", q0 = " << q0 << ", p-value = " << pValue
                      << ", Z = " << significance << " σ\n";
        }

        delete dsData;
    }

    // =========================================================================
    // ПОИСК MU_95 (где p-value = 0.05)
    // =========================================================================
    double mu95 = -1.0;
    for (size_t i = 0; i < muValues.size() - 1; ++i) {
        if ((pValues[i] >= SCAN_CL_TARGET && pValues[i + 1] <= SCAN_CL_TARGET) ||
            (pValues[i] <= SCAN_CL_TARGET && pValues[i + 1] >= SCAN_CL_TARGET)) {
            // Линейная интерполяция
            double frac = (SCAN_CL_TARGET - pValues[i]) / (pValues[i + 1] - pValues[i]);
            mu95 = muValues[i] + frac * (muValues[i + 1] - muValues[i]);
            break;
        }
    }

    std::cout << "\n[Scan] ====================================================\n";
    if (mu95 > 0) {
        std::cout << "[Scan] Найдено mu_95 = " << mu95 << " (95% CL)\n";
    } else {
        std::cout << "[Scan] mu_95 не найдено в заданном диапазоне\n";
    }
    std::cout << "[Scan] ====================================================\n";

    // =========================================================================
    // ВИЗУАЛИЗАЦИЯ
    // =========================================================================
    TCanvas *cScan = new TCanvas("cMuScan", "Muon Scan - p-value vs mu", 900, 700);
    cScan->SetLeftMargin(0.13);
    cScan->SetRightMargin(0.05);
    cScan->SetBottomMargin(0.12);
    cScan->SetTopMargin(0.08);

    // Создаём график
    TGraph *grPValue = new TGraph(muValues.size(), muValues.data(), pValues.data());
    grPValue->SetTitle(";#mu;p-value");
    grPValue->SetMarkerStyle(20);
    grPValue->SetMarkerSize(1.0);
    grPValue->SetLineColor(kBlue);
    grPValue->SetLineWidth(2);

    // Устанавливаем диапазон Y
    grPValue->GetYaxis()->SetRangeUser(0.0, 1.0);
    grPValue->GetXaxis()->SetTitle("#mu");
    grPValue->GetYaxis()->SetTitle("p-value");

    grPValue->Draw("APL");

    // Горизонтальная линия на уровне p-value = 0.05
    TLine *lineCL = new TLine(SCAN_MU_MIN, SCAN_CL_TARGET, SCAN_MU_MAX, SCAN_CL_TARGET);
    lineCL->SetLineColor(kRed);
    lineCL->SetLineStyle(kDashed);
    lineCL->SetLineWidth(2);
    lineCL->Draw("SAME");

    // Вертикальная линия на mu_95 (если найдено)
    TLine *lineMu95 = nullptr;
    if (mu95 > 0 && mu95 >= SCAN_MU_MIN && mu95 <= SCAN_MU_MAX) {
        lineMu95 = new TLine(mu95, 0.0, mu95, 1.0);
        lineMu95->SetLineColor(kGreen + 1);
        lineMu95->SetLineStyle(kDashed);
        lineMu95->SetLineWidth(2);
        lineMu95->Draw("SAME");
    }

    // Текст с результатами
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);

    std::string resultText = mu95 > 0 ? Form("#mu_{95} = %.2f", mu95) : "#mu_{95} not found";
    latex.DrawLatex(0.15, 0.85, resultText.c_str());
    latex.DrawLatex(0.15, 0.80, Form("95%% CL: p-value = %.2f", SCAN_CL_TARGET));

    // Легенда
    TLegend *leg = new TLegend(0.6, 0.8, 0.98, 0.98);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(grPValue, "p-value vs #mu", "pl");
    leg->AddEntry(lineCL, "p-value = 0.05 (95% CL)", "l");
    if (lineMu95)
        leg->AddEntry(lineMu95, Form("#mu_{95} = %.2f", mu95), "l");
    leg->Draw();

    cScan->SaveAs((outputPath + "/pvalue_vs_mu.pdf").c_str());
    std::cout << "[Scan] График сохранен: " << outputPath << "/pvalue_vs_mu.pdf\n";

    // Очистка
    delete leg;
    delete lineCL;
    if (lineMu95)
        delete lineMu95;
    delete grPValue;
    delete cScan;
    delete dsSignal;
    delete dsBkg;
}

// =============================================================================
// ФУНКЦИЯ runMrecoilScanAdapt
// =============================================================================
//
// Выполняет сканирование параметров адаптивности RooKeysPdf и числа бинов
// для шаблонов сигнала и фона.
//
// ВХОДНЫЕ ДАННЫЕ:
//   vSignal    - MC-события процесса qqHinvi, пары (M_recoil, вес).
//   vBkg       - MC-события всех фоновых процессов, пары (M_recoil, вес).
//   v_qqHX     - MC-события процесса qqHX, пары (M_recoil, вес).
//   outputPath - путь для сохранения результатов (CSV и PDF).
//
// АЛГОРИТМ:
//   1. Для каждого значения числа бинов nBins в диапазоне [SCAN_BINS_MIN, SCAN_BINS_MAX]:
//      a. Для каждого значения адаптивности сигнала adaptS в диапазоне [SCAN_ADAPT_SIGNAL_MIN,
//      SCAN_ADAPT_SIGNAL_MAX]:
//         i. Для каждого значения адаптивности фона adaptB в диапазоне [SCAN_ADAPT_BGD_MIN,
//         SCAN_ADAPT_BGD_MAX]:
//            - Строим шаблоны с данными параметрами
//            - Вычисляем chi2/NDF для сигнала и фона
//            - Сохраняем результаты в CSV
//            - Сохраняем графики для комбинаций с chi2 в диапазоне [1, 2]
//
// =============================================================================
void runMrecoilScanAdapt(const std::vector<std::pair<double, double>> &vSignal,
                         const std::vector<std::pair<double, double>> &vBkg,
                         const std::vector<std::pair<double, double>> &v_qqHX,
                         const std::string &outputPath) {

    std::cout << "\n[ScanAdapt] ====================================================\n";
    std::cout << "[ScanAdapt] Запуск сканирования параметров адаптивности и бинов\n";
    std::cout << "[ScanAdapt] Диапазон бинов: [" << SCAN_BINS_MIN << ", " << SCAN_BINS_MAX
              << "] шаг " << SCAN_BINS_STEP << "\n";
    std::cout << "[ScanAdapt] Адаптивность сигнала: [" << SCAN_ADAPT_SIGNAL_MIN << ", "
              << SCAN_ADAPT_SIGNAL_MAX << "] шаг " << SCAN_ADAPT_SIGNAL_STEP << "\n";
    std::cout << "[ScanAdapt] Адаптивность фона: [" << SCAN_ADAPT_BGD_MIN << ", "
              << SCAN_ADAPT_BGD_MAX << "] шаг " << SCAN_ADAPT_BGD_STEP << "\n";

    // Директории для резульатов
    std::string backgroundOutput = (fs::path(outputPath) / "scanFitParams" / "background").string();
    std::string signalOutput = (fs::path(outputPath) / "scanFitParams" / "signal").string();
    fs::create_directories(backgroundOutput);
    fs::create_directories(signalOutput);

    // Подготовка: наблюдаемая переменная
    RooRealVar Mrecoil("Mrecoil", "M_{recoil} [GeV]", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);
    Mrecoil.setRange("fitRange", FIT_MRECOIL_MIN, FIT_MRECOIL_MAX);

    RooRealVar wVar("eventWeight", "Event weight", 1e-9, 1e9);
    RooArgSet argSet(Mrecoil, wVar);

    // Шаблон сигнала
    RooDataSet *dsSignal =
        new RooDataSet("dsSignal", "Signal template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vSignal) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsSignal->add(argSet, entry.second);
    }

    // Шаблон фона (с учётом qqHX и вычитанием чистого сигнала)
    RooDataSet *dsBkg =
        new RooDataSet("dsBkg", "Background template MC", argSet, RooFit::WeightVar(wVar));
    for (const auto &entry : vBkg) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second);
    }
    for (const auto &entry : v_qqHX) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, entry.second);
    }
    for (const auto &entry : vSignal) {
        if (entry.first < FIT_MRECOIL_MIN || entry.first > FIT_MRECOIL_MAX)
            continue;
        Mrecoil.setVal(entry.first);
        dsBkg->add(argSet, -entry.second);
    }

    double sumW_S = dsSignal->sumEntries();
    double sumW_B = dsBkg->sumEntries();

    std::cout << "[ScanAdapt] sumW_S (сигнал): " << sumW_S << "\n";
    std::cout << "[ScanAdapt] sumW_B (фон): " << sumW_B << "\n";

    if (sumW_S <= 0 || sumW_B <= 0) {
        std::cerr << "[ScanAdapt] Ошибка: суммарные веса равны нулю.\n";
        delete dsSignal;
        delete dsBkg;
        return;
    }

    int totalIterations = 0;

    // Скан по параметрам
    for (int nBins = SCAN_BINS_MIN; nBins <= SCAN_BINS_MAX; nBins += SCAN_BINS_STEP) {
        for (double adaptS = SCAN_ADAPT_SIGNAL_MIN; adaptS <= SCAN_ADAPT_SIGNAL_MAX;
             adaptS += SCAN_ADAPT_SIGNAL_STEP) {
            for (double adaptB = SCAN_ADAPT_BGD_MIN; adaptB <= SCAN_ADAPT_BGD_MAX;
                 adaptB += SCAN_ADAPT_BGD_STEP) {

                totalIterations++;

                // Построение PDF с текущими параметрами
                auto mirrorOpt = FIT_KEYSPDF_MIRROR ? RooKeysPdf::MirrorBoth : RooKeysPdf::NoMirror;
                RooKeysPdf pdfSignal("pdfSignal_scan", "Signal PDF", Mrecoil, *dsSignal, mirrorOpt,
                                     adaptS);
                RooKeysPdf pdfBkg("pdfBkg_scan", "Background PDF", Mrecoil, *dsBkg, mirrorOpt,
                                  adaptB);
                pdfSignal.setNormRange("fitRange");
                pdfBkg.setNormRange("fitRange");

                // Проверка качества шаблона сигнала
                RooPlot *frameSig = Mrecoil.frame(RooFit::Title("Signal Template Quality Check"));
                dsSignal->plotOn(frameSig, RooFit::Binning(nBins), RooFit::MarkerStyle(22),
                                 RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
                pdfSignal.plotOn(frameSig, RooFit::LineColor(kRed), RooFit::LineWidth(2));
                double chi2S = frameSig->chiSquare();

                // Проверка качества шаблона фона
                RooPlot *frameBkg =
                    Mrecoil.frame(RooFit::Title("Background Template Quality Check"));
                dsBkg->plotOn(frameBkg, RooFit::Binning(nBins), RooFit::MarkerStyle(21),
                              RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
                pdfBkg.plotOn(frameBkg, RooFit::LineColor(kBlue), RooFit::LineWidth(2));
                double chi2B = frameBkg->chiSquare();

                // Сохраняем графики
                std::string sigPlotName =
                    Form("template_signal_nBins%d_aS%.2f_aB%.2f.pdf", nBins, adaptS, adaptB);
                TCanvas *cSig = new TCanvas("cSig_scan", "Signal Template", 800, 600);
                frameSig->Draw();
                TLatex latexS;
                latexS.SetNDC();
                latexS.SetTextFont(42);
                latexS.SetTextSize(0.04);
                latexS.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.3f", chi2S));
                latexS.DrawLatex(0.15, 0.80, Form("n_{bins} = %d", nBins));
                latexS.DrawLatex(0.15, 0.75, Form("adapt = %.1f", adaptS));
                cSig->SaveAs((signalOutput + "/" + sigPlotName).c_str());

                std::string bkgPlotName =
                    Form("template_background_nBins%d_aS%.2f_aB%.2f.pdf", nBins, adaptS, adaptB);
                TCanvas *cBkg = new TCanvas("cBkg_scan", "Background Template", 800, 600);
                frameBkg->Draw();
                TLatex latexB;
                latexB.SetNDC();
                latexB.SetTextFont(42);
                latexB.SetTextSize(0.04);
                latexB.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.3f", chi2B));
                latexS.DrawLatex(0.15, 0.80, Form("n_{bins} = %d", nBins));
                latexS.DrawLatex(0.15, 0.75, Form("adapt = %.1f", adaptB));
                cBkg->SaveAs((backgroundOutput + "/" + bkgPlotName).c_str());

                delete cSig;
                delete cBkg;
                delete frameSig;
                delete frameBkg;

                if (totalIterations % 20 == 0) {
                    std::cout << "[ScanAdapt] Прогресс: " << totalIterations << " итераций."
                              << "\n";
                }
            }
        }
    }

    delete dsSignal;
    delete dsBkg;

    std::cout << "\n[ScanAdapt] ====================================================\n";
    std::cout << "[ScanAdapt] Сканирование завершено\n";
    std::cout << "[ScanAdapt] Всего итераций: " << totalIterations << "\n";
    std::cout << "[ScanAdapt] Графики сохранены в: " << outputPath << "/template_*.pdf\n";
    std::cout << "[ScanAdapt] ====================================================\n";
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
    std::cout
        << "Многофайловый анализ ee -> ZH -> qq + invisible\n\n"
        << "Использование:\n"
        << "  " << progName << " file1.root file2.root ... [options]\n\n"
        << "Опции:\n"
        << "  -h, --help                   Показать эту справку\n"
        << "  -o, --output-dir DIR         Базовая директория результатов (по умолчанию: "
           "../pdf_results)\n"
        << "  --export-csv DIR             Экспорт данных в CSV после предотборов (без "
           "гистограмм)\n"
        << "  -use-bdt                     Использовать BDT вместо основных отборов\n"
        << "  -s, --scan-mu                Запустить сканирование по mu\n"
        << "  --template-fit               Выполнить шаблонный фит\n"
        << "  --analytical-fit             Выполнить аналитический фит (Gauss + Polynomial)\n"
        << "  --analytical-scan-mu         Запустить сканирование mu для аналитического фита "
           "(поиск 5 #sigma)\n"
        << "  -a, --scan-fit-params        Запустить сканирование параметров адаптивности и бинов\n"
        << "\nПример:\n"
        << "  " << progName
        << " merged_E240_qqHX.root merged_E240_qq.root merged_E240_qqHinvi.root\n";
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

    std::vector<std::string> inputFiles;
    std::string outputBaseDir = OUTPUT_BASE_DIR;
    std::string exportCsvDir = "../ml";
    std::string bdtModelPath = DEFAULT_BDT_MODEL_PATH;
    double bdtThreshold = DEFAULT_BDT_THRESHOLD;
    bool runScanMu = false;
    bool runAnalyticalScanMu = false;
    bool runScanFitParams = false;
    bool runFit = false;
    bool runAnalyticalFit = false;
    bool exportCsv = false;
    bool useBdt = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 < argc)
                outputBaseDir = argv[++i];
        } else if (arg == "--export-csv") {
            exportCsv = true;
        } else if (arg == "--use-bdt") {
            useBdt = true;
        } else if (arg == "--template-fit") {
            runFit = true;
        } else if (arg == "--analytical-fit") {
            runAnalyticalFit = true;
        } else if (arg == "--analytical-scan-mu") {
            runAnalyticalScanMu = true;
        } else if (arg == "-s" || arg == "--scan-mu") {
            runScanMu = true;
        } else if (arg == "-a" || arg == "--scan-fit-params") {
            runScanFitParams = true;
        } else if (arg[0] != '-') {
            inputFiles.push_back(arg);
        } else {
            std::cerr << "Неизвестный аргумент: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    if (inputFiles.empty()) {
        std::cerr << "Ошибка: не указано ни одного ROOT-файла\n";
        printUsage(argv[0]);
        return 1;
    }

    // Проверка взаимоисключения режимов
    if (exportCsv && useBdt) {
        std::cerr << "Ошибка: Режимы --export-csv и --use-bdt взаимоисключающие.\n";
        return 1;
    }

    // Файл с дампом данных для обучения
    std::ofstream csvFile;
    bool csvInitialized = false;
    if (exportCsv) {
        fs::create_directories(exportCsvDir);
        std::string csvPath = (fs::path(exportCsvDir) / "ml_data.csv").string();
        csvFile.open(csvPath, std::ios::out | std::ios::trunc); // trunc очищает файл при старте
        if (!csvFile.is_open()) {
            std::cerr << "Ошибка: не удалось создать CSV файл " << csvPath << std::endl;
            return 1;
        }
        // Заголовок CSV
        csvFile << "process,is_signal,weight,"
                << "invMass,cosThetaZ,deltaR,"
                << "cosTheta1,cosTheta2,jet1_E,jet2_E,"
                << "met_jet,dijetEnergy,deltaTheta,deltaPhi,"
                << "met_pfo,pmiss_mag,cosThetaPmiss\n";
        csvInitialized = true;
        std::cout << "[CSV] Файл для выгрузки открыт: " << csvPath << std::endl;
    }

    // Загрузка модели
    BoosterHandle h_booster = nullptr;
    if (useBdt) {
        int ret = XGBoosterCreate(nullptr, 0, &h_booster);
        if (ret != 0) {
            std::cerr << "Ошибка: Не удалось создать Booster: " << XGBGetLastError() << "\n";
            return 1;
        }

        ret = XGBoosterLoadModel(h_booster, bdtModelPath.c_str());
        if (ret != 0) {
            std::cerr << "Ошибка: Не удалось загрузить модель " << bdtModelPath << ": "
                      << XGBGetLastError() << "\n";
            XGBoosterFree(h_booster);
            return 1;
        }
        std::cout << "XGBoost модель загружена: " << bdtModelPath << "\n";
        std::cout << "  Порог отбора: " << bdtThreshold << "\n";
    }

    // Контейнер для отдельных гистограмм каждого процесса + их метаданные
    std::map<std::string, std::pair<TH1F *, ProcessInfo>> processRecoilHists;

    // Контейнеры для накопления Mrecoil после всех отборов
    std::vector<std::pair<double, double>> vMrecoil_Signal_Weighted;
    std::vector<std::pair<double, double>> vMrecoil_Bkg_Weighted;
    // Для qqHX (сигнал + фон) нужен отдельный контейнер, см. функцию фита
    std::vector<std::pair<double, double>> vMrecoil_qqHX_Weighted;

    // Сравнительные гистограммы предотборов (общие для всех процессов)
    std::map<std::string, std::map<std::string, TH1F *>> preselectionHists;

    // Суммарные гистограммы для основных отборов (общие для всех процессов)
    std::map<std::string, std::map<std::string, TH1F *>> mainSelHists;

    // Цикл по входным файлам
    auto processDB = getProcessDatabase();
    for (const auto &inputRootFile : inputFiles) {

        // Ищем такой процесс в базе процессов
        auto it = processDB.find(fs::path(inputRootFile).filename().string());
        ProcessInfo proc;
        if (it != processDB.end()) {
            proc = it->second;
        } else {
            std::cout << "Предупреждение: файл " << inputRootFile
                      << " не найден в базе. Используется weight = 1.0" << std::endl;
            proc.legendName = extractProcessName(inputRootFile);
            proc.weight = 1.0;
            proc.color = kGray;
        }

        // Извлечение имени процесса и подготовка путей
        std::string processName = extractProcessName(inputRootFile);
        if (inputRootFile.empty()) {
            std::cerr << "Ошибка: не указан входной файл" << std::endl;
            printUsage(argv[0]);
            return 1;
        }
        std::cout << "Процесс: " << processName << std::endl;
        std::cout << "Входной файл: " << inputRootFile << std::endl;

        fs::path processOutputDir = fs::path(OUTPUT_BASE_DIR) / processName;
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

        if (APPLY_PRE_LEPTON_VETO || APPLY_PRE_HIGH_E_PHOTON_VETO ||
            APPLY_PRE_ISOLATED_PHOTON_VETO) {
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
            new TH1F("hRecoilMass", "Recoil Mass Against Two Jets;M_{recoil} [GeV];Events",
                     RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);

        TH2F *h2D_Correlation = new TH2F(
            "h2D_Correlation", "Invariant Mass vs Recoil Mass;M_{jj} [GeV];M_{recoil} [GeV]",
            MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);

        TH1F *hCosThetaZ =
            new TH1F("hCosThetaZ ", "cos#theta of Z Boson (Two-Jet System);cos#theta_{Z};Events ",
                     COS_THETA_Z_BINS, COS_THETA_Z_MIN, COS_THETA_Z_MAX);

        TH1F *hDeltaR =
            new TH1F("hDeltaR",
                     "Distance #Delta R Between Two Jets;#Delta R = #sqrt{#Delta#eta^{2} + "
                     "#Delta#phi^{2}};Events",
                     DELTA_R_BINS, DELTA_R_MIN, DELTA_R_MAX);

        TH1F *hCosThetaJet = new TH1F("hCosThetaJet", "cos#theta of Jets;cos#theta;Events",
                                      COS_THETA_JET_BINS, COS_THETA_JET_MIN, COS_THETA_JET_MAX);

        TH1F *hDeltaTheta =
            new TH1F("hDeltaTheta", "#Delta#theta between two jets;#Delta#theta [rad];Events",
                     DELTA_THETA_BINS, DELTA_THETA_MIN, DELTA_THETA_MAX);

        TH1F *hDeltaPhi =
            new TH1F("hDeltaPhi", "#Delta#phi between two jets;#Delta#phi [rad];Events",
                     DELTA_PHI_BINS, DELTA_PHI_MIN, DELTA_PHI_MAX);

        TH1F *hMETpfo = new TH1F("hMETpfo", "MET from all PFOs;MET_{PFO} [GeV];Events",
                                 MET_PFO_BINS, MET_PFO_MIN, MET_PFO_MAX);
        TH1F *hMETjet = new TH1F("hMETjet", "MET from Two Jets;MET_{jet} [GeV];Events",
                                 MET_JET_BINS, MET_JET_MIN, MET_JET_MAX);
        TH1F *hPmissMag =
            new TH1F("hPmissMag", "Magnitude of Missing 3-Momentum;|P_{miss}| [GeV];Events",
                     PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV);

        TH1F *hCosThetaPmiss =
            new TH1F("hCosThetaPmiss", "cos#theta of Missing 3-Momentum;cos#theta_{miss};Events",
                     COS_THETA_PMISS_BINS, COS_THETA_PMISS_MIN, COS_THETA_PMISS_MAX);

        TH2F *h2D_Mrecoil_vs_MET = new TH2F(
            "h2D_Mrecoil_vs_MET", "M_{recoil} vs MET_{jet};MET_{jet} [GeV];M_{recoil} [GeV]",
            MET_JET_BINS, MET_JET_MIN, MET_JET_MAX, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);
        TH2F *h2D_Mrecoil_vs_Pmiss = new TH2F(
            "h2D_Mrecoil_vs_Pmiss", "M_{recoil} vs |P_{miss}|;|P_{miss}| [GeV];M_{recoil} [GeV]",
            PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV, RECOIL_BINS, RECOIL_MIN_GEV, RECOIL_MAX_GEV);
        TH2F *h2D_MET_vs_Pmiss = new TH2F(
            "h2D_MET_vs_Pmiss", "MET_{jet} vs |P_{miss}|;|P_{miss}| [GeV];MET_{jet} [GeV]",
            PMISS_BINS, PMISS_MIN_GEV, PMISS_MAX_GEV, MET_JET_BINS, MET_JET_MIN, MET_JET_MAX);
        TH2F *h2D_Mjj_vs_MET =
            new TH2F("h2D_Mjj_vs_MET", "M_{jj} vs MET_{jet};MET_{jet} [GeV];M_{jj} [GeV]",
                     MET_JET_BINS, MET_JET_MIN, MET_JET_MAX, MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV);
        TH2F *h2D_Mjj_vs_Pmiss = new TH2F(
            "h2D_Mjj_vs_Pmiss", "M_{jj} vs |P_{miss}|;|P_{miss}| [GeV];M_{jj} [GeV]", PMISS_BINS,
            PMISS_MIN_GEV, PMISS_MAX_GEV, MASS_BINS, MASS_MIN_GEV, MASS_MAX_GEV);
        TH2F *h2D_CosThetaZ_vs_CosThetaPmiss =
            new TH2F("h2D_CosThetaZ_vs_CosThetaPmiss",
                     "cos#theta_{Z} vs cos#theta_{miss};cos#theta_{miss};cos#theta_{Z}",
                     COS_THETA_PMISS_BINS, COS_THETA_PMISS_MIN, COS_THETA_PMISS_MAX,
                     COS_THETA_Z_BINS, COS_THETA_Z_MIN, COS_THETA_Z_MAX);
        TH1F *hDijetEnergy =
            new TH1F("hDijetEnergy", "Dijet System Energy;E_{jj} [GeV];Events", DIJET_ENERGY_BINS,
                     DIJET_ENERGY_MIN_GEV, DIJET_ENERGY_MAX_GEV);

        // Создаем гистограммы для предотборов
        TH1F *hPhotonEnergy = new TH1F(Form("hPhotonEnergy_%s", processName.c_str()),
                                       "Photon Energy;E_{#gamma} [GeV];Events", PHOTON_ENERGY_BINS,
                                       PHOTON_ENERGY_MIN, PHOTON_ENERGY_MAX);
        hPhotonEnergy->SetDirectory(0);

        TH1F *hNJets = new TH1F(Form("hNJets_%s", processName.c_str()),
                                "Number of Jets;N_{jets};Events", NJETS_BINS, NJETS_MIN, NJETS_MAX);
        hNJets->SetDirectory(0);

        TH1F *hNConstituents = new TH1F(Form("hNConstituents_%s", processName.c_str()),
                                        "Number of Constituents per Jet;N_{const};Events",
                                        NCONST_BINS, NCONST_MIN, NCONST_MAX);
        hNConstituents->SetDirectory(0);

        // Сохраняем гистограммы предотборов в контейнер
        preselectionHists[processName]["photonEnergy"] = hPhotonEnergy;
        preselectionHists[processName]["nJets"] = hNJets;
        preselectionHists[processName]["nConstituents"] = hNConstituents;

        // Гистограммы основных отборов
        TH1F *hMainMET = new TH1F(Form("hMainMET_%s", processName.c_str()),
                                  "MET from Two Jets;MET_{jet} [GeV];Events", MET_JET_BINS,
                                  MET_JET_MIN, MET_JET_MAX);
        hMainMET->SetDirectory(0);

        TH1F *hMainDeltaPhi = new TH1F(Form("hMainDeltaPhi_%s", processName.c_str()),
                                       "#Delta#phi between jets;#Delta#phi [rad];Events",
                                       DELTA_PHI_BINS, DELTA_PHI_MIN, DELTA_PHI_MAX);
        hMainDeltaPhi->SetDirectory(0);

        TH1F *hMainCosThetaZ = new TH1F(Form("hMainCosThetaZ_%s", processName.c_str()),
                                        "cos#theta_{Z};cos#theta_{Z};Events", COS_THETA_Z_BINS,
                                        COS_THETA_Z_MIN, COS_THETA_Z_MAX);
        hMainCosThetaZ->SetDirectory(0);

        TH1F *hMainDijetMass = new TH1F(Form("hMainDijetMass_%s", processName.c_str()),
                                        "Invariant Mass of Two Jets;M_{jj} [GeV];Events", MASS_BINS,
                                        MASS_MIN_GEV, MASS_MAX_GEV);
        hMainDijetMass->SetDirectory(0);

        TH1F *hMainPmiss = new TH1F(Form("hMainPmiss_%s", processName.c_str()),
                                    "Missing Momentum;|P_{miss}| [GeV];Events", PMISS_BINS,
                                    PMISS_MIN_GEV, PMISS_MAX_GEV);
        hMainPmiss->SetDirectory(0);

        mainSelHists[processName]["met"] = hMainMET;
        mainSelHists[processName]["deltaPhi"] = hMainDeltaPhi;
        mainSelHists[processName]["cosThetaZ"] = hMainCosThetaZ;
        mainSelHists[processName]["dijetMass"] = hMainDijetMass;
        mainSelHists[processName]["pmiss"] = hMainPmiss;

        // Взвешенная гистограмма массы отдачи для текущего процесса
        TH1F *hRecoilMassWeight = new TH1F(
            ("hRecoil_" + processName).c_str(), "Weight Recoil Mass;M_{recoil} [GeV];Events",
            RECOIL_STACK_BINS, RECOIL_STACK_MIN_GEV, RECOIL_STACK_MAX_GEV);
        hRecoilMassWeight->SetDirectory(0); // Отключаем автоматическое удаление при закрытии файла
        processRecoilHists[processName] = {hRecoilMassWeight, proc};

        // Статистики
        CutStatistics stats;
        IsoElectronStats elecStats;

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

            // ======== ЗАПОЛНЕНИЕ СРАВНИТЕЛЬНЫХ ГИСТОГРАММ ПРЕДОТБОРОВ ========
            // Распределение энергии фотонов (максимальная энергия фотона в событии)
            if (particleType && pfoE) {
                double maxPhotonE = 0.0;
                for (size_t k = 0; k < particleType->size(); ++k) {
                    if (std::abs(particleType->at(k)) == PDG_PHOTON) {
                        if (pfoE->at(k) > maxPhotonE) {
                            maxPhotonE = pfoE->at(k);
                        }
                    }
                }
                if (maxPhotonE > 0) {
                    preselectionHists[processName]["photonEnergy"]->Fill(maxPhotonE);
                }
            }

            // Распределение числа джетов
            if (inclJetE) {
                int nJets = static_cast<int>(inclJetE->size());
                preselectionHists[processName]["nJets"]->Fill(nJets);
            }

            // Распределение числа конституентов в джетах
            if (inclJetSize) {
                for (size_t j = 0; j < inclJetSize->size(); ++j) {
                    int nConst = static_cast<int>(inclJetSize->at(j));
                    preselectionHists[processName]["nConstituents"]->Fill(nConst);
                }
            }

            // ==================== ПРЕДОТБОРЫ ====================
            if (APPLY_PRE_LEPTON_VETO &&
                hasIsolatedLeptonROOT_FSR(particleType, pfoE, pfoPx, pfoPy, pfoPz))
                continue;
            stats.afterPreLeptonVeto++;

            if (APPLY_PRE_HIGH_E_PHOTON_VETO &&
                hasHighEnergyPhoton(particleType, pfoE, PHOTON_ENERGY_CUT_GEV))
                continue;
            stats.afterPreHighEPhotonVeto++;

            if (APPLY_PRE_ISOLATED_PHOTON_VETO &&
                hasIsolatedPhoton(particleType, pfoE, pfoPx, pfoPy, pfoPz,
                                  PHOTON_ISO_MIN_ENERGY_GEV, PHOTON_ISO_COS_CONE_ANGLE,
                                  PHOTON_ISO_MAX_CONE_ENERGY_GEV))
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

            // ==================== КИНЕМАТИКА И ЗАПОЛНЕНИЕ ГИСТОГРАММ ====================
            // Событие прошло все предотборы, строим распределения
            TLorentzVector jet1(inclJetPx->at(0), inclJetPy->at(0), inclJetPz->at(0),
                                inclJetE->at(0));
            TLorentzVector jet2(inclJetPx->at(1), inclJetPy->at(1), inclJetPz->at(1),
                                inclJetE->at(1));
            TLorentzVector dijet = jet1 + jet2;

            double invMass = dijet.M();
            double recoilMass = calculateRecoilMass(dijet, SQRT_S_GEV);
            double cosThetaZ = std::cos(calculatePolarAngle(dijet));
            double deltaR = calculateDeltaR(jet1, jet2);
            double cosTheta1 = (jet1.P() > 1e-9) ? jet1.Pz() / jet1.P() : 0.0;
            double cosTheta2 = (jet2.P() > 1e-9) ? jet2.Pz() / jet2.P() : 0.0;
            double met_jet = dijet.Pt();
            double dijetEnergy = dijet.E();
            double theta1 = calculatePolarAngle(jet1);
            double theta2 = calculatePolarAngle(jet2);
            double deltaTheta = std::abs(theta1 - theta2);
            double phi1 = jet1.Phi();
            double phi2 = jet2.Phi();
            double deltaPhi = std::abs(phi1 - phi2);
            if (deltaPhi > M_PI)
                deltaPhi = 2 * M_PI - deltaPhi;

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

            // Детерминированное разбиение на Analysis и ML samples
            // Используем хеш от номера события, чтобы избежать корреляций внутри батчей MC
            // Результат всегда одинаков для одного и того же события
            uint64_t eventHash = std::hash<Long64_t>{}(tree->GetReadEntry());
            bool isAnalysisSample = (eventHash % 2 == 0); // 50/50 split

            // Экспорт в CSV если включен режим экспорта
            if (exportCsv && csvInitialized) {
                bool isSignal = (processName.find("qqHinvi") != std::string::npos);

                // Экспортируем половину всех событий для обучения
                if (!isAnalysisSample) {
                    csvFile << processName << "," << (isSignal ? 1 : 0) << "," << proc.weight << ","
                            << invMass << "," << cosThetaZ << "," << deltaR << "," << cosTheta1
                            << "," << cosTheta2 << "," << inclJetE->at(0) << "," << inclJetE->at(1)
                            << "," << met_jet << "," << dijetEnergy << "," << deltaTheta << ","
                            << deltaPhi << "," << met_pfo << "," << pmiss_mag << ","
                            << cosThetaPmiss << "\n";
                }

                // Пропускаем основные отборы и гистограммы
                continue;
            }

            // В режиме анализа обрабатываем только те события, которые не участвовали в
            // обучении. Это только в случае применения bdt, в обычном режиме обрабатываем
            // все события.
            if (useBdt && !isAnalysisSample)
                continue;

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
            hDijetEnergy->Fill(dijetEnergy);
            hDeltaTheta->Fill(deltaTheta);
            hDeltaPhi->Fill(deltaPhi);

            h2D_Mrecoil_vs_MET->Fill(met_jet, recoilMass);
            h2D_Mrecoil_vs_Pmiss->Fill(pmiss_mag, recoilMass);
            h2D_MET_vs_Pmiss->Fill(pmiss_mag, met_jet);
            h2D_Mjj_vs_MET->Fill(met_jet, invMass);
            h2D_Mjj_vs_Pmiss->Fill(pmiss_mag, invMass);
            h2D_CosThetaZ_vs_CosThetaPmiss->Fill(cosThetaPmiss, cosThetaZ);

            // Заполнение суммарных гистограмм для основных отборов (после предотборов)
            mainSelHists[processName]["met"]->Fill(met_jet);
            mainSelHists[processName]["deltaPhi"]->Fill(deltaPhi);
            mainSelHists[processName]["cosThetaZ"]->Fill(cosThetaZ);
            mainSelHists[processName]["dijetMass"]->Fill(invMass);
            mainSelHists[processName]["pmiss"]->Fill(pmiss_mag);

            // ==================== ОСНОВНЫЕ ОТБОРЫ ИЛИ BDT ====================
            if (useBdt) {
                // 1. Формируем вектор фичей
                std::vector<float> features(NUM_BDT_FEATURES);
                features[0] = static_cast<float>(invMass);
                features[1] = static_cast<float>(cosThetaZ);
                features[2] = static_cast<float>(deltaR);
                features[3] = static_cast<float>(cosTheta1);
                features[4] = static_cast<float>(cosTheta2);
                features[5] = static_cast<float>(inclJetE->at(0));
                features[6] = static_cast<float>(inclJetE->at(1));
                features[7] = static_cast<float>(met_jet);
                features[8] = static_cast<float>(dijetEnergy);
                features[9] = static_cast<float>(deltaTheta);
                features[10] = static_cast<float>(deltaPhi);
                features[11] = static_cast<float>(met_pfo);
                features[12] = static_cast<float>(pmiss_mag);
                features[13] = static_cast<float>(cosThetaPmiss);

                // 2. Создаем временный DMatrix из массива float
                // Аргументы: data, nrow, ncol, missing_value, &handle
                DMatrixHandle dmat = nullptr;
                int ret = XGDMatrixCreateFromMat(features.data(), 1, NUM_BDT_FEATURES, 0.0f, &dmat);

                if (ret != 0) {
                    std::cerr << "Ошибка создания DMatrix: " << XGBGetLastError() << "\n";
                    continue;
                }

                // 3. Делаем предсказание
                bst_ulong out_len = 0;
                const float *result = nullptr;
                ret = XGBoosterPredict(h_booster, dmat, 0, 0, 0, &out_len, &result);

                // 4. Освобождаем память DMatrix сразу после использования
                XGDMatrixFree(dmat);

                if (ret != 0) {
                    std::cerr << "Ошибка предсказания XGBoost: " << XGBGetLastError() << "\n";
                    continue;
                }

                // 5. Получаем скор
                // result указывает на внутренний буфер XGBoost, копировать не нужно, пока не
                // освободим booster
                float bdtScore = result[0];

                // 6. Применяем порог
                if (bdtScore < bdtThreshold)
                    continue;
                stats.afterBdt++;
            } else {
                if (APPLY_MAIN_MET_CUT && (met_jet < MET_CUT_MIN_GEV || met_jet > MET_CUT_MAX_GEV))
                    continue;
                stats.afterMetCut++;

                if (APPLY_MAIN_DELTA_PHI_CUT && deltaPhi >= DELTA_PHI_CUT_MAX)
                    continue;
                stats.afterDeltaPhiCut++;

                if (APPLY_MAIN_COS_THETA_Z_CUT && std::abs(cosThetaZ) >= COS_THETA_Z_CUT)
                    continue;
                stats.afterCosThetaZCut++;

                if (APPLY_MAIN_DIJET_MASS_WINDOW &&
                    (invMass < DIJET_MASS_WINDOW_MIN_GEV || invMass > DIJET_MASS_WINDOW_MAX_GEV))
                    continue;
                stats.afterDijetMassWindow++;

                if (APPLY_MAIN_PMISS_CUT &&
                    (pmiss_mag < PMISS_CUT_MIN_GEV || pmiss_mag > PMISS_CUT_MAX_GEV))
                    continue;
                stats.afterPmissCut++;

                if (APPLY_MAIN_RECOIL_MASS_WINDOW && (recoilMass < RECOIL_MASS_WINDOW_MIN_GEV ||
                                                      recoilMass > RECOIL_MASS_WINDOW_MAX_GEV))
                    continue;
                stats.afterRecoilMassWindow++;

                if (APPLY_MAIN_ELLIPSE_CUT &&
                    !isInsideEllipse(invMass, recoilMass, ELLIPSE_CX_GEV, ELLIPSE_CY_GEV,
                                     ELLIPSE_A_GEV, ELLIPSE_B_GEV, ELLIPSE_THETA))
                    continue;
                stats.afterEllipseCut++;
            }

            stats.finalSelected++;

            // Если событие прошло все отборы, то заполняем гистограмму текущего процесса. Веса
            // применим уже при построении стековой гистограммы
            hRecoilMassWeight->Fill(recoilMass);

            // Вспомогательная лямбда функция для отбора только тех процессов, у которых есть
            // достаточная статистика для дальнейшего анализа.
            auto isProcessAllowed = [](const std::string &name) {
                return std::find(RECOIL_STACK_ORDER.begin(), RECOIL_STACK_ORDER.end(), name) !=
                       RECOIL_STACK_ORDER.end();
            };

            // Определяем, сигнал это или фон (по имени процесса). Заполняем соотвествующий
            // контейнер с весом. Добавляем в данные для фита только процессы с достаточной
            // статистикой после отборов. Процесс qqHX сохраняем в отдельный контейнер потому
            // что он содержит и фон и сигнал.
            bool isSignal = (processName.find("qqHinvi") != std::string::npos);
            bool is_qqHX = (processName.find("qqHX") != std::string::npos);

            if (isSignal)
                vMrecoil_Signal_Weighted.emplace_back(recoilMass, proc.weight);
            else if (is_qqHX)
                vMrecoil_qqHX_Weighted.emplace_back(recoilMass, proc.weight);
            else if (isProcessAllowed(processName))
                vMrecoil_Bkg_Weighted.emplace_back(recoilMass, proc.weight);
        }

        // Итоговая статистика
        auto endTime = std::chrono::high_resolution_clock::now();
        auto totalSec =
            std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

        std::cout << "\n=======================================================" << std::endl;
        std::cout << "Обработка завершена!" << std::endl;
        std::cout << "Всего событий: " << nEntries << std::endl;
        std::cout << "Прошло времени: " << totalSec << " с (" << totalSec / 60.0 << " мин)"
                  << std::endl;

        stats.print(processName, useBdt);
        elecStats.print();

        // Пропускаем отрисовку гистограмм, если включен режим дампа или bdt
        if (exportCsv) {
            std::cout << "Режим экспорта CSV: гистограммы не строятся" << std::endl;
        } else if (useBdt) {
            std::cout << "Режим BDT: гистограммы не строятся" << std::endl;
        } else {
            // Отрисовка гистограмм с условным отображением основных отборов
            std::vector<std::pair<double, std::string>> invMassMarks = {{MZ_GEV, "M_{Z}"}};
#if APPLY_MAIN_DIJET_MASS_WINDOW
            invMassMarks.emplace_back(DIJET_MASS_WINDOW_MIN_GEV, "M_{jj}^{min}");
            invMassMarks.emplace_back(DIJET_MASS_WINDOW_MAX_GEV, "M_{jj}^{max}");
#endif
            drawHistogram1D(hInvMass, "cInvMass", "M_{jj} [GeV]", OUTPUT_INV_MASS, invMassMarks,
                            kRed, 2);

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
                            ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV, ELLIPSE_B_GEV,
                            ELLIPSE_THETA, true
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

            drawHistogram1D(hCosThetaJet, "cCosThetaJet", "cos#theta", OUTPUT_COS_THETA_JET, {},
                            kCyan, 2);

            drawHistogram1D(hMETpfo, "cMETpfo", "MET_{PFO} [GeV]", OUTPUT_MET_PFO, {}, kOrange + 1,
                            2);

            std::vector<std::pair<double, std::string>> metMarks;
#if APPLY_MAIN_MET_CUT
            metMarks.emplace_back(MET_CUT_MIN_GEV, "MET_{min}");
#endif
            drawHistogram1D(hMETjet, "cMETjet", "MET_{jet} [GeV]", OUTPUT_MET_JET, metMarks,
                            kViolet, 2);

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
            drawHistogram2D(h2D_MET_vs_Pmiss, "c2D_MET_vs_Pmiss", "|P_{miss}| [GeV]",
                            "MET_{jet} [GeV]", makeOutputPath("2D_MET_vs_Pmiss"),
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
            drawHistogram2D(h2D_Mjj_vs_Pmiss, "c2D_Mjj_vs_Pmiss", "|P_{miss}| [GeV]",
                            "M_{jj} [GeV]", makeOutputPath("2D_Mjj_vs_Pmiss"));

            // h2D_CosThetaZ_vs_CosThetaPmiss: горизонтальные линии |cosθ_Z| < 0.98
            drawHistogram2D(h2D_CosThetaZ_vs_CosThetaPmiss, "c2D_CosThetaZ_vs_CosThetaPmiss",
                            "cos#theta_{miss}", "cos#theta_{Z}",
                            makeOutputPath("2D_CosThetaZ_vs_CosThetaPmiss"),
#if APPLY_MAIN_COS_THETA_Z_CUT
                            -1, COS_THETA_Z_CUT, "", "|cos#theta|^{cut}");
#else
                            -1, -1, "", "");
#endif

            std::vector<std::pair<double, std::string>> deltaPhiMarks;
#if APPLY_MAIN_DELTA_PHI_CUT
            deltaPhiMarks.emplace_back(DELTA_PHI_CUT_MAX, "#Delta#phi^{cut}");
#endif
            drawHistogram1D(hDeltaPhi, "cDeltaPhi", "#Delta#phi [rad]",
                            makeOutputPath("deltaPhi_jets"), deltaPhiMarks, kGreen + 2, 2);

            // Остальные гистограммы без линий отборов
            drawHistogram1D(hPmissMag, "cPmissMag", "|P_{miss}| [GeV]", OUTPUT_PMISS_MAG, {},
                            kOrange, 2);
            drawHistogram1D(hCosThetaPmiss, "cCosThetaPmiss", "cos#theta_{miss}",
                            OUTPUT_COS_THETA_PMISS, {}, kMagenta, 2);
            drawHistogram1D(hDijetEnergy, "cDijetEnergy", "E_{jj} [GeV]",
                            makeOutputPath("dijet_energy"), {}, kAzure + 1, 2);
            drawHistogram1D(hDeltaTheta, "cDeltaTheta", "#Delta#theta [rad]",
                            makeOutputPath("deltaTheta_jets"), {}, kOrange + 1, 2);
        }

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
        delete hDijetEnergy;
        delete hDeltaTheta;
        delete hDeltaPhi;
        inputFile->Close();
        delete inputFile;

        if (!useBdt && !exportCsv) {
            std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
                      << std::endl;
        }
    }

    if (!exportCsv && !useBdt) {
        // =============================================================================
        // ОТРИСОВКА ГИСТОГРАММ ПРЕДОТБОРОВ
        // =============================================================================
        std::cout << "\nОтрисовка гистограмм предотборов...\n";

        // Создаем временные карты для отрисовки
        std::map<std::string, std::pair<TH1F *, ProcessInfo>> histsPhotonEnergy;
        std::map<std::string, std::pair<TH1F *, ProcessInfo>> histsNJets;
        std::map<std::string, std::pair<TH1F *, ProcessInfo>> histsNConstituents;

        for (const auto &procEntry : preselectionHists) {
            const std::string &procName = procEntry.first;

            // Находим ProcessInfo из processRecoilHists
            ProcessInfo info;
            auto it = processRecoilHists.find(procName);
            if (it != processRecoilHists.end()) {
                info = it->second.second;
            } else {
                info = ProcessInfo{procName, 1.0, kBlack, 1001};
            }

            if (procEntry.second.count("photonEnergy")) {
                histsPhotonEnergy[procName] = {procEntry.second.at("photonEnergy"), info};
            }
            if (procEntry.second.count("nJets")) {
                histsNJets[procName] = {procEntry.second.at("nJets"), info};
            }
            if (procEntry.second.count("nConstituents")) {
                histsNConstituents[procName] = {procEntry.second.at("nConstituents"), info};
            }
        }

        // Отрисовка
        std::string preselectionDir = (fs::path(outputBaseDir) / "preselection").string();
        fs::create_directories(preselectionDir);

        // 1. Энергия фотонов с отметкой PHOTON_ENERGY_CUT_GEV
        drawPreselectionHistograms(histsPhotonEnergy, "photonEnergy",
                                   "Max Photon Energy Distribution", "E_{#gamma}^{max} [GeV]",
                                   (fs::path(preselectionDir) / "photon_energy.pdf").string(),
#if APPLY_PRE_HIGH_E_PHOTON_VETO
                                   PHOTON_ENERGY_CUT_GEV
#else
                                   -1
#endif
        );

        // 2. Число джетов с отметкой требования ровно 2 джета
        drawPreselectionHistograms(histsNJets, "nJets", "Number of Jets Distribution", "N_{jets}",
                                   (fs::path(preselectionDir) / "n_jets.pdf").string(),
#if APPLY_PRE_TWO_JETS_REQUIREMENT
                                   2
#else
                                   -1
#endif
        );

        // 3. Число конституентов с отметкой MIN_CONSTITUENTS_PER_JET
        drawPreselectionHistograms(histsNConstituents, "nConstituents",
                                   "Number of Constituents per Jet", "N_{constituents}",
                                   (fs::path(preselectionDir) / "n_constituents.pdf").string(),
#if APPLY_PRE_CONSTITUENTS_REQUIREMENT
                                   MIN_CONSTITUENTS_PER_JET
#else
                                   -1
#endif
        );

        std::cout << "Гистограммы предотборов сохранены в: " << preselectionDir << "\n";

        // =============================================================================
        // ОТРИСОВКА СУММАРНЫХ ГИСТОГРАММ ОСНОВНЫХ ОТБОРОВ
        // =============================================================================
        std::cout << "\nОтрисовка гистограмм основных отборов...\n";

        std::map<std::string, std::pair<TH1F *, ProcessInfo>> mainMETmap, mainDPhiMap, mainCosTMap,
            mainMjjMap, mainPmissMap;

        for (const auto &procEntry : mainSelHists) {
            const std::string &procName = procEntry.first;
            ProcessInfo info;
            auto it = processRecoilHists.find(procName);
            if (it != processRecoilHists.end())
                info = it->second.second;
            else
                info = ProcessInfo{procName, 1.0, kBlack, 1001};

            if (procEntry.second.count("met"))
                mainMETmap[procName] = {procEntry.second.at("met"), info};
            if (procEntry.second.count("deltaPhi"))
                mainDPhiMap[procName] = {procEntry.second.at("deltaPhi"), info};
            if (procEntry.second.count("cosThetaZ"))
                mainCosTMap[procName] = {procEntry.second.at("cosThetaZ"), info};
            if (procEntry.second.count("dijetMass"))
                mainMjjMap[procName] = {procEntry.second.at("dijetMass"), info};
            if (procEntry.second.count("pmiss"))
                mainPmissMap[procName] = {procEntry.second.at("pmiss"), info};
        }

        std::string mainSelDir = (fs::path(outputBaseDir) / "main_selection").string();
        fs::create_directories(mainSelDir);

        // MET: окно [min, max] 2 стрелки
        drawMainSelectionHistograms(mainMETmap, "MET Distribution", "MET_{jet} [GeV]",
                                    (fs::path(mainSelDir) / "main_met.pdf").string(),
#if APPLY_MAIN_MET_CUT
                                    { MET_CUT_MIN_GEV, MET_CUT_MAX_GEV }
#else
                                    {}
#endif
        );

        // deltaPhi: один верхний предел 1 стрелка
        drawMainSelectionHistograms(mainDPhiMap, "#Delta#phi Distribution", "#Delta#phi [rad]",
                                    (fs::path(mainSelDir) / "main_deltaPhi.pdf").string(),
#if APPLY_MAIN_DELTA_PHI_CUT
                                    { DELTA_PHI_CUT_MAX }
#else
                                    {}
#endif
        );

        // |cos theta_Z|: симметричный порог — 2 стрелки
        drawMainSelectionHistograms(mainCosTMap, "cos#theta_{Z} Distribution", "cos#theta_{Z}",
                                    (fs::path(mainSelDir) / "main_cosThetaZ.pdf").string(),
#if APPLY_MAIN_COS_THETA_Z_CUT
                                    { -COS_THETA_Z_CUT, COS_THETA_Z_CUT }
#else
                                    {}
#endif
        );

        // M_jj: окно [min, max] 2 стрелки
        drawMainSelectionHistograms(mainMjjMap, "Dijet Mass Distribution", "M_{jj} [GeV]",
                                    (fs::path(mainSelDir) / "main_dijetMass.pdf").string(),
#if APPLY_MAIN_DIJET_MASS_WINDOW
                                    { DIJET_MASS_WINDOW_MIN_GEV, DIJET_MASS_WINDOW_MAX_GEV }
#else
                                    {}
#endif
        );

        // Pmiss: окно [min, max] 2 стрелки
        drawMainSelectionHistograms(mainPmissMap, "Missing Momentum Distribution",
                                    "|P_{miss}| [GeV]",
                                    (fs::path(mainSelDir) / "main_pmiss.pdf").string(),
#if APPLY_MAIN_PMISS_CUT
                                    { PMISS_CUT_MIN_GEV, PMISS_CUT_MAX_GEV }
#else
                                    {}
#endif
        );

        std::cout << "Гистограммы основных отборов сохранены в: " << mainSelDir << "\n";
    }

    if (exportCsv && csvFile.is_open()) {
        csvFile.close();
        std::cout << "[CSV] Выгрузка завершена." << std::endl;
    }

    // Построение сравнительной гистограммы массы отдачи (qqHX и qqHinvi)
    std::string compOutput =
        (fs::path(outputBaseDir) / "recoil_comparison_qqHX_vs_signal.pdf").string();
    drawRecoilComparison(processRecoilHists, compOutput);

    // Построение стек гистограммы
    std::string stackOutput = (fs::path(outputBaseDir) / "recoil_stack_weighted.pdf").string();
    drawRecoilStack(processRecoilHists, RECOIL_STACK_ORDER, stackOutput);

    // Запуск шаблонного фита на накопленных данных
    if (runFit) {
        runMrecoilTemplateFit(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted,
                              vMrecoil_qqHX_Weighted, fs::path(outputBaseDir).string());
    }

    if (runAnalyticalFit) {
        runMrecoilAnalyticalFit(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted,
                                vMrecoil_qqHX_Weighted, fs::path(outputBaseDir).string());
    }

    // Запуск сканирования параметров адаптивности и числа бинов
    if (runScanFitParams) {
        runMrecoilScanAdapt(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted, vMrecoil_qqHX_Weighted,
                            fs::path(outputBaseDir).string());
    }

    // Запуск сканирования по mu для определения чувствительности на 95% CL
    if (runScanMu) {
        runMrecoilScanMu(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted, vMrecoil_qqHX_Weighted,
                         fs::path(outputBaseDir).string());
    }

    if (runAnalyticalScanMu) {
        runMrecoilAnalyticalScanMu(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted,
                                   vMrecoil_qqHX_Weighted, fs::path(outputBaseDir).string());
    }

    // Очистка
    for (auto &p : processRecoilHists)
        delete p.second.first;

    for (auto &procEntry : preselectionHists) {
        for (auto &histEntry : procEntry.second) {
            delete histEntry.second;
        }
    }

    for (auto &procEntry : mainSelHists)
        for (auto &histEntry : procEntry.second)
            delete histEntry.second;

    if (h_booster)
        XGBoosterFree(h_booster);

    return 0;
}
