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
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooKeysPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <TCanvas.h>
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
    TLegend *leg = new TLegend(0.75, 0.65, 0.95, 0.88);
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

            // Применяем стиль и добавляем в стек
            hForStack->SetFillColor(info.color);
            hForStack->SetFillStyle(info.fillStyle);
            hForStack->SetMarkerStyle(21);
            hForStack->SetMarkerColor(info.color);
            hForStack->SetLineWidth(1);
            stack->Add(hForStack);
            leg->AddEntry(hForStack, info.legendName.c_str(), "f");

            // Клонируем сигнал чтобы отрисовать его от низа стековой гистограммы и настраиваем
            // параметры отображения
            if (info.legendName.find("signal") != std::string::npos) {
                signalHistBottom = (TH1F *)hForStack->Clone("signal_bottom");
                signalHistBottom->SetFillStyle(0);
                signalHistBottom->SetLineWidth(3);
                signalHistBottom->SetLineColor(info.color);
                leg->AddEntry(signalHistBottom, "Signal (overlay)", "L");
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
                           const std::string &outputPath) {

    if (vSignal.empty() || vBkg.empty()) {
        std::cerr << "[Fit] Ошибка: входные данные пусты!\n";
        return;
    }

    std::cout << "\n[Fit] =====================================================\n";
    std::cout << "[Fit] Запуск шаблонного фита M_recoil\n";
    std::cout << "[Fit] Сигнальных MC-событий: " << vSignal.size() << "\n";
    std::cout << "[Fit] Фоновых MC-событий:    " << vBkg.size() << "\n";
    std::cout << "[Fit] mu (доля сигнала):     " << FIT_PSEUDO_MU << "\n";

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

    // Суммарные взвешенные числа событий в шаблонах.
    // sumEntries() с весами возвращает сумму весов всех событий.
    // Это ожидаемое число событий при данной светимости.
    double sumW_S = dsSignal->sumEntries();
    double sumW_B = dsBkg->sumEntries();

    std::cout << "[Fit] sumW_S (ожидаемый сигнал, полный):  " << sumW_S << "\n";
    std::cout << "[Fit] sumW_B (ожидаемый фон):             " << sumW_B << "\n";
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
    // ШАГ 4: ПАРАМЕТРЫ НОРМИРОВКИ (число событий)
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
    // ШАГ 5: КОМБИНИРОВАННАЯ МОДЕЛЬ
    // =========================================================================
    //
    // RooAddPdf с двумя аргументами-нормировками создаёт расширенный PDF:
    //   model(x) = nS * fs(x) + nB * fb(x)
    //
    RooAddPdf model("model", "nS*fs + nB*fb", RooArgList(pdfSignal, pdfBkg), RooArgList(nS, nB));
    model.setNormRange("fitRange");

    // =========================================================================
    // ШАГ 6: ПСЕВДОДАННЫЕ
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
    // ШАГ 7: ФИТ 1. НУЛЕВАЯ ГИПОТЕЗА H0 (только фон, nS = 0)
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
    // ШАГ 8: ФИТ 2. АЛЬТЕРНАТИВНАЯ ГИПОТЕЗА H1 (сигнал + фон)
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
    // ШАГ 9: ОЦЕНКА ЗНАЧИМОСТИ
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
    double significance_lrt = std::sqrt(std::max(0.0, q0));

    std::cout << "\n[Fit] =====================================================\n";
    std::cout << "[Fit] РЕЗУЛЬТАТЫ ФИТА (mu = " << FIT_PSEUDO_MU << ")\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   nS (fitted)  = " << fit_nS << " ± " << fit_nS_err << "\n";
    std::cout << "[Fit]   nB (fitted)  = " << fit_nB << " ± " << fit_nB_err << "\n";
    std::cout << "[Fit]   nS (expected)= " << expectedNS << "\n";
    std::cout << "[Fit]   nB (expected)= " << expectedNB << "\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   NLL (H0, фон):      " << nll_b << "\n";
    std::cout << "[Fit]   NLL (H1, сигн+фон): " << nll_sb << "\n";
    std::cout << "[Fit]   q0 = 2*(NLL_b - NLL_sb) = " << q0 << "\n";
    std::cout << "[Fit] -----------------------------------------------------\n";
    std::cout << "[Fit]   Значимость (простая):  Z = " << significance_approx << " σ\n";
    std::cout << "[Fit]   Значимость (LRT):      Z = " << significance_lrt << " σ\n";
    std::cout << "[Fit] =====================================================\n\n";

    // =========================================================================
    // ШАГ 10: ВИЗУАЛИЗАЦИЯ
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
    TPaveText *info = new TPaveText(0.76, 0.7, 0.98, 0.83, "NDC NB");
    info->SetFillColor(0);
    info->SetFillStyle(1001);
    info->SetBorderSize(1);
    info->SetTextAlign(12);
    info->SetTextSize(0.020);
    info->AddText(Form("mu = %.2f", FIT_PSEUDO_MU));
    info->AddText(Form("nS = %.1f #pm %.1f", fit_nS, fit_nS_err));
    info->AddText(Form("nB = %.1f #pm %.1f", fit_nB, fit_nB_err));
    info->AddText(Form("Z (simple) = %.2f #sigma", significance_approx));
    info->AddText(Form("Z (LRT)    = %.2f #sigma", significance_lrt));
    info->Draw();

    cFit->SaveAs(outputPath.c_str());
    std::cout << "[Fit] График сохранён: " << outputPath << "\n";

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
    std::cout << "Многофайловый анализ ee -> ZH -> qq + invisible\n\n"
              << "Использование:\n"
              << "  " << progName << " file1.root file2.root ... [options]\n\n"
              << "Опции:\n"
              << "  -h, --help              Показать эту справку\n"
              << "  -o, --output-dir DIR    Базовая директория результатов (по умолчанию: "
                 "../pdf_results)\n"
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

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 < argc)
                outputBaseDir = argv[++i];
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

    // Контейнер для отдельных гистограмм каждого процесса + их метаданные
    std::map<std::string, std::pair<TH1F *, ProcessInfo>> processRecoilHists;

    // Контейнеры для накопления Mrecoil после всех отборов
    std::vector<std::pair<double, double>> vMrecoil_Signal_Weighted;
    std::vector<std::pair<double, double>> vMrecoil_Bkg_Weighted;

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

            // ==================== ОСНОВНЫЕ ОТБОРЫ ====================
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
                !isInsideEllipse(invMass, recoilMass, ELLIPSE_CX_GEV, ELLIPSE_CY_GEV, ELLIPSE_A_GEV,
                                 ELLIPSE_B_GEV, ELLIPSE_THETA))
                continue;
            stats.afterEllipseCut++;

            stats.finalSelected++;

            // Если событие прошло все отборы, то заполняем гистограмму текущего процесса. Веса
            // применим уже при построении стековой гистограммы
            hRecoilMassWeight->Fill(recoilMass);

            // Определяем, сигнал это или фон (по имени процесса). Заполняем соотвествующий
            // контейнер с весом.
            bool isSignal = (processName.find("qqHinvi") != std::string::npos);
            if (isSignal)
                vMrecoil_Signal_Weighted.emplace_back(recoilMass, proc.weight);
            else
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

        std::vector<std::pair<double, std::string>> deltaPhiMarks;
#if APPLY_MAIN_DELTA_PHI_CUT
        deltaPhiMarks.emplace_back(DELTA_PHI_CUT_MAX, "#Delta#phi^{cut}");
#endif
        drawHistogram1D(hDeltaPhi, "cDeltaPhi", "#Delta#phi [rad]", makeOutputPath("deltaPhi_jets"),
                        deltaPhiMarks, kGreen + 2, 2);

        // Остальные гистограммы без линий отборов
        drawHistogram1D(hPmissMag, "cPmissMag", "|P_{miss}| [GeV]", OUTPUT_PMISS_MAG, {}, kOrange,
                        2);
        drawHistogram1D(hCosThetaPmiss, "cCosThetaPmiss", "cos#theta_{miss}",
                        OUTPUT_COS_THETA_PMISS, {}, kMagenta, 2);
        drawHistogram1D(hDijetEnergy, "cDijetEnergy", "E_{jj} [GeV]",
                        makeOutputPath("dijet_energy"), {}, kAzure + 1, 2);
        drawHistogram1D(hDeltaTheta, "cDeltaTheta", "#Delta#theta [rad]",
                        makeOutputPath("deltaTheta_jets"), {}, kOrange + 1, 2);

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

        std::cout << "\nГотово. Результаты сохранены в: " << fs::absolute(processOutputDir)
                  << std::endl;
    }

    // Построение сравнительной гистограммы массы отдачи (qqHX и qqHinvi)
    std::string compOutput =
        (fs::path(outputBaseDir) / "recoil_comparison_qqHX_vs_signal.pdf").string();
    drawRecoilComparison(processRecoilHists, compOutput);

    // Построение стек гистограммы
    std::string stackOutput = (fs::path(outputBaseDir) / "recoil_stack_weighted.pdf").string();
    drawRecoilStack(processRecoilHists, RECOIL_STACK_ORDER, stackOutput);

    // Запуск шаблонного фита на накопленных данных
    std::string fitOutput = (fs::path(outputBaseDir) / "template_fit_recoil.pdf").string();
    runMrecoilTemplateFit(vMrecoil_Signal_Weighted, vMrecoil_Bkg_Weighted, fitOutput);

    // Очистка
    for (auto &p : processRecoilHists)
        delete p.second.first;

    return 0;
}
