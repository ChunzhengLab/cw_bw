#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include "BlastWaveFitter.C"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"

// 全局：已知的 Blast-Wave 参数
const int kNspecies = 3;
const double kMass[kNspecies] = {0.13957, 0.49368, 0.93827};
static constexpr double Tkin = 0.108;                                      // [GeV] 平均温度
static constexpr double betaT = 0.837;                                     // 无量纲 平均流速
static constexpr double n_flow = 1.0;                                      //
static constexpr double rho2 = 0.0854;                                     //
static constexpr double Rx = 8.305;                                        //
static constexpr double Ry = 10.0;                                         //
static constexpr double A_const[kNspecies] = {7815.43, 10806.9, 90246.1};  // 来自对谱的fit的参数A
static constexpr int NClass = 3;
const double rho2_class[NClass] = {0.04, 0.08, 0.12};

std::array<std::array<double, 3>, NClass> gIntegralN{};  // 对pt积分的粒子的个数
std::array<std::array<double, 3>, NClass> gMT{};         // MT    N'T(pt)dpt
std::array<std::array<double, 3>, NClass> gMB{};         // MB    N'Beta(pt)dpt
// 全局变量：积分得到的 A 和 B
std::array<double, 6> gA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 6> gB = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
std::array<double, 6> gIntegralV2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// 计算积分 A 与 B 的 pT 范围
static constexpr double kPTMin = 0.2;  /// 影响不大
static constexpr double kPTMax = 7.0;

TF1* f_v0_pion[NClass];

struct IntegralResults {
  double Ni, IPT, DT, DB, MT, MB, IV2;
};
// ——————————————————————————
/// Step 1: 数值计算 A 与 B
/// A ≡ ∫ pT * ∂_T N(pT) dpT
/// B ≡ ∫ pT * ∂_β N(pT) dpT
IntegralResults ComputeOneSpecies(int i, double i_rho2) {
  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intB(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intC(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intD(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intE(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);

  // 具体的 lambda 取决于你的 BW_* 函数签名
  std::function<double(double)> lamA = [&](double x) { return x * BW_dSpectrum_dT(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry, A_const[i]); };
  std::function<double(double)> lamB = [&](double x) { return x * BW_dSpectrum_dBeta(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry, A_const[i]); };
  std::function<double(double)> lamC = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry); };
  std::function<double(double)> lamD = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry); };
  std::function<double(double)> lamE = [&](double x) {
    return A_const[i] * 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry) * BW_v2shape(x, kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry);
  };

  intA.SetFunction(lamA);
  intB.SetFunction(lamB);
  intC.SetFunction(lamC);
  intD.SetFunction(lamD);
  intE.SetFunction(lamE);

  IntegralResults R;
  R.Ni = intC.Integral(kPTMin, kPTMax);
  R.IPT = intD.Integral(kPTMin, kPTMax);
  R.DT = intA.Integral(kPTMin, kPTMax);
  R.DB = intB.Integral(kPTMin, kPTMax);
  R.MT = BW_MT(kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry, A_const[i]);
  R.MB = BW_MB(kMass[i], Tkin, betaT, n_flow, i_rho2, Rx, Ry, A_const[i]);
  R.IV2 = intE.Integral(kPTMin, kPTMax);
  // std::cout << "Species " << i << "  Ni=" << R.Ni << "  IPT=" << R.IPT << "  DT=" << R.DT << "  DB=" << R.DB << "  MT=" << R.MT << "  MB=" << R.MB << std::endl;
  return R;
}
void Compute_AB() {
  // 累加变量

  // 对每个 species 一次循环完成积分和累加
  for (int iclass = 0; iclass < NClass; ++iclass) {
    double sumNi = 0, sumIPT = 0, sumDT = 0, sumDB = 0, sumMT = 0, sumMB = 0;
    for (int i = 0; i < kNspecies; ++i) {
      IntegralResults R;
      if (i == 0) {
        R = ComputeOneSpecies(i, rho2_class[iclass]);
        gIntegralV2[iclass] = R.IV2 / R.Ni;
      } else {
        R = ComputeOneSpecies(i, rho2);
      }
      sumNi += R.Ni;
      sumIPT += R.IPT;
      sumDT += R.DT;
      sumDB += R.DB;
      sumMT += R.MT;
      sumMB += R.MB;
      gIntegralN[iclass][i] = R.Ni;
      gMT[iclass][i] = R.MT;
      gMB[iclass][i] = R.MB;
      std::cout << " Species " << i << "  Ni=" << R.Ni << "  IPT=" << R.IPT << "  DT=" << R.DT << "  DB=" << R.DB << "  MT=" << R.MT << "  MB=" << R.MB
                << " V2=" << gIntegralV2[iclass] << std::endl;
    }
    // 最终的 A, B
    gA[iclass] = sumDT / sumNi - sumIPT * sumMT / (sumNi * sumNi);
    gB[iclass] = sumDB / sumNi - sumIPT * sumMB / (sumNi * sumNi);
    //   std::cout << "Computed global parameters:\n"
    //             << "  mass_class = " << iclass << "\n"
    //             << "  gA = " << gA[iclass] << "\n"
    //             << "  gB = " << gB[iclass] << "\n";
  }
}

void run_v0pT_rho2() {
  // 1) 创建一个画布
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  auto c = new TCanvas("c_v0_10", "10 条 v0(pt) 曲线比较", 900, 700);
  // c->SetGrid();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  // 2) 先定义 pt 的范围、以及 slope、sT2、sB2 所需的参数
  const double ptMin = 0.2;
  const double ptMax = 5.0;
  // 假设我们只用 sT2 = 0.0, sB2 = 4.975e-05 这对值。如果每条线也有不同的 sT2、sB2，
  // 只需要把它们做成两个向量，然后在 lambda 里针对 index i 使用不同值。
  const double sT2_val = 1.411e-06;
  // const double sT2_val = 0.0;
  const double sB2_val = 9.579e-05;

  const int iCase = 0;

  // 3) 定义颜色调色板 —— 10 种颜色
  //    下面用 TColor::GetColorPalette() 也可以，或者手动挑 10 个 ROOT 预定义颜色
  int colors[10] = {kBlue, kRed, kGreen + 2, kMagenta, kCyan, kOrange + 1, kViolet + 1, kAzure + 3, kPink + 7, kSpring + 5};

  // 4) 先来一个“基准”TF1，让它把坐标轴范围定好，并且设置好 Y 轴上下限
  //    这里我们假设 y 在 -0.1 到 0.3 之间（根据实际情况调整）。如果不确定，可以
  //    按照第一个 mass (0.1) 画出来后，读数再手动修改 SetMinimum/SetMaximum。
  TH1D* f_base = new TH1D("f_base", "", 10, ptMin, ptMax);
  f_base->SetMinimum(-0.1);  // 下限根据实际数据调整
  f_base->SetMaximum(0.4);   // 上限根据实际数据调整
  f_base->SetTitle("");
  f_base->GetXaxis()->SetTitle("#it{P}_{T} (Gev/c)");
  f_base->GetYaxis()->SetTitle("#it{v}_{0}(#it{p}_{T})");
  f_base->GetYaxis()->SetTitleSize(0.056);
  f_base->GetXaxis()->SetTitleSize(0.056);
  f_base->GetYaxis()->SetTitleOffset(1.05);
  f_base->GetXaxis()->SetLabelSize(0.055);
  f_base->GetYaxis()->SetLabelSize(0.045);
  f_base->GetYaxis()->SetNdivisions(410);
  f_base->GetXaxis()->SetNdivisions(410);
  f_base->Draw();  // 只画坐标系

  // 5) 循环 10 个质量值，从 0.1 到 1.0
  Compute_AB();
  for (int idx = 0; idx < NClass; ++idx) {
    // 每一条用一个 lambda 把 v0 函数包起来
    // 用 std::function 包装 lambda，确保 ROOT 不会把它当“字符串公式”解析
    std::function<double(double*, double*)> fun_v0 = [=](double* x, double* p) {
      double pT = x[0];
      double sT2 = p[0];
      double sB2 = p[1];
      if (pT <= 1e-6) return 0.0;
      double v0 = BW_V0_pT(pT, kMass[0], Tkin, betaT, n_flow, rho2_class[idx], Rx, Ry, A_const[0], gA[idx], gB[idx], sT2, sB2, gMT[idx][0], gMB[idx][0], gIntegralN[idx][0]);
      return std::isfinite(v0) ? v0 : 0.0;
    };
    f_v0_pion[idx] = new TF1(Form("v0_mass_%.2f", rho2_class[idx]), fun_v0, ptMin, ptMax,
                             2  // 参数个数
    );
    f_v0_pion[idx]->SetParameters(sT2_val, sB2_val);
    // 画线：指定不同颜色，不用 "same"（因为我们已经画好了基准面）
    f_v0_pion[idx]->SetLineColor(colors[idx]);
    f_v0_pion[idx]->SetLineStyle(1);
    f_v0_pion[idx]->SetLineWidth(2);
  }
  for (int idx = 0; idx < NClass; ++idx) {
    cout << "v0  = " << f_v0_pion[idx]->Eval(1.0) << endl;
    // if (idx == 0) {
    //   f_v0_pion[0]->Draw("L");
    // } else {
    f_v0_pion[idx]->Draw("L same");
    // }
  }
  /////
  // // 6) 加一个图例，把 10 条质量标注出来
  auto leg = new TLegend(0.15, 0.55, 0.43, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int idx = 0; idx < NClass; ++idx) {
    TString entryName = Form("#rho_{2} = %.2f, v_{2} = %.3f", rho2_class[idx], gIntegralV2[idx]);
    // 注意：Legend 中用 TF1 的名字来索引
    leg->AddEntry(Form("v0_mass_%.2f", rho2_class[idx]), entryName.Data(), "l");
  }
  leg->Draw("same");

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);

  // 7) 最后保存成 PDF
  c->SaveAs("v0_10_rho2.pdf");
}
