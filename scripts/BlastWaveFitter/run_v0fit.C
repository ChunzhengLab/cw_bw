// V0SigmaFit.C
// 拟合 Blast-Wave v0(pT) 测得 sigma_T^2 与 sigma_beta^2

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include "BlastWaveFitter.C"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"

// ——————————————————————————
// 全局：已知的 Blast-Wave 参数
const int kNspecies = 3;
const double kMass[kNspecies] = {0.13957, 0.49368, 0.93827};
////static constexpr double Tkin = 0.107;                                     // [GeV] 平均温度
////static constexpr double betaT = 0.855;                                    // 无量纲 平均流速
////static constexpr double n_flow = 1.310;                                   //
////static constexpr double rho2 = 0.0942;                                    //
////static constexpr double Rx = 8.29055;                                     //
////static constexpr double Ry = 10.0;                                        //
////static constexpr double A_const[kNspecies] = {108200, 8709.45, 12292.3};  // 来自对谱的fit的参数A

static constexpr double Tkin = 0.108;                                      // [GeV] 平均温度
static constexpr double betaT = 0.837;                                     // 无量纲 平均流速
static constexpr double n_flow = 1.0;                                      //
static constexpr double rho2 = 0.0854;                                     //
static constexpr double Rx = 8.305;                                        //
static constexpr double Ry = 10.0;                                         //
static constexpr double A_const[kNspecies] = {7815.43, 10806.9, 90246.1};  // 来自对谱的fit的参数A
std::array<double, 3> gIntegralN = {0.0, 0.0, 0.0};                        // 对pt积分的粒子的个数
std::array<double, 3> gMT = {0.0, 0.0, 0.0};                               // MT    N'T(pt)dpt
std::array<double, 3> gMB = {0.0, 0.0, 0.0};                               // MB    N'Beta(pt)dpt
//// std::array<double, 3> gIntegralPt = {0.0, 0.0, 0.0};                       // 每种粒子的<[pT]>
//// std::array<double, 3> gDT = {0.0, 0.0, 0.0};
//// std::array<double, 3> gDB = {0.0, 0.0, 0.0};
//// double gSumN = 0.0;
//// double gSumDT = 0.0;
//// double gSumDB = 0.0;
//// double gSumMT = 0.0;
//// double gSumMB = 0.0;
//// double gSumIntegralPt = 0.0;
double gA = 0.0;
double gB = 0.0;
// 计算积分 A 与 B 的 pT 范围
static constexpr double kPTMin = 0.2;  /// 影响不大
static constexpr double kPTMax = 7.0;

/// // 全局变量：积分得到的 A 和 B
/// std::array<double, 3> gIntegralA = {0.0, 0.0, 0.0};
/// std::array<double, 3> gIntegralB = {0.0, 0.0, 0.0};
///  拟合范围
static constexpr double PTMIN_SPEC_PION_V0 = 0.5;
static constexpr double PTMAX_SPEC_PION_V0 = 1.5;
static constexpr double PTMIN_SPEC_KAON_V0 = 0.5;
static constexpr double PTMAX_SPEC_KAON_V0 = 2.0;
static constexpr double PTMIN_SPEC_PROTON_V0 = 0.5;
static constexpr double PTMAX_SPEC_PROTON_V0 = 2.5;

using namespace ROOT::Minuit2;

// ——————————————————————————
/// Step 0: 数值计算 DT DB MT MB Ni gIntegralPt
/// DT ≡ ∫ pT * ∂_T N(pT) dpT
/// DB ≡ ∫ pT * ∂_β N(pT) dpT
/// MT = ∫ ∂_T N(pT) dpT
/// MB = ∫ ∂_β N(pT) dpT
/// Ni = ∫ Ni(pT) dpT
/// gIntegralPt = ∫ pT*Ni(pT) dpT
////void Compute_IntegralConst(int i) {
////  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
////  ROOT::Math::Integrator intB(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
////  ROOT::Math::Integrator intC(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
////  ROOT::Math::Integrator intD(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
////
////  auto lamA = [&](double x) { return x * BW_dSpectrum_dT(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
////  auto lamB = [&](double x) { return x * BW_dSpectrum_dBeta(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
////  auto lamC = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry); };      // Ni
////  auto lamD = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry); };  // Ni
////
////  // 包装成 Functor，第二个参数 “1” 表示一维函数
////  std::function<double(double)> funA(lamA);
////  std::function<double(double)> funB(lamB);
////  std::function<double(double)> funC(lamC);
////  std::function<double(double)> funD(lamD);
////
////  // 再传给积分器
////  intA.SetFunction(funA);
////  intB.SetFunction(funB);
////  intC.SetFunction(funC);
////  intD.SetFunction(funD);
////
////  gIntegralN[i] = intC.Integral(kPTMin, kPTMax);
////  gIntegralPt[i] = intD.Integral(kPTMin, kPTMax);
////  gDT[i] = intA.Integral(kPTMin, kPTMax);
////  gDB[i] = intB.Integral(kPTMin, kPTMax);
////  gMT[i] = BW_MT(kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]);
////  gMB[i] = BW_MB(kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]);
////  std::cout << ">>> Computed integrals: A=" << gIntegralA[i] << "  B=" << gIntegralB[i] << "  N_i =" << gIntegralN[i] < < "  MT=" << gMT[i] << "  MB=" << gMB[i] << "\n";
////}
/////// Step 1: 计算常数 A B
////void Compute_AB() {
////  for (int i = 0; i < kNspecies; ++i) { Compute_IntegralConst(i); }
////  for (int i = 0; i < kNspecies; ++i) {
////    gSumN += gIntegralN[i];
////    gSumDT += gDT[i];
////    gSumDB += gDB[i];
////    gSumMT += gMT[i];
////    gSumMB += gMB[i];
////    gSumIntegralPt += gIntegralPt[i];
////  }
////  gA = gSumDT / gSumN - gSumIntegralPt * gSumMT / pow(gSumN, 2);
////  gB = gSumDB / gSumN - gSumIntegralPt * gSumMB / pow(gSumN, 2);
////}
// DT ≡ ∫ pT * ∂_T N(pT) dpT
// DB ≡ ∫ pT * ∂_β N(pT) dpT
// MT = ∫ ∂_T N(pT) dpT
// MB = ∫ ∂_β N(pT) dpT
// Ni = ∫ Ni(pT) dpT
// gIntegralPt = ∫ pT*Ni(pT) dpT
struct IntegralResults {
  double Ni, IPT, DT, DB, MT, MB;
};

IntegralResults ComputeOneSpecies(int i) {
  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intB(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intC(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intD(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);

  // 具体的 lambda 取决于你的 BW_* 函数签名
  std::function<double(double)> lamA = [&](double x) { return x * BW_dSpectrum_dT(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
  std::function<double(double)> lamB = [&](double x) { return x * BW_dSpectrum_dBeta(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
  std::function<double(double)> lamC = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry); };
  std::function<double(double)> lamD = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry); };

  intA.SetFunction(lamA);
  intB.SetFunction(lamB);
  intC.SetFunction(lamC);
  intD.SetFunction(lamD);

  IntegralResults R;
  R.Ni = intC.Integral(kPTMin, kPTMax);
  R.IPT = intD.Integral(kPTMin, kPTMax);
  R.DT = intA.Integral(kPTMin, kPTMax);
  R.DB = intB.Integral(kPTMin, kPTMax);
  R.MT = BW_MT(kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]);
  R.MB = BW_MB(kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]);
  // std::cout << "Species " << i << "  Ni=" << R.Ni << "  IPT=" << R.IPT << "  DT=" << R.DT << "  DB=" << R.DB << "  MT=" << R.MT << "  MB=" << R.MB << std::endl;
  return R;
}

void Compute_AB() {
  // 累加变量
  double sumNi = 0, sumIPT = 0, sumDT = 0, sumDB = 0, sumMT = 0, sumMB = 0;

  // 对每个 species 一次循环完成积分和累加
  for (int i = 0; i < kNspecies; ++i) {
    IntegralResults R = ComputeOneSpecies(i);
    sumNi += R.Ni;
    sumIPT += R.IPT;
    sumDT += R.DT;
    sumDB += R.DB;
    sumMT += R.MT;
    sumMB += R.MB;
    gIntegralN[i] = R.Ni;
    gMT[i] = R.MT;
    gMB[i] = R.MB;
    std::cout << " Species " << i << "  Ni=" << R.Ni << "  IPT=" << R.IPT << "  DT=" << R.DT << "  DB=" << R.DB << "  MT=" << R.MT << "  MB=" << R.MB << std::endl;
  }

  // 最终的 A, B
  gA = sumDT / sumNi - sumIPT * sumMT / (sumNi * sumNi);
  gB = sumDB / sumNi - sumIPT * sumMB / (sumNi * sumNi);

  std::cout << "Computed global parameters:\n"
            << "  gA = " << gA << "\n"
            << "  gB = " << gB << "\n";
}

// ——————————————————————————
// 全局指针：指向实验输入的 v0(pT) 带误差的 TGraph
TFile* fMain = TFile::Open("finalExp_3040.root");
static TGraph* gExpV0_pion = (TGraph*)fMain->Get("h_v0_pion");
static TGraph* gExpV0_kaon = (TGraph*)fMain->Get("h_v0_kaon");
static TGraph* gExpV0_proton = (TGraph*)fMain->Get("h_v0_proton");

// ——————————————————————————
// 拟合函数：Minuit2 用于计算 chi2
class V0FitterFunc : public ROOT::Minuit2::FCNBase {
 public:
  // par[0] = sigmaT2, par[1] = sigmaBeta2
  double operator()(const std::vector<double>& par) const override {
    double sigmaT2 = par[0];
    double sigmaBeta2 = par[1];
    double chi2 = 0.0;

    // 遍历所有数据点
    // Spectrum chi2 for pion
    for (int i = 0; i < gExpV0_pion->GetN(); ++i) {
      double pT, v0meas, err_lo, err_hi;
      gExpV0_pion->GetPoint(i, pT, v0meas);
      if (pT < PTMIN_SPEC_PION_V0 || pT > PTMAX_SPEC_PION_V0) continue;
      // err_lo = gExpV0_pion->GetErrorYlow(i);
      // err_hi = gExpV0_pion->GetErrorYhigh(i);
      err_lo = 0.005;
      err_hi = 0.005;

      double v0err = 0.5 * (err_lo + err_hi);

      // 计算模型 v0 模板
      double v0model = BW_V0_pT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0], gA, gB, sigmaT2, sigmaBeta2, gMT[0], gMB[0], gIntegralN[0]);
      //  double Np = A_const[0] * BW_Spectrum(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry);          // N(pT)
      //  double dNdT = BW_dSpectrum_dT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);     // ∂ₜ N
      //  double dNdB = BW_dSpectrum_dBeta(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);  // ∂_β N
      //  double numerator = gIntegralA[0] * dNdT * sigmaT2 + gIntegralB[0] * dNdB * sigmaBeta2;
      //  double denominator = Np * std::sqrt(gIntegralA[0] * gIntegralA[0] * sigmaT2 + gIntegralB[0] * gIntegralB[0] * sigmaBeta2);
      //  double v0model = numerator / denominator;

      chi2 += std::pow((v0meas - v0model) / v0err, 2);
    }
    // Spectrum chi2 for kaon
    for (int i = 0; i < gExpV0_kaon->GetN(); ++i) {
      double pT, v0meas, err_lo, err_hi;
      gExpV0_kaon->GetPoint(i, pT, v0meas);
      if (pT < PTMIN_SPEC_KAON_V0 || pT > PTMAX_SPEC_KAON_V0) continue;
      // err_lo = gExpV0_kaon->GetErrorYlow(i);
      // err_hi = gExpV0_kaon->GetErrorYhigh(i);
      err_lo = 0.005;
      err_hi = 0.005;
      double v0err = 0.5 * (err_lo + err_hi);

      // 计算模型 v0 模板
      double v0model = BW_V0_pT(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1], gA, gB, sigmaT2, sigmaBeta2, gMT[1], gMB[1], gIntegralN[1]);
      //  double Np = A_const[1] * BW_Spectrum(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry);          // N(pT)
      //  double dNdT = BW_dSpectrum_dT(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1]);     // ∂ₜ N
      //  double dNdB = BW_dSpectrum_dBeta(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1]);  // ∂_β N
      //  double numerator = gIntegralA[1] * dNdT * sigmaT2 + gIntegralB[1] * dNdB * sigmaBeta2;
      //  double denominator = Np * std::sqrt(gIntegralA[1] * gIntegralA[1] * sigmaT2 + gIntegralB[1] * gIntegralB[1] * sigmaBeta2);
      //  double v0model = numerator / denominator;

      chi2 += std::pow((v0meas - v0model) / v0err, 2);
    }
    // Spectrum chi2 for proton
    for (int i = 0; i < gExpV0_proton->GetN(); ++i) {
      double pT, v0meas, err_lo, err_hi;
      gExpV0_proton->GetPoint(i, pT, v0meas);
      if (pT < PTMIN_SPEC_PROTON_V0 || pT > PTMAX_SPEC_PROTON_V0) continue;
      // err_lo = gExpV0_proton->GetErrorYlow(i);
      // err_hi = gExpV0_proton->GetErrorYhigh(i);
      err_lo = 0.005;
      err_hi = 0.005;
      double v0err = 0.5 * (err_lo + err_hi);

      // 计算模型 v0 模板
      double v0model = BW_V0_pT(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2], gA, gB, sigmaT2, sigmaBeta2, gMT[2], gMB[2], gIntegralN[2]);

      //  double Np = A_const[2] * BW_Spectrum(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry);          // N(pT)
      //  double dNdT = BW_dSpectrum_dT(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2]);     // ∂ₜ N
      //  double dNdB = BW_dSpectrum_dBeta(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2]);  // ∂_β N
      //  double numerator = gIntegralA[2] * dNdT * sigmaT2 + gIntegralB[2] * dNdB * sigmaBeta2;
      //  double denominator = Np * std::sqrt(gIntegralA[2] * gIntegralA[2] * sigmaT2 + gIntegralB[2] * gIntegralB[2] * sigmaBeta2);
      //  double v0model = numerator / denominator;

      chi2 += std::pow((v0meas - v0model) / v0err, 2);
    }
    return chi2;
  }

  double Up() const override {
    return 1.0;  // chi2 的单位
  }
};

// ——————————————————————————
// 主执行接口
void run_v0fit() {
  if (!gExpV0_pion) {
    std::cerr << "[Error] Experimental V0 graph not set! Call SetExperimentalV0(...) first.\n";
    return;
  }

  // 1) 先算 A,B
  // for (int i = 0; i < kNspecies; ++i) { Compute_IntegralAB(i); }
  Compute_AB();

  // 2) 初始化 Minuit2 参数
  ROOT::Minuit2::MnUserParameters upp;
  // upp.Add("sigmaT2", 0.0, 1e-5, 0.0, 1e-2);
  // upp.Fix("sigmaT2");
  upp.Add("sigmaT2", 1e-6, 1e-6, 0.0, 1e-2);
  upp.Add("sigmaBeta2", 1e-4, 1e-5, 0.0, 1e-2);

  // 3) 最小化
  V0FitterFunc fitter;
  ROOT::Minuit2::MnMigrad migrad(fitter, upp);
  auto result = migrad();

  // 4) 打印与检验
  bool converged = result.IsValid();
  double chi2 = result.Fval();
  int ndf = gExpV0_pion->GetN() - 2;
  std::cout << "\n=== V0(Sigma^2) Fit Results ===\n"
            << " Converged: " << converged << "  chi2/ndf = " << chi2 << "/" << ndf << " = " << chi2 / ndf << "\n";

  double sigmaT2 = result.UserState().Value(0);
  double sigmaT2_err = result.UserState().Error(0);
  double sigmaB2 = result.UserState().Value(1);
  double sigmaB2_err = result.UserState().Error(1);

  std::cout << Form(
      "  sigma_T^2    = %.3e ± %.3e\n"
      "  sigma_Beta^2 = %.3e ± %.3e\n",
      sigmaT2, sigmaT2_err, sigmaB2, sigmaB2_err);

  // 5) 画图对比：带实验点与模型曲线
  TCanvas c("cV0Fit", "v0(pT) Fit", 800, 600);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gExpV0_pion->SetMarkerStyle(21);
  gExpV0_pion->GetYaxis()->SetRangeUser(-0.08, 0.35);
  gExpV0_pion->SetTitle("");
  gExpV0_pion->GetXaxis()->SetTitle("#it{P}_{T} (Gev/c)");
  gExpV0_pion->GetYaxis()->SetTitle("#it{v}_{0}(#it{p}_{T})");
  gExpV0_pion->GetYaxis()->SetTitleSize(0.056);
  gExpV0_pion->GetXaxis()->SetTitleSize(0.056);
  gExpV0_pion->GetYaxis()->SetTitleOffset(1.05);
  gExpV0_pion->GetXaxis()->SetLabelSize(0.055);
  gExpV0_pion->GetYaxis()->SetLabelSize(0.045);
  gExpV0_pion->GetYaxis()->SetNdivisions(410);
  gExpV0_pion->GetXaxis()->SetNdivisions(410);
  gExpV0_pion->SetMarkerSize(1.5);
  gExpV0_pion->SetMarkerColor(kGreen + 2);
  gExpV0_kaon->SetMarkerSize(1.5);
  gExpV0_kaon->SetMarkerColor(kRed);
  gExpV0_kaon->SetMarkerStyle(22);
  gExpV0_proton->SetMarkerSize(1.5);
  gExpV0_proton->SetMarkerColor(kBlue);
  gExpV0_pion->Draw("AP");
  gExpV0_kaon->Draw("P SAME");
  gExpV0_proton->Draw("P SAME");

  // pion
  TF1 fModel_pion(
      "fV0Model_pion",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_pT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0], gA, gB, sT2, sB2, gMT[0], gMB[0], gIntegralN[0]);
        // double Np = A_const[0] * BW_Spectrum(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry);
        // double dNdT = BW_dSpectrum_dT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);
        // double dNdB = BW_dSpectrum_dBeta(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);
        // double num = gIntegralA[0] * dNdT * sT2 + gIntegralB[0] * dNdB * sB2;
        // double den = Np * std::sqrt(gIntegralA[0] * gIntegralA[0] * sT2 + gIntegralB[0] * gIntegralB[0] * sB2);
        // return num / den;
      },
      0.5, 6., /*Npar=*/2);

  fModel_pion.SetParameters(sigmaT2, sigmaB2);
  fModel_pion.SetLineColor(kGreen + 2);
  fModel_pion.SetLineWidth(2);
  fModel_pion.Draw("Same");
  /// kaon
  TF1 fModel_kaon(
      "fV0Model_kaon",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_pT(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1], gA, gB, sT2, sB2, gMT[1], gMB[1], gIntegralN[1]);
        //  double Np = A_const[1] * BW_Spectrum(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry);
        //  double dNdT = BW_dSpectrum_dT(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1]);
        //  double dNdB = BW_dSpectrum_dBeta(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1]);
        //  double num = gIntegralA[1] * dNdT * sT2 + gIntegralB[1] * dNdB * sB2;
        //  double den = Np * std::sqrt(gIntegralA[1] * gIntegralA[1] * sT2 + gIntegralB[1] * gIntegralB[1] * sB2);
        //  return num / den;
      },
      0.5, 6., /*Npar=*/2);

  fModel_kaon.SetParameters(sigmaT2, sigmaB2);
  fModel_kaon.SetLineColor(kRed);
  fModel_kaon.SetLineWidth(2);
  fModel_kaon.Draw("Same");
  // proton
  TF1 fModel_proton(
      "fV0Model",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_pT(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2], gA, gB, sT2, sB2, gMT[2], gMB[2], gIntegralN[2]);
        // double Np = A_const[2] * BW_Spectrum(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry);
        // double dNdT = BW_dSpectrum_dT(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2]);
        // double dNdB = BW_dSpectrum_dBeta(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2]);
        // double num = gIntegralA[2] * dNdT * sT2 + gIntegralB[2] * dNdB * sB2;
        // double den = Np * std::sqrt(gIntegralA[2] * gIntegralA[2] * sT2 + gIntegralB[2] * gIntegralB[2] * sB2);
        // return num / den;
      },
      0.5, 6., /*Npar=*/2);

  fModel_proton.SetParameters(sigmaT2, sigmaB2);
  fModel_proton.SetLineColor(kBlue);
  fModel_proton.SetLineWidth(2);
  fModel_proton.Draw("Same");

  auto leg = new TLegend(0.15, 0.65, 0.51, 0.88);
  leg->SetFillStyle(0);
  leg->SetHeader("ALICE     BW fit");
  leg->SetNColumns(2);

  leg->AddEntry(gExpV0_pion, "       ", "p");
  leg->AddEntry(&fModel_pion, "#pi^{#pm}", "l");
  leg->AddEntry(gExpV0_kaon, "       ", "p");
  leg->AddEntry(&fModel_kaon, "k^{#pm}", "l");
  leg->AddEntry(gExpV0_proton, "       ", "p");
  leg->AddEntry(&fModel_proton, "p (#bar{p})", "l");
  leg->SetBorderSize(0);
  leg->Draw("same");

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);

  c.SaveAs("V0SigmaFit.pdf");
}
