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
/// static constexpr double Tkin = 0.107;                                     // [GeV] 平均温度
/// static constexpr double betaT = 0.855;                                    // 无量纲 平均流速
///                                                                           // static constexpr double n_flow = 1.0;
/// static constexpr double n_flow = 1.310;                                   //                                                                        //
/// static constexpr double rho2 = 0.0942;                                    //
/// static constexpr double Rx = 8.29055;                                     //
/// static constexpr double Ry = 10.0;                                        //
/// static constexpr double A_const[kNspecies] = {108200, 8709.45, 12292.3};  //

static constexpr double Tkin = 0.108;                                      // [GeV] 平均温度
static constexpr double betaT = 0.837;                                     // 无量纲 平均流速
static constexpr double n_flow = 1.0;                                      //
static constexpr double rho2 = 0.0854;                                     //
static constexpr double Rx = 8.305;                                        //
static constexpr double Ry = 10.0;                                         //
static constexpr double A_const[kNspecies] = {90246.1, 7815.43, 10806.9};  // 来自对谱的fit的参数A
std::array<double, 3> gIntegralN = {0.0, 0.0, 0.0};                        // 对pt积分的粒子的个数

// 计算积分 A 与 B 的 pT 范围
static constexpr double kPTMin = 0.2;
static constexpr double kPTMax = 7.0;

// 全局变量：积分得到的 A 和 B
std::array<double, 3> gIntegralA = {0.0, 0.0, 0.0};
std::array<double, 3> gIntegralB = {0.0, 0.0, 0.0};

// ——————————————————————————
/// Step 1: 数值计算 A 与 B
/// A ≡ ∫ pT * ∂_T N(pT) dpT
/// B ≡ ∫ pT * ∂_β N(pT) dpT
void Compute_IntegralAB(int i, double i_betaT) {
  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intB(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  ROOT::Math::Integrator intC(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);

  auto lamA = [&](double x) { return x * BW_dSpectrum_dT(x, kMass[i], Tkin, i_betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
  auto lamB = [&](double x) { return x * BW_dSpectrum_dBeta(x, kMass[i], Tkin, i_betaT, n_flow, rho2, Rx, Ry, A_const[i]); };
  auto lamC = [&](double x) { return A_const[i] * 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry); };  // Ni

  // 包装成 Functor，第二个参数 “1” 表示一维函数
  std::function<double(double)> funA(lamA);
  std::function<double(double)> funB(lamB);
  std::function<double(double)> funC(lamC);

  // 再传给积分器
  intA.SetFunction(funA);
  intB.SetFunction(funB);
  intC.SetFunction(funC);

  gIntegralN[i] = intC.Integral(kPTMin, kPTMax);
  gIntegralA[i] = intA.Integral(kPTMin, kPTMax) / gIntegralN[i];
  gIntegralB[i] = intB.Integral(kPTMin, kPTMax) / gIntegralN[i];
  std::cout << "betaT = " << i_betaT << ">>> Computed integrals: A=" << gIntegralA[i] << "  B=" << gIntegralB[i] << "\n";
}

void test_rum() {
  TCanvas* c = new TCanvas("c1", "v0 Model", 800, 600);
  c->SetGrid();
  c->cd();
  const int iCase = 0;
  TF1 fModel_proton(
      "fV0Model",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        double Np = A_const[0] * BW_Spectrum(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry);
        double dNdT = BW_dSpectrum_dT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);
        double dNdB = BW_dSpectrum_dBeta(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);
        double num = gIntegralA[0] * dNdT * sT2 + gIntegralB[0] * dNdB * sB2;
        double den = Np * std::sqrt(gIntegralA[0] * gIntegralA[0] * sT2 + gIntegralB[0] * gIntegralB[0] * sB2);
        return num / (den);
      },
      0., 10., /*Npar=*/2);
  fModel_proton.SetParameters(0.00, 4.975e-05);
  fModel_proton.SetLineColor(kBlue);
  fModel_proton.SetLineWidth(2);
  // fModel_proton.Draw();
  //  TF1 f("f", "x", 1e-3, 5.);
  //  f.SetLineColor(kBlue);
  //  f.Draw();

  // TLegend leg(0.62, 0.65, 0.88, 0.88);
  // leg.AddEntry(&fModel_proton, "Model", "l");
  // leg.Draw("same");

  TF1 fModel_proton_slope(
      "fModel_proton_slope",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        double Np = A_const[0] * BW_Spectrum(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry);
        double dNdB = BW_dSpectrum_dBeta(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0]);
        double slope_A = sqrt(sB2) * dNdB * (1 - betaT * betaT) * Tkin / pT;
        double slope_B = Np;
        double slope = slope_A / (slope_B * (1 - betaT * betaT) * Tkin);
        return slope;
      },
      1., 10., /*Npar=*/2);
  fModel_proton_slope.SetParameters(0.00, 4.975e-05);
  fModel_proton_slope.SetLineColor(kBlue);
  fModel_proton_slope.SetLineWidth(2);
  // fModel_proton_slope.Draw();

  TF1* f_v0_pion = new TF1(
      "f_v0_pion",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_pT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0], gIntegralA, gIntegralB, sT2, sB2);
      },
      0.,  // pT 下限
      5.,  // pT 上限
      2);
  f_v0_pion->SetParameters(0.00, 2.330e-04);
  f_v0_pion->SetLineColor(kBlue);
  f_v0_pion->SetLineWidth(2);
  f_v0_pion->Draw();

  /////////For slop
  Compute_IntegralAB(iCase, 0.837);
  TF1* f_v0_pion_slope = new TF1(
      "f_v0_pion_slope",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_dV0_dpT(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0], gIntegralA, gIntegralB, sT2, sB2);
      },
      0.2,  // pT 下限
      5,    // pT 上限
      2);
  f_v0_pion_slope->SetParameters(0.00, 2.330e-04);
  f_v0_pion_slope->SetLineColor(kBlue);
  f_v0_pion_slope->SetLineWidth(2);
  f_v0_pion_slope->SetMinimum(0.04);  // y 轴下限
  f_v0_pion_slope->SetMaximum(0.08);  // y 轴上限
                                      // f_v0_pion_slope->Draw("L");
  ////with different mass
  Compute_IntegralAB(iCase, 0.5);
  TF1* f_v0_pion_slope_2 = new TF1(
      "f_v0_pion_slope_2",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_dV0_dpT(pT, kMass[0], Tkin, 0.30, n_flow, rho2, Rx, Ry, A_const[0], gIntegralA, gIntegralB, sT2, sB2);
      },
      0.5,  // pT 下限
      5,    // pT 上限
      2);
  f_v0_pion_slope_2->SetParameters(0.00, 2.330e-04);
  f_v0_pion_slope_2->SetLineColor(kRed);
  f_v0_pion_slope_2->SetLineWidth(2);
  // f_v0_pion_slope_2->Draw("L same");

  TLegend leg(0.62, 0.65, 0.88, 0.88);
  leg.AddEntry(f_v0_pion_slope, "#beta_{T} = 0.855", "l");
  leg.AddEntry(f_v0_pion_slope_2, "#beta_{T} = 0.8", "l");
  leg.Draw("same");

  cout << "slope1 @ pT=2.0 = " << f_v0_pion_slope->Eval(5.0) << endl;
  cout << "slope2 @ pT=2.0 = " << f_v0_pion_slope_2->Eval(5.0) << endl;
  c->SaveAs("test.pdf");
}
