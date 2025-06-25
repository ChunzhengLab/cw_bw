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
static constexpr double A_const[kNspecies] = {7815.43, 10806.9, 90246.1};  // 来自对谱的fit的参数A
std::array<double, 3> gIntegralN = {0.0, 0.0, 0.0};                        // 对pt积分的粒子的个数
std::array<double, 3> gMT = {0.0, 0.0, 0.0};                               // 对pt积分的粒子的个数

// 计算积分 A 与 B 的 pT 范围
static constexpr double kPTMin = 0.2;
static constexpr double kPTMax = 7.0;

// 全局变量：积分得到的 A 和 B
std::array<double, 3> gIntegralA = {0.0, 0.0, 0.0};
std::array<double, 3> gIntegralB = {0.0, 0.0, 0.0};

// ——————————————————————————
// 全局指针：指向实验输入的 v0(pT) 带误差的 TGraph
TFile* fMain = TFile::Open("finalExp_3040.root");
static TGraph* gExpV0_pion = (TGraph*)fMain->Get("h_v0_pion");
static TGraph* gExpV0_kaon = (TGraph*)fMain->Get("h_v0_kaon");
static TGraph* gExpV0_proton = (TGraph*)fMain->Get("h_v0_proton");

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
  gMT[i] = BW_MT(kMass[i], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[i]);

  std::cout << "betaT = " << i_betaT << ">>> Computed integrals: A=" << gIntegralA[i] << "  B=" << gIntegralB[i] << "\n";
}
double Compute_IntegralN(int i, double T) {
  ROOT::Math::Integrator intD(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  auto lamD = [&](double x) { return 2 * TMath::Pi() * x * BW_Spectrum(x, kMass[i], T, betaT, n_flow, rho2, Rx, Ry); };  // Ni
  std::function<double(double)> funD(lamD);
  intD.SetFunction(funD);
  double N = intD.Integral(kPTMin, kPTMax);
  std::cout << "T = " << T << ">>> Computed integrals: N=" << N << "\n";
  return N;
}
double N_factor[2] = {0.0, 0.0};

void run_st2_sb2() {
  TCanvas* c = new TCanvas("c1", "v0 Model", 800, 600);
  c->SetGrid();
  c->cd();
  /// 1) 计算A B
  for (int i = 0; i < kNspecies; ++i) { Compute_IntegralAB(i, 0.837); }
  ////2）画实验图
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
  gExpV0_proton->SetMarkerSize(1.5);
  gExpV0_proton->SetMarkerColor(kBlue);
  gExpV0_pion->Draw("AP");
  gExpV0_kaon->Draw("P SAME");
  gExpV0_proton->Draw("P SAME");

  TF1* f_v0_pion = new TF1(
      "f_v0_pion",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_sT2(pT, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[0], sT2, gMT[0], gIntegralN[0]);
      },
      0.,  // pT 下限
      5.,  // pT 上限
      2);
  f_v0_pion->SetParameters(3e-5, 0.0);
  f_v0_pion->SetLineColor(kGreen + 2);
  f_v0_pion->SetLineWidth(2);
  f_v0_pion->Draw("same");

  TF1* fModel_kaon = new TF1(
      "fV0Model_kaon",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_pT(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1], gIntegralA, gIntegralB, sT2, sB2);
      },
      0.5, 6., /*Npar=*/2);
  fModel_kaon->SetParameters(3e-5, 0.0);
  fModel_kaon->SetLineColor(kRed + 2);
  fModel_kaon->SetLineWidth(2);
  // fModel_kaon->Draw("same");

  TF1* f_v0_kaon = new TF1(
      "f_v0_kaon",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_sT2(pT, kMass[1], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[1], sT2, gMT[1], gIntegralN[1]);
      },
      0.,  // pT 下限
      5.,  // pT 上限
      2);
  f_v0_kaon->SetParameters(3e-5, 0.0);
  f_v0_kaon->SetLineColor(kRed);
  f_v0_kaon->SetLineWidth(2);
  f_v0_kaon->Draw("same");

  TF1* f_v0_proton = new TF1(
      "f_v0_proton",
      [&](double* x, double* p) {
        double pT = x[0];
        double sT2 = p[0], sB2 = p[1];
        return BW_V0_sT2(pT, kMass[2], Tkin, betaT, n_flow, rho2, Rx, Ry, A_const[2], sT2, gMT[2], gIntegralN[2]);
      },
      0.,  // pT 下限
      5.,  // pT 上限
      2);
  f_v0_proton->SetParameters(3e-5, 0.0);
  f_v0_proton->SetLineColor(kBlue);
  f_v0_proton->SetLineWidth(2);
  f_v0_proton->Draw("same");

  TCanvas* c2 = new TCanvas("c2", "", 800, 600);
  c2->cd();
  N_factor[0] = Compute_IntegralN(0, Tkin);
  cout << "N_factor= " << N_factor[0] << endl;
  TF1* f_test = new TF1(
      "f_test",
      // 注意：下面这个 lambda 接受两个参数 (x,p)，
      //       但是我们不用 p[...]，所以直接忽略它
      [&](double* x, double* /*p*/) {
        double pt = x[0];  // 自变量就是质量
        double N_pt = 2 * TMath::Pi() * pt * BW_Spectrum(pt, kMass[0], Tkin, betaT, n_flow, rho2, Rx, Ry);
        return N_pt / N_factor[0];
      },
      0.1,  // mass 下限
      2.0,  // mass 上限
      0     // 拟合参数个数
  );
  f_test->SetLineColor(kBlack);
  f_test->SetLineWidth(2);
  f_test->Draw();

  N_factor[1] = Compute_IntegralN(0, 0.118);
  cout << "N_factor= " << N_factor[1] << endl;
  TF1* f_test2 = new TF1(
      "f_test2",
      // 注意：下面这个 lambda 接受两个参数 (x,p)，
      //       但是我们不用 p[...]，所以直接忽略它
      [&](double* x, double* /*p*/) {
        double pt = x[0];  // 自变量就是质量
        return (1. / N_factor[1]) * 2 * TMath::Pi() * pt * BW_Spectrum(pt, kMass[0], 0.118, betaT, n_flow, rho2, Rx, Ry);
      },
      0.1,  // mass 下限
      2.0,  // mass 上限
      0     // 拟合参数个数
  );
  f_test2->SetLineColor(kBlue);
  f_test2->SetLineWidth(2);
  f_test2->Draw("same");

  ////////////////第二轮//////////////

  c->SaveAs("st2_exp.pdf");
}
