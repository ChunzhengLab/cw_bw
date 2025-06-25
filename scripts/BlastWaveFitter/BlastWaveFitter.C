// BlastWaveFit.C
// ROOT macro for blast-wave fit using Minuit2
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include "Minuit2/FCNBase.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

// ---- Output-file control ---------------------------------
static std::string gSaveName = "BlastWaveFit.pdf";
static std::string gCentLabel = "10-20%";
void SetSaveName(const char* filename) {
  gSaveName = filename;
}
void SetCentralityLabel(const char* label) {
  gCentLabel = label;
}

// ---- Fit‑summary struct & accessor -----------------------
struct BWFitResult {
  bool valid;        // fit converged?
  double chi2;       // minimum χ²
  int ndf;           // degrees of freedom
  double params[8];  // β_T, T_kin, n_flow, ρ₂, A_p, A_pion,A_kaon, R_x
};
static BWFitResult gLastFitResult;       // filled by BlastWaveFitter()
const BWFitResult& GetLastFitResult() {  // simple accessor
  return gLastFitResult;
}

using namespace ROOT::Minuit2;

// Setter so the driving macro can supply any graphs it likes
void SetDataGraphs(TGraphAsymmErrors* spec_proton, TGraphAsymmErrors* spec_pion, TGraphAsymmErrors* spec_kaon, TGraphAsymmErrors* v2_proton, TGraphAsymmErrors* v2_pion,
                   TGraphAsymmErrors* v2_kaon);
// Global pointers to data graphs (10-20% centrality)
TGraphAsymmErrors *gSpec_proton, *gSpec_pion, *gSpec_kaon, *gV2_proton, *gV2_pion, *gV2_kaon;

void SetDataGraphs(TGraphAsymmErrors* spec_proton, TGraphAsymmErrors* spec_pion, TGraphAsymmErrors* spec_kaon, TGraphAsymmErrors* v2_proton, TGraphAsymmErrors* v2_pion,
                   TGraphAsymmErrors* v2_kaon) {
  gSpec_proton = spec_proton;
  gSpec_pion = spec_pion;
  gSpec_kaon = spec_kaon;
  gV2_proton = v2_proton;
  gV2_pion = v2_pion;
  gV2_kaon = v2_kaon;
};

const double M_PROTON = 0.93827;
const double M_PION = 0.1395;
const double M_KAON = 0.493;

// Fit ranges (GeV/c)
const double PTMIN_SPEC_PROTON = 0.3;
const double PTMAX_SPEC_PROTON = 3.0;
const double PTMIN_SPEC_PION = 0.5;
const double PTMAX_SPEC_PION = 1.0;
const double PTMIN_SPEC_KAON = 0.2;
const double PTMAX_SPEC_KAON = 1.5;
const double PTMIN_V2_PROTON = 0.5;
const double PTMAX_V2_PROTON = 1.5;
const double PTMIN_V2_PION = 0.5;  //??
const double PTMAX_V2_PION = 1.0;  //??
const double PTMIN_V2_KAON = 0.5;  //??
const double PTMAX_V2_KAON = 1.0;  //??

// Declarations for integrator functions
double BW_Spectrum(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry);
double BW_v2shape(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry);

namespace {
struct BWParams {
  double pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry;
  int bessel_order;  // 0 for spectrum, 2 for v2 numerator
};

static double bw_integrand(double r, double phi_hat, void* vp) {
  auto* p = static_cast<BWParams*>(vp);
  double pT = p->pT, m = p->mass, T = p->Tkin, betaT = p->betaT, n = p->n_flow, rho2 = p->rho2, Rx = p->Rx, Ry = p->Ry;
  // Compute boost angle phi_b from phi_hat
  double phi_b = std::atan2(Rx * std::sin(phi_hat), Ry * std::cos(phi_hat));
  // radial profile in terms of transverse rapidity (exact)
  double rho0 = std::atanh(betaT);                                    // outer-edge rapidity
  double rho = std::pow(r, n) * (rho0 + rho2 * std::cos(2 * phi_b));  // Eq. (7) with optional power n
                                                                      // double betar = r * betaT;        // outer-edge rapidity
                                                                      // double rho = std::atanh(betar);  // Eq. (7) with optional power n

  // transverse mass
  double mT = std::sqrt(pT * pT + m * m);
  double argK = mT * TMath::CosH(rho) / T;
  if (argK > 700) return 0.0;
  double argI = pT * TMath::SinH(rho) / T;
  if (argI > 700) return 0.0;  // analogous to the argK guard
  if (argI < 0) argI = -argI;
  if (argI < 1e-8) argI = 1e-8;
  double K1 = ROOT::Math::cyl_bessel_k(1, argK);
  double In = ROOT::Math::cyl_bessel_i(p->bessel_order, argI);
  double weight = (p->bessel_order == 2 ? std::cos(2 * phi_b) : 1.0);  // v2 weight uses φ_b
  return r * mT * K1 * In * weight;                                    // Jacobian factor dA = R_x R_y r dr dφ;  if remove Rx*Ry 应该不会对其他参数有影响
}
}  // anonymous namespace

double BW_Spectrum(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry) {
  const int Nphi = 64;
  // For simplified spectrum: ignore anisotropy and Jacobian
  BWParams params{pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, 0};
  // BWParams params{pT, mass, Tkin, betaT, n_flow, 0.0, 10.0, 10.0, 0};
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS, 64);
  ig.SetAbsTolerance(1e-4);
  // For fixed-node Gauss-Legendre, tolerance refers to internal convergence.
  double sum = 0;
  for (int i = 0; i < Nphi; ++i) {
    double phi_hat = 2 * M_PI * i / Nphi;
    std::function<double(double)> ffun = [=, &params](double r) { return bw_integrand(r, phi_hat, &params); };
    ig.SetFunction(ffun);
    sum += ig.Integral(0, 1.0);
  }
  return sum * (2 * M_PI / Nphi);
}

double BW_v2shape(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry) {
  const int Nphi = 64;
  BWParams p0{pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, 0};
  BWParams p2{pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, 2};
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS, 64);
  ig.SetAbsTolerance(1e-6);
  double sum0 = 0, sum2 = 0;
  for (int i = 0; i < Nphi; ++i) {
    double phi_hat = 2 * M_PI * i / Nphi;
    std::function<double(double)> ffun0 = [=, &p0](double r) { return bw_integrand(r, phi_hat, &p0); };
    ig.SetFunction(ffun0);
    sum0 += ig.Integral(0, 1.0);
    std::function<double(double)> ffun2 = [=, &p2](double r) { return bw_integrand(r, phi_hat, &p2); };
    ig.SetFunction(ffun2);
    sum2 += ig.Integral(0, 1.0);
  }
  sum0 *= (2 * M_PI / Nphi);
  sum2 *= (2 * M_PI / Nphi);
  return (sum0 > 1e-15 ? sum2 / sum0 : 0);
}
/////FC for dN/dT
double BW_dSpectrum_dT(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A) {
  const double h = 1e-6;
  double N_plus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin + h, betaT, n_flow, rho2, Rx, Ry);
  double N_minus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin - h, betaT, n_flow, rho2, Rx, Ry);
  return (N_plus - N_minus) / (2 * h);
  // return (std::log(N_plus) - std::log(N_minus)) / (2 * h);
}
/////FC for dN/dBetaT
double BW_dSpectrum_dBeta(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A) {
  const double h = 1e-6;
  double N_plus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin, betaT + h, n_flow, rho2, Rx, Ry);
  double N_minus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin, betaT - h, n_flow, rho2, Rx, Ry);
  return (N_plus - N_minus) / (2 * h);
}

//////FC for dNdT
double BW_dNdT(double pT, double mass, double T, double betaT, double n_flow, double rho2, double Rx, double Ry, double A) {
  // 预计算
  const double prefac = -A * 2.0 * M_PI * pT / (T * T);
  const double mT = std::sqrt(pT * pT + mass * mass);

  // 内层：对 r 进行积分
  auto integrand_r = [&](double r, double phi_hat) {
    // 计算 phi_b
    double phi_b = std::atan2(Rx * std::sin(phi_hat), Ry * std::cos(phi_hat));
    // 计算局部 rapidity rho
    double rho0 = std::atanh(betaT);
    double rho = std::pow(r, n_flow) * (rho0 + rho2 * std::cos(2 * phi_b));

    // 计算 Bessel 参数
    double argI = pT * std::sinh(rho) / T;
    double argK = mT * std::cosh(rho) / T;

    if (argK > 700) return 0.0;
    if (argI > 700) return 0.0;  // analogous to the argK guard
    if (argI < 0) argI = -argI;
    if (argI < 1e-8) argI = 1e-8;

    // I1, I0, K1 and derivative of K1
    double I0 = ROOT::Math::cyl_bessel_i(0, argI);
    double I1 = ROOT::Math::cyl_bessel_i(1, argI);
    double K1 = ROOT::Math::cyl_bessel_k(1, argK);
    // K1' = d/dx K1(x) = -K0(x) - K1(x)/x  but here we use
    double K1p = -ROOT::Math::cyl_bessel_k(0, argK)  // K0
                 - K1 / argK;

    // the bracket [ ... ]
    double term1 = pT * std::sinh(rho) * I1 * K1;
    double term2 = mT * std::cosh(rho) * I0 * K1p;
    return r * (term1 + term2);
  };

  // Gauss–Legendre for phi and r
  const int Nphi = 128;
  ROOT::Math::Integrator ig_r(ROOT::Math::IntegrationOneDim::kGAUSS, 64);
  ig_r.SetAbsTolerance(1e-6);
  ig_r.SetRelTolerance(1e-6);

  double doubleInt = 0.0;
  // loop over phi
  for (int i = 0; i < Nphi; ++i) {
    double phi = 2 * M_PI * i / double(Nphi);
    // wrap integrand only in r
    std::function<double(double)> f_r = [&](double r) { return integrand_r(r, phi); };
    ig_r.SetFunction(f_r);
    doubleInt += ig_r.Integral(0.0, 1.0);
  }

  // 完成 phi 积分
  doubleInt *= (2 * M_PI / double(Nphi));

  // 最终结果
  return prefac * mT * doubleInt;
}

////FC for MT
double BW_MT(double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A) {
  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  auto lamA = [&](double x) { return BW_dSpectrum_dT(x, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A); };
  std::function<double(double)> funA(lamA);
  intA.SetFunction(funA);
  return intA.Integral(0.2, 7.0);
}
////FC for MB
double BW_MB(double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A) {
  ROOT::Math::Integrator intA(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-4);
  auto lamA = [&](double x) { return BW_dSpectrum_dBeta(x, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A); };
  std::function<double(double)> funA(lamA);
  intA.SetFunction(funA);
  return intA.Integral(0.2, 7.0);
}
double BW_V0_sT2(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A, double sT2, double MT, double Ntotal) {
  // const double h = 1e-4;
  // double N_plus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin + h, betaT, n_flow, rho2, Rx, Ry);
  // double N_minus = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin - h, betaT, n_flow, rho2, Rx, Ry);
  // double dIn_NdT = (std::log(N_plus) - std::log(N_minus)) / (2 * h);

  double Np = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry);
  double dNdT = BW_dSpectrum_dT(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  // double dNdT = BW_dNdT(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  return (dNdT * sqrt(sT2) - Np * MT * sqrt(sT2) / Ntotal) / Np;
}
/////   FC for v0(pt)
///   double BW_V0_pT(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A, const std::array<double, 3>& gIntegralA,
///                   const std::array<double, 3>& gIntegralB, double sT2, double sB2, double MT, double MB, double Ntotal) {
///     double Np = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry);
///     double dNdT = BW_dSpectrum_dT(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
///     double dNdB = BW_dSpectrum_dBeta(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
///     double gIntegralA_sum = 0.0;
///     double gIntegralB_sum = 0.0;
///     for (int i = 0; i < 3; ++i) {
///       gIntegralA_sum += gIntegralA[i];
///       gIntegralB_sum += gIntegralB[i];
///     }
///     double num = gIntegralA_sum * sT2 * dNdT + gIntegralB_sum * sB2 * dNdB - Np * (gIntegralA_sum * MT * sT2 + gIntegralB_sum * MB * sB2) / Ntotal;
///     double den = sqrt(pow(gIntegralA_sum, 2) * sT2 + pow(gIntegralB_sum, 2) * sB2);
///
///     return num / (den * Np);
///   }
double BW_V0_pT(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A, double gA, double gB, double sT2, double sB2,
                double MT, double MB, double Ntotal) {
  double Np = A * 2 * TMath::Pi() * pT * BW_Spectrum(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry);
  double dNdT = BW_dSpectrum_dT(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  double dNdB = BW_dSpectrum_dBeta(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  double num = gA * sT2 * dNdT + gB * sB2 * dNdB - Np * (gA * MT * sT2 + gB * MB * sB2) / Ntotal;
  double den = sqrt(pow(gA, 2) * sT2 + pow(gB, 2) * sB2);

  return num / (den * Np);
}
/////FC for d v0(pt) /dpT
double BW_dV0_dpT(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A, double gA, double gB, double sT2, double sB2,
                  double MT, double MB, double Ntotal) {
  const double h = 1e-4;
  double v0_plus = BW_V0_pT(pT + h, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A, gA, gB, sT2, sB2, MT, MB, Ntotal);
  double v0_minus = BW_V0_pT(pT - h, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A, gA, gB, sT2, sB2, MT, MB, Ntotal);
  return (v0_plus - v0_minus) / (2 * h);
}

// FCN for Minuit2
class BWFCN : public FCNBase {
 public:
  BWFCN() {
  }
  double operator()(const std::vector<double>& par) const override {
    double betaT = par[0], Tkin = par[1], n_flow = par[2];
    double rho2 = par[3];
    double A_p = par[4], A_pion = par[5], A_kaon = par[6], Rx = par[7];
    double Ry = 10.0;  // fixed

    double chi2 = 0;
    // Spectrum chi2 for proton
    for (int i = 0; i < gSpec_proton->GetN(); ++i) {
      double x, y;
      gSpec_proton->GetPoint(i, x, y);
      if (x < PTMIN_SPEC_PROTON || x > PTMAX_SPEC_PROTON) continue;
      double err = 0.5 * (gSpec_proton->GetErrorYlow(i) + gSpec_proton->GetErrorYhigh(i));
      double mshape = BW_Spectrum(x, M_PROTON, Tkin, betaT, n_flow, rho2, Rx, Ry);  // 修改 rho2
      double m = A_p * mshape;
      chi2 += pow((y - m) / err, 2);
    }
    // Spectrum chi2 for pion
    for (int i = 0; i < gSpec_pion->GetN(); ++i) {
      double x, y;
      gSpec_pion->GetPoint(i, x, y);
      if (x < PTMIN_SPEC_PION || x > PTMAX_SPEC_PION) continue;
      double err = 0.5 * (gSpec_pion->GetErrorYlow(i) + gSpec_pion->GetErrorYhigh(i));
      double mshape = BW_Spectrum(x, M_PION, Tkin, betaT, n_flow, rho2, Rx, Ry);
      double m = A_pion * mshape;
      chi2 += pow((y - m) / err, 2);
    }
    // Spectrum chi2 for kaon
    for (int i = 0; i < gSpec_kaon->GetN(); ++i) {
      double x, y;
      gSpec_kaon->GetPoint(i, x, y);
      if (x < PTMIN_SPEC_KAON || x > PTMAX_SPEC_KAON) continue;
      double err = 0.5 * (gSpec_kaon->GetErrorYlow(i) + gSpec_kaon->GetErrorYhigh(i));
      double mshape = BW_Spectrum(x, M_KAON, Tkin, betaT, n_flow, rho2, Rx, Ry);
      double m = A_kaon * mshape;
      chi2 += pow((y - m) / err, 2);
    }
    // v2 chi2 for proton
    for (int i = 0; i < gV2_proton->GetN(); ++i) {
      double x, y;
      gV2_proton->GetPoint(i, x, y);
      if (x < PTMIN_V2_PROTON || x > PTMAX_V2_PROTON) continue;
      double err = 0.5 * (gV2_proton->GetErrorYlow(i) + gV2_proton->GetErrorYhigh(i));
      double m = BW_v2shape(x, M_PROTON, Tkin, betaT, n_flow, rho2, Rx, Ry);
      chi2 += pow((y - m) / err, 2);
    }
    // v2 chi2 for pion
    for (int i = 0; i < gV2_pion->GetN(); ++i) {
      double x, y;
      gV2_pion->GetPoint(i, x, y);
      if (x < PTMIN_V2_PION || x > PTMAX_V2_PION) continue;
      double err = 0.5 * (gV2_pion->GetErrorYlow(i) + gV2_pion->GetErrorYhigh(i));
      double m = BW_v2shape(x, M_PION, Tkin, betaT, n_flow, rho2, Rx, Ry);
      chi2 += pow((y - m) / err, 2);
    }
    // v2 chi2 for kaon
    for (int i = 0; i < gV2_kaon->GetN(); ++i) {
      double x, y;
      gV2_kaon->GetPoint(i, x, y);
      if (x < PTMIN_V2_KAON || x > PTMAX_V2_KAON) continue;
      double err = 0.5 * (gV2_kaon->GetErrorYlow(i) + gV2_kaon->GetErrorYhigh(i));
      double m = BW_v2shape(x, M_KAON, Tkin, betaT, n_flow, rho2, Rx, Ry);
      chi2 += pow((y - m) / err, 2);
    }
    return chi2;
  }
  double Up() const override {
    return 1.0;
  }
};

void BlastWaveFitter() {
  // Colors
  int ci[6];
  TColor* color[6];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0], 0 / 255., 24 / 255., 113 / 255.);  // dark blue
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 65 / 255., 182 / 255., 230 / 255.);  // light blue
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 255 / 255., 88 / 255., 93 / 255.);  // red
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 255 / 255., 181 / 255., 73 / 255.);  // yellow
  ci[4] = TColor::GetFreeColorIndex();
  color[4] = new TColor(ci[4],
                        157.0 / 255.0,   // R = 157
                        78.0 / 255.0,    // G = 78
                        221.0 / 255.0);  // B = 221 :contentReference[oaicite:5]{index=5}
  ci[5] = TColor::GetFreeColorIndex();
  color[5] = new TColor(ci[5],
                        86.0 / 255.0,    // R = 86
                        78.0 / 255.0,    // G = 78
                        221.0 / 255.0);  // B = 221 :contentReference[oaicite:2]{index=2}

  // Make sure external macro already set the data graphs
  if (!gSpec_proton || !gSpec_pion || !gSpec_kaon || !gV2_proton || !gV2_pion || !gV2_kaon) {
    std::cerr << "BlastWaveFitter: ERROR – data graphs not set!\n"
              << "Call SetDataGraphs(...) before BlastWaveFitter().\n";
    return;
  }

  MnUserParameters params;
  params.Add("betaT", 0.6, 0.01, 0.2, 0.9);    // 上限设为0.999避免atanh溢出
  params.Add("Tkin", 0.13, 0.005, 0.01, 0.5);  // 温度下限
                                               // params.Add("n_flow", 0.5, 0.1, 0.0, 10.0);   // 非负
  params.Add("n_flow", 1.0, 0.1, 0.0, 10.0);
  params.Fix("n_flow");
  params.Add("rho2", 0.0, 0.01, 0.0, 0.5);
  params.Add("A_p", 1e3, 1e2);
  params.Add("A_pion", 1e3, 1e2);
  params.Add("A_kaon", 1e3, 1e2);
  params.Add("Rx", 10.0, 0.5, 0.1, 20.0);
  // Fix Ry=10 fm by adding as constant inside FCN

  BWFCN fcn;
  // Perform minimization with MnMigrad
  ROOT::Minuit2::MnMigrad migrad(fcn, params);
  auto result = migrad();
  // Detailed fit output
  std::cout << "\n===== Fit Results =====\n";
  bool fit_ok = result.IsValid();
  double chi2_min = result.Fval();
  double edm = result.Edm();
  int n_pts_total = gSpec_proton->GetN() + gSpec_pion->GetN() + gSpec_kaon->GetN() + gV2_proton->GetN() + gV2_pion->GetN() + gV2_kaon->GetN();
  int n_par_var = 8;
  int n_dof = n_pts_total - n_par_var;
  std::cout << " Fit Valid: " << (fit_ok ? "Yes" : "No") << "\n";
  std::cout << " Minimum Chi2 = " << chi2_min << "\n";
  std::cout << " Estimated Distance to Minimum (EDM) = " << edm << "\n";
  std::cout << " Number of Data Points = " << n_pts_total << "\n";
  std::cout << " Number of Free Parameters = " << n_par_var << "\n";
  if (n_dof > 0) {
    std::cout << " Number of Degrees of Freedom (NDF) = " << n_dof << "\n";
    std::cout << " Chi2 / NDF = " << chi2_min / n_dof << "\n";
  } else {
    std::cout << " NDF <= 0, Chi2/NDF not calculated.\n";
  }
  std::cout << "\n--- Fitted Parameters ---\n";
  const char* parNames[8] = {"betaT", "Tkin", "n_flow", "rho2", "A_p", "A_pion", "A_kaon", "Rx"};
  for (int i = 0; i < n_par_var; ++i) {
    double val = result.UserState().Value(i);
    double err = result.UserState().Error(i);
    std::cout << " " << std::setw(12) << parNames[i] << ": " << std::setw(10) << val << " ± " << std::setw(8) << err << "\n";
  }
  std::cout << "=======================\n\n";

  // Store the fit outcome so that external code can fetch it
  gLastFitResult.valid = fit_ok;
  gLastFitResult.chi2 = chi2_min;
  gLastFitResult.ndf = n_dof;
  for (int i = 0; i < n_par_var; ++i) gLastFitResult.params[i] = result.UserState().Value(i);

  // Create canvas with two panels
  TCanvas* c1 = new TCanvas("c1", "Joint Blast-Wave Fit Results", 1200, 500);
  c1->Divide(2, 1);

  // Extract Rx and rho2 from fit results for spectrum and v2 functions
  double Rx_fit = result.UserState().Value(7);
  double rho2_fit = result.UserState().Value(3);

  // Define fit functions for spectra and v2
  // Proton spectrum
  TF1* fSpec_p_fit = new TF1(
      "fSpec_p_fit",
      [rho2_fit, Rx_fit](double* x, double* p) {
        // p[0]=A_p, p[1]=betaT, p[2]=Tkin, p[3]=n_flow
        return p[0] * BW_Spectrum(x[0], M_PROTON,
                                  p[2],      // Tkin
                                  p[1],      // betaT
                                  p[3],      // n_flow
                                  rho2_fit,  // rho2
                                  Rx_fit,    // Rx
                                  10.0);     // Ry fixed
      },
      gSpec_proton->GetX()[0], gSpec_proton->GetX()[gSpec_proton->GetN() - 1], 4);
  fSpec_p_fit->SetParameters(result.UserState().Value(4), result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2));
  // Limit display of fit function to proton spectrum range
  // fSpec_p_fit->SetRange(PTMIN_SPEC_PROTON, PTMAX_SPEC_PROTON);
  fSpec_p_fit->SetLineColor(ci[3]);
  fSpec_p_fit->SetLineWidth(4);
  fSpec_p_fit->SetNpx(500);

  // pion spectrum
  TF1* fSpec_pion_fit = new TF1(
      "fSpec_pion_fit",
      [rho2_fit, Rx_fit](double* x, double* p) {
        // p[0]=A_L, p[1]=betaT, p[2]=Tkin, p[3]=n_flow
        return p[0] * BW_Spectrum(x[0], M_PION,
                                  p[2],      // Tkin
                                  p[1],      // betaT
                                  p[3],      // n_flow
                                  rho2_fit,  // rho2
                                  Rx_fit,    // Rx
                                  10.0);     // Ry fixed
      },
      gSpec_pion->GetX()[0], gSpec_pion->GetX()[gSpec_pion->GetN() - 1], 4);
  fSpec_pion_fit->SetParameters(result.UserState().Value(5), result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2));
  // Limit display to lambda spectrum range
  // fSpec_pion_fit->SetRange(PTMIN_SPEC_LAMBDA, PTMAX_SPEC_LAMBDA);
  fSpec_pion_fit->SetLineColor(ci[1]);
  fSpec_pion_fit->SetLineWidth(4);
  fSpec_pion_fit->SetNpx(500);

  // Kaon spectrum
  TF1* fSpec_kaon_fit = new TF1(
      "fSpec_kaon_fit",
      [rho2_fit, Rx_fit](double* x, double* p) {
        // p[0]=A_L, p[1]=betaT, p[2]=Tkin, p[3]=n_flow
        return p[0] * BW_Spectrum(x[0], M_KAON,
                                  p[2],      // Tkin
                                  p[1],      // betaT
                                  p[3],      // n_flow
                                  rho2_fit,  // rho2
                                  Rx_fit,    // Rx
                                  10.0);     // Ry fixed
      },
      gSpec_kaon->GetX()[0], gSpec_kaon->GetX()[gSpec_kaon->GetN() - 1], 4);
  fSpec_kaon_fit->SetParameters(result.UserState().Value(6), result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2));
  // Limit display to lambda spectrum range
  // fSpec_kaon_fit->SetRange(PTMIN_SPEC_LAMBDA, PTMAX_SPEC_LAMBDA);
  fSpec_kaon_fit->SetLineColor(ci[4]);
  fSpec_kaon_fit->SetLineWidth(4);
  fSpec_kaon_fit->SetNpx(500);

  // Proton v2
  TF1* fV2_p_fit = new TF1(
      "fV2_p_fit",
      [Rx_fit](double* x, double* p) {
        // p[0]=betaT, p[1]=Tkin, p[2]=n_flow, p[3]=rho2_p
        return BW_v2shape(x[0], M_PROTON,
                          p[1],    // Tkin
                          p[0],    // betaT
                          p[2],    // n_flow
                          p[3],    // rho2
                          Rx_fit,  // fitted Rx
                          10.0);   // fixed Ry
      },
      gV2_proton->GetX()[0], gV2_proton->GetX()[gV2_proton->GetN() - 1], 4);
  fV2_p_fit->SetParameters(result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2), result.UserState().Value(3));
  // Limit display to proton v2 range
  // fV2_p_fit->SetRange(PTMIN_V2_PROTON, PTMAX_V2_PROTON);
  fV2_p_fit->SetLineColor(ci[3]);
  fV2_p_fit->SetLineWidth(4);
  fV2_p_fit->SetNpx(500);

  // pion v2
  TF1* fV2_pion_fit = new TF1(
      "fV2_pion_fit",
      [Rx_fit](double* x, double* p) {
        // p[0]=betaT, p[1]=Tkin, p[2]=n_flow, p[3]=rho2_L
        return BW_v2shape(x[0], M_PION,
                          p[1],    // Tkin
                          p[0],    // betaT
                          p[2],    // n_flow
                          p[3],    // rho2
                          Rx_fit,  // fitted Rx
                          10.0);   // fixed Ry
      },
      gV2_pion->GetX()[0], gV2_pion->GetX()[gV2_pion->GetN() - 1], 4);
  fV2_pion_fit->SetParameters(result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2), result.UserState().Value(3));
  // Limit display to lambda v2 range
  // fV2_pion_fit->SetRange(PTMIN_V2_LAMBDA, PTMAX_V2_LAMBDA);
  fV2_pion_fit->SetLineColor(ci[1]);
  fV2_pion_fit->SetLineWidth(4);
  fV2_pion_fit->SetNpx(500);

  // kaon v2
  TF1* fV2_kaon_fit = new TF1(
      "fV2_kaon_fit",
      [Rx_fit](double* x, double* p) {
        // p[0]=betaT, p[1]=Tkin, p[2]=n_flow, p[3]=rho2_L
        return BW_v2shape(x[0], M_KAON,
                          p[1],    // Tkin
                          p[0],    // betaT
                          p[2],    // n_flow
                          p[3],    // rho2
                          Rx_fit,  // fitted Rx
                          10.0);   // fixed Ry
      },
      gV2_kaon->GetX()[0], gV2_kaon->GetX()[gV2_kaon->GetN() - 1], 4);
  fV2_kaon_fit->SetParameters(result.UserState().Value(0), result.UserState().Value(1), result.UserState().Value(2), result.UserState().Value(3));
  // Limit display to lambda v2 range
  // fV2_kaon_fit->SetRange(PTMIN_V2_LAMBDA, PTMAX_V2_LAMBDA);
  fV2_kaon_fit->SetLineColor(ci[4]);
  fV2_kaon_fit->SetLineWidth(4);
  fV2_kaon_fit->SetNpx(500);

  // Draw spectra in upper panel
  std::string title1 = "p, pion Spectra " + gCentLabel + ";p_{T} [GeV/c];#frac{1}{2#pi p_{T}} d^{2}N/dp_{T}dy";
  c1->cd(1)->DrawFrame(0.2, 1e-4, 10, 1e3, title1.c_str());
  gSpec_proton->SetMarkerStyle(20);
  gSpec_proton->SetMarkerColor(ci[2]);
  gSpec_proton->SetLineColor(ci[2]);
  gSpec_proton->Draw("Same P");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  gPad->SetLogy();
  gSpec_pion->SetMarkerStyle(21);
  gSpec_pion->SetMarkerColor(kBlue + 1);
  gSpec_pion->SetLineColor(kBlue + 1);
  gSpec_pion->Draw("P SAME");
  gSpec_kaon->SetMarkerStyle(22);
  gSpec_kaon->SetMarkerColor(kBlack + 1);
  gSpec_kaon->SetLineColor(kBlack + 1);
  gSpec_kaon->Draw("P SAME");
  // Set x-axis limits to fit range
  fSpec_p_fit->SetRange(PTMIN_SPEC_PROTON, PTMAX_SPEC_PROTON);
  fSpec_pion_fit->SetRange(PTMIN_SPEC_PION, PTMAX_SPEC_PION);
  fSpec_kaon_fit->SetRange(PTMIN_SPEC_KAON, PTMAX_SPEC_KAON);
  fSpec_p_fit->Draw("SAME");
  fSpec_pion_fit->Draw("SAME");
  fSpec_kaon_fit->Draw("SAME");
  TLegend* leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg1->AddEntry(gSpec_proton, "Proton Data", "p");
  leg1->AddEntry(gSpec_pion, "Pion Data", "p");
  leg1->AddEntry(gSpec_kaon, "Kaon Data", "p");
  leg1->AddEntry(fSpec_p_fit, "Proton Fit", "l");
  leg1->AddEntry(fSpec_pion_fit, "Pion Fit", "l");
  leg1->AddEntry(fSpec_kaon_fit, "Kaon Fit", "l");
  leg1->Draw();
  // Force x-axis to display 0 - 5 GeV

  // Draw v2 in lower panel
  std::string title2 = "p, pion v_{2} " + gCentLabel + ";p_{T} [GeV/c];v_{2}";
  c1->cd(2)->DrawFrame(0.2, 0, 5, 0.25, title2.c_str());
  gV2_proton->SetMarkerStyle(20);
  gV2_proton->SetMarkerColor(ci[2]);
  gV2_proton->SetLineColor(ci[2]);
  gV2_proton->Draw("Same P");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  // Set x-axis limits to v2 fit range
  // fV2_p_fit->SetRange(PTMIN_V2_PROTON, PTMAX_V2_PROTON);
  // fV2_L_fit->SetRange(PTMIN_V2_LAMBDA, PTMAX_V2_LAMBDA);
  gV2_pion->SetMarkerStyle(21);
  gV2_pion->SetMarkerColor(ci[0]);
  gV2_pion->SetLineColor(ci[0]);
  gV2_pion->Draw("P SAME");
  gV2_kaon->SetMarkerStyle(22);
  gV2_kaon->SetMarkerColor(ci[5]);
  gV2_kaon->SetLineColor(ci[5]);
  gV2_kaon->Draw("P SAME");
  fV2_p_fit->Draw("SAME");
  fV2_pion_fit->Draw("SAME");
  fV2_kaon_fit->Draw("SAME");
  fV2_p_fit->SetRange(PTMIN_V2_PROTON, PTMAX_V2_PROTON);
  fV2_pion_fit->SetRange(PTMIN_V2_PION, PTMAX_V2_PION);
  fV2_kaon_fit->SetRange(PTMIN_V2_KAON, PTMAX_V2_KAON);
  TLegend* leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg2->AddEntry(gV2_proton, "Proton Data", "p");
  leg2->AddEntry(gV2_pion, "Pion Data", "p");
  leg2->AddEntry(gV2_kaon, "Kaon Data", "p");
  leg2->AddEntry(fV2_p_fit, "Proton Fit", "l");
  leg2->AddEntry(fV2_pion_fit, "Pion Fit", "l");
  leg2->AddEntry(fV2_kaon_fit, "Kaon Fit", "l");
  leg2->Draw();
  // Force x-axis to display 0 - 5 GeV

  c1->SaveAs(gSaveName.c_str());

  TFile* fout = TFile::Open("spec_v2_number.root", "RECREATE");

  fout->cd();
  fSpec_p_fit->Write();
  fSpec_pion_fit->Write();
  fSpec_kaon_fit->Write();
  fV2_p_fit->Write();
  fV2_pion_fit->Write();
  fV2_kaon_fit->Write();
  fout->Close();
}
