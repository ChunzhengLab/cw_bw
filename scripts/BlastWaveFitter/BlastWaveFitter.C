// BlastWaveFit.C
// ROOT macro for blast-wave fit using Minuit2
#include "Minuit2/FCNBase.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include <cmath>
#include <vector>

#include <iomanip>

#include <string>

// ---- Output-file control ---------------------------------
static std::string gSaveName = "BlastWaveFit.pdf";
static std::string gCentLabel = "10-20%";
void SetSaveName(const char *filename) { gSaveName = filename; }
void SetCentralityLabel(const char *label) { gCentLabel = label; }

// ---- Fit‑summary struct & accessor -----------------------
struct BWFitResult {
  bool valid;       // fit converged?
  double chi2;      // minimum χ²
  int ndf;          // degrees of freedom
  double params[8]; // β_T, T_kin, n_flow, ρ₂_p, ρ₂_Λ, A_p, A_Λ, R_x
};
static BWFitResult gLastFitResult;      // filled by BlastWaveFitter()
const BWFitResult &GetLastFitResult() { // simple accessor
  return gLastFitResult;
}

using namespace ROOT::Minuit2;

// Setter so the driving macro can supply any graphs it likes
void SetDataGraphs(TGraphAsymmErrors *spec_p, TGraphAsymmErrors *spec_L,
                   TGraphAsymmErrors *v2_p, TGraphAsymmErrors *v2_L);
// Global pointers to data graphs (10-20% centrality)
TGraphAsymmErrors *gSpec_p, *gSpec_L, *gV2_p, *gV2_L;

void SetDataGraphs(TGraphAsymmErrors *spec_p, TGraphAsymmErrors *spec_L,
                   TGraphAsymmErrors *v2_p, TGraphAsymmErrors *v2_L) {
  gSpec_p = spec_p;
  gSpec_L = spec_L;
  gV2_p = v2_p;
  gV2_L = v2_L;
}

const double M_PROTON = 0.93827;
const double M_LAMBDA = 1.11568;

// Fit ranges (GeV/c)
const double PTMIN_SPEC_PROTON = 0.3;
const double PTMAX_SPEC_PROTON = 3.0;
const double PTMIN_SPEC_LAMBDA = 0.5;
const double PTMAX_SPEC_LAMBDA = 2.5;
const double PTMIN_V2_PROTON = 0.5;
const double PTMAX_V2_PROTON = 1.5;
const double PTMIN_V2_LAMBDA = 0.5;
const double PTMAX_V2_LAMBDA = 1.5;

// Declarations for integrator functions
double BW_Spectrum(double pT, double mass, double Tkin, double betaT,
                   double n_flow, double rho2, double Rx, double Ry);
double BW_v2shape(double pT, double mass, double Tkin, double betaT,
                  double n_flow, double rho2, double Rx, double Ry);

namespace {
struct BWParams {
  double pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry;
  int bessel_order; // 0 for spectrum, 2 for v2 numerator
};

static double bw_integrand(double r, double phi_s, void *vp) {
  auto *p = static_cast<BWParams *>(vp);
  double pT = p->pT, m = p->mass, T = p->Tkin, betaT = p->betaT, n = p->n_flow,
         rho2 = p->rho2, Rx = p->Rx, Ry = p->Ry;
  // Compute boost angle phi_b from phi_s
  double phi_b =
      std::atan2(Rx * Rx * std::sin(phi_s), Ry * Ry * std::cos(phi_s));
  // radial profile in terms of transverse rapidity (exact)
  double rho0 = std::atanh(betaT); // outer-edge rapidity
  double rho =
      std::pow(r, n) *
      (rho0 + rho2 * std::cos(2 * phi_b)); // Eq. (7) with optional power n
  // transverse mass
  double mT = std::sqrt(pT * pT + m * m);
  double argK = mT * TMath::CosH(rho) / T;
  if (argK > 700)
    return 0.0;
  double argI = pT * TMath::SinH(rho) / T;
  if (argI > 700)
    return 0.0; // analogous to the argK guard
  if (argI < 0)
    argI = -argI;
  if (argI < 1e-8)
    argI = 1e-8;
  double K1 = ROOT::Math::cyl_bessel_k(1, argK);
  double In = ROOT::Math::cyl_bessel_i(p->bessel_order, argI);
  double weight =
      (p->bessel_order == 2 ? std::cos(2 * phi_b) : 1.0); // v2 weight uses φ_b
  return r * mT * K1 * In * weight * Rx *
         Ry; // include Jacobian factor dA = R_x R_y r dr dφ
}
} // anonymous namespace

double BW_Spectrum(double pT, double mass, double Tkin, double betaT,
                   double n_flow, double rho2, double Rx, double Ry) {
  const int Nphi = 64;
  // For simplified spectrum: ignore anisotropy and Jacobian
  BWParams params{pT, mass, Tkin, betaT, n_flow, 0.0, 10., 10., 0};
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS, 64);
  ig.SetAbsTolerance(1e-4);
  // For fixed-node Gauss-Legendre, tolerance refers to internal convergence.
  double sum = 0;
  for (int i = 0; i < Nphi; ++i) {
    double phi_s = 2 * M_PI * i / Nphi;
    std::function<double(double)> ffun = [=, &params](double r) {
      return bw_integrand(r, phi_s, &params);
    };
    ig.SetFunction(ffun);
    sum += ig.Integral(0, 1.0);
  }
  return sum * (2 * M_PI / Nphi);
}

double BW_v2shape(double pT, double mass, double Tkin, double betaT,
                  double n_flow, double rho2, double Rx, double Ry) {
  const int Nphi = 64;
  BWParams p0{pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, 0};
  BWParams p2{pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, 2};
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS, 64);
  ig.SetAbsTolerance(1e-6);
  double sum0 = 0, sum2 = 0;
  for (int i = 0; i < Nphi; ++i) {
    double phi_s = 2 * M_PI * i / Nphi;
    std::function<double(double)> ffun0 = [=, &p0](double r) {
      return bw_integrand(r, phi_s, &p0);
    };
    ig.SetFunction(ffun0);
    sum0 += ig.Integral(0, 1.0);
    std::function<double(double)> ffun2 = [=, &p2](double r) {
      return bw_integrand(r, phi_s, &p2);
    };
    ig.SetFunction(ffun2);
    sum2 += ig.Integral(0, 1.0);
  }
  sum0 *= (2 * M_PI / Nphi);
  sum2 *= (2 * M_PI / Nphi);
  return (sum0 > 1e-15 ? sum2 / sum0 : 0);
}

// FCN for Minuit2
class BWFCN : public FCNBase {
public:
  BWFCN() {}
  double operator()(const std::vector<double> &par) const override {
    double betaT = par[0], Tkin = par[1], n_flow = par[2];
    double rho2_p = par[3], rho2_L = par[4];
    double A_p = par[5], A_L = par[6], Rx = par[7];
    double Ry = 10.0; // fixed

    double chi2 = 0;
    // Spectrum chi2 for p
    for (int i = 0; i < gSpec_p->GetN(); ++i) {
      double x, y;
      gSpec_p->GetPoint(i, x, y);
      if (x < PTMIN_SPEC_PROTON || x > PTMAX_SPEC_PROTON)
        continue;
      double err = 0.5 * (gSpec_p->GetErrorYlow(i) + gSpec_p->GetErrorYhigh(i));
      double mshape =
          BW_Spectrum(x, M_PROTON, Tkin, betaT, n_flow, rho2_p, Rx, Ry);
      double m = A_p * mshape;
      chi2 += pow((y - m) / err, 2);
    }
    // Spectrum chi2 for Lambda
    for (int i = 0; i < gSpec_L->GetN(); ++i) {
      double x, y;
      gSpec_L->GetPoint(i, x, y);
      if (x < PTMIN_SPEC_LAMBDA || x > PTMAX_SPEC_LAMBDA)
        continue;
      double err = 0.5 * (gSpec_L->GetErrorYlow(i) + gSpec_L->GetErrorYhigh(i));
      double mshape =
          BW_Spectrum(x, M_LAMBDA, Tkin, betaT, n_flow, rho2_L, Rx, Ry);
      double m = A_L * mshape;
      chi2 += pow((y - m) / err, 2);
    }
    // v2 chi2 for p
    for (int i = 0; i < gV2_p->GetN(); ++i) {
      double x, y;
      gV2_p->GetPoint(i, x, y);
      if (x < PTMIN_V2_PROTON || x > PTMAX_V2_PROTON)
        continue;
      double err = 0.5 * (gV2_p->GetErrorYlow(i) + gV2_p->GetErrorYhigh(i));
      double m = BW_v2shape(x, M_PROTON, Tkin, betaT, n_flow, rho2_p, Rx, Ry);
      chi2 += pow((y - m) / err, 2);
    }
    // v2 chi2 for Lambda
    for (int i = 0; i < gV2_L->GetN(); ++i) {
      double x, y;
      gV2_L->GetPoint(i, x, y);
      if (x < PTMIN_V2_LAMBDA || x > PTMAX_V2_LAMBDA)
        continue;
      double err = 0.5 * (gV2_L->GetErrorYlow(i) + gV2_L->GetErrorYhigh(i));
      double m = BW_v2shape(x, M_LAMBDA, Tkin, betaT, n_flow, rho2_L, Rx, Ry);
      chi2 += pow((y - m) / err, 2);
    }
    return chi2;
  }
  double Up() const override { return 1.0; }
};

void BlastWaveFitter() {
  // Colors
  int ci[4];
  TColor *color[4];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0], 0 / 255., 24 / 255., 113 / 255.); // dark blue
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 65 / 255., 182 / 255., 230 / 255.); // light blue
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 255 / 255., 88 / 255., 93 / 255.); // red
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 255 / 255., 181 / 255., 73 / 255.); // yellow

  // Make sure external macro already set the data graphs
  if (!gSpec_p || !gSpec_L || !gV2_p || !gV2_L) {
    std::cerr << "BlastWaveFitter: ERROR – data graphs not set!\n"
              << "Call SetDataGraphs(...) before BlastWaveFitter().\n";
    return;
  }

  MnUserParameters params;
  params.Add("betaT", 0.6, 0.01, 0.2, 0.9);   // 上限设为0.999避免atanh溢出
  params.Add("Tkin", 0.13, 0.005, 0.01, 0.5); // 温度下限
  params.Add("n_flow", 0.5, 0.1, 0.0, 10.0);  // 非负
  params.Add("rho2_p", 0.0, 0.01, 0.0, 0.5);
  params.Add("rho2_L", 0.0, 0.01, 0.0, 0.5);
  params.Add("A_p", 1e3, 1e2);
  params.Add("A_L", 1e3, 1e2);
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
  int n_pts_total =
      gSpec_p->GetN() + gSpec_L->GetN() + gV2_p->GetN() + gV2_L->GetN();
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
  const char *parNames[8] = {"betaT",  "Tkin", "n_flow", "rho2_p",
                             "rho2_L", "A_p",  "A_L",    "Rx"};
  for (int i = 0; i < n_par_var; ++i) {
    double val = result.UserState().Value(i);
    double err = result.UserState().Error(i);
    std::cout << " " << std::setw(12) << parNames[i] << ": " << std::setw(10)
              << val << " ± " << std::setw(8) << err << "\n";
  }
  std::cout << "=======================\n\n";

  // Store the fit outcome so that external code can fetch it
  gLastFitResult.valid = fit_ok;
  gLastFitResult.chi2 = chi2_min;
  gLastFitResult.ndf = n_dof;
  for (int i = 0; i < n_par_var; ++i)
    gLastFitResult.params[i] = result.UserState().Value(i);

  // Create canvas with two panels
  TCanvas *c1 = new TCanvas("c1", "Joint Blast-Wave Fit Results", 1200, 500);
  c1->Divide(2, 1);

  // Extract Rx and rho2 from fit results for spectrum and v2 functions
  double Rx_fit = result.UserState().Value(7);
  double rho2_p_fit = result.UserState().Value(3);
  double rho2_L_fit = result.UserState().Value(4);

  // Define fit functions for spectra and v2
  // Proton spectrum
  TF1 *fSpec_p_fit = new TF1(
      "fSpec_p_fit",
      [rho2_p_fit, Rx_fit](double *x, double *p) {
        // p[0]=A_p, p[1]=betaT, p[2]=Tkin, p[3]=n_flow
        return p[0] * BW_Spectrum(x[0], M_PROTON,
                                  p[2],       // Tkin
                                  p[1],       // betaT
                                  p[3],       // n_flow
                                  rho2_p_fit, // rho2
                                  Rx_fit,     // Rx
                                  10.0);      // Ry fixed
      },
      gSpec_p->GetX()[0], gSpec_p->GetX()[gSpec_p->GetN() - 1], 4);
  fSpec_p_fit->SetParameters(
      result.UserState().Value(5), result.UserState().Value(0),
      result.UserState().Value(1), result.UserState().Value(2));
  // Limit display of fit function to proton spectrum range
  // fSpec_p_fit->SetRange(PTMIN_SPEC_PROTON, PTMAX_SPEC_PROTON);
  fSpec_p_fit->SetLineColor(ci[3]);
  fSpec_p_fit->SetLineWidth(4);
  fSpec_p_fit->SetNpx(500);

  // Lambda spectrum
  TF1 *fSpec_L_fit = new TF1(
      "fSpec_L_fit",
      [rho2_L_fit, Rx_fit](double *x, double *p) {
        // p[0]=A_L, p[1]=betaT, p[2]=Tkin, p[3]=n_flow
        return p[0] * BW_Spectrum(x[0], M_LAMBDA,
                                  p[2],       // Tkin
                                  p[1],       // betaT
                                  p[3],       // n_flow
                                  rho2_L_fit, // rho2
                                  Rx_fit,     // Rx
                                  10.0);      // Ry fixed
      },
      gSpec_L->GetX()[0], gSpec_L->GetX()[gSpec_L->GetN() - 1], 4);
  fSpec_L_fit->SetParameters(
      result.UserState().Value(6), result.UserState().Value(0),
      result.UserState().Value(1), result.UserState().Value(2));
  // Limit display to lambda spectrum range
  fSpec_L_fit->SetRange(PTMIN_SPEC_LAMBDA, PTMAX_SPEC_LAMBDA);
  fSpec_L_fit->SetLineColor(ci[1]);
  fSpec_L_fit->SetLineWidth(4);
  fSpec_L_fit->SetNpx(500);

  // Proton v2
  TF1 *fV2_p_fit = new TF1(
      "fV2_p_fit",
      [Rx_fit](double *x, double *p) {
        // p[0]=betaT, p[1]=Tkin, p[2]=n_flow, p[3]=rho2_p
        return BW_v2shape(x[0], M_PROTON,
                          p[1],   // Tkin
                          p[0],   // betaT
                          p[2],   // n_flow
                          p[3],   // rho2_p
                          Rx_fit, // fitted Rx
                          10.0);  // fixed Ry
      },
      gV2_p->GetX()[0], gV2_p->GetX()[gV2_p->GetN() - 1], 4);
  fV2_p_fit->SetParameters(
      result.UserState().Value(0), result.UserState().Value(1),
      result.UserState().Value(2), result.UserState().Value(3));
  // Limit display to proton v2 range
  // fV2_p_fit->SetRange(PTMIN_V2_PROTON, PTMAX_V2_PROTON);
  fV2_p_fit->SetLineColor(ci[3]);
  fV2_p_fit->SetLineWidth(4);
  fV2_p_fit->SetNpx(500);

  // Lambda v2
  TF1 *fV2_L_fit = new TF1(
      "fV2_L_fit",
      [Rx_fit](double *x, double *p) {
        // p[0]=betaT, p[1]=Tkin, p[2]=n_flow, p[3]=rho2_L
        return BW_v2shape(x[0], M_LAMBDA,
                          p[1],   // Tkin
                          p[0],   // betaT
                          p[2],   // n_flow
                          p[3],   // rho2_L
                          Rx_fit, // fitted Rx
                          10.0);  // fixed Ry
      },
      gV2_L->GetX()[0], gV2_L->GetX()[gV2_L->GetN() - 1], 4);
  fV2_L_fit->SetParameters(
      result.UserState().Value(0), result.UserState().Value(1),
      result.UserState().Value(2), result.UserState().Value(4));
  // Limit display to lambda v2 range
  // fV2_L_fit->SetRange(PTMIN_V2_LAMBDA, PTMAX_V2_LAMBDA);
  fV2_L_fit->SetLineColor(ci[1]);
  fV2_L_fit->SetLineWidth(4);
  fV2_L_fit->SetNpx(500);

  // Draw spectra in upper panel
  std::string title1 = "p, #Lambda Spectra " + gCentLabel +
                       ";p_{T} [GeV/c];#frac{1}{2#pi p_{T}} d^{2}N/dp_{T}dy";
  c1->cd(1)->DrawFrame(0.2, 1e-4, 10, 1e2, title1.c_str());
  gSpec_p->SetMarkerStyle(20);
  gSpec_p->SetMarkerColor(ci[2]);
  gSpec_p->SetLineColor(ci[2]);
  gSpec_p->Draw("Same P");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  gPad->SetLogy();
  gSpec_L->SetMarkerStyle(21);
  gSpec_L->SetMarkerColor(kBlue + 1);
  gSpec_L->SetLineColor(kBlue + 1);
  gSpec_L->Draw("P SAME");
  // Set x-axis limits to fit range
  fSpec_p_fit->SetRange(PTMIN_SPEC_PROTON, PTMAX_SPEC_PROTON);
  fSpec_L_fit->SetRange(PTMIN_SPEC_LAMBDA, PTMAX_SPEC_LAMBDA);
  fSpec_p_fit->Draw("SAME");
  fSpec_L_fit->Draw("SAME");
  TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg1->AddEntry(gSpec_p, "Proton Data", "p");
  leg1->AddEntry(gSpec_L, "Lambda Data", "p");
  leg1->AddEntry(fSpec_p_fit, "Proton Fit", "l");
  leg1->AddEntry(fSpec_L_fit, "Lambda Fit", "l");
  leg1->Draw();
  // Force x-axis to display 0 - 5 GeV

  // Draw v2 in lower panel
  std::string title2 =
      "p, #Lambda v_{2} " + gCentLabel + ";p_{T} [GeV/c];v_{2}";
  c1->cd(2)->DrawFrame(0.2, 0, 5, 0.25, title2.c_str());
  gV2_p->SetMarkerStyle(20);
  gV2_p->SetMarkerColor(ci[2]);
  gV2_p->SetLineColor(ci[2]);
  gV2_p->Draw("Same P");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  // Set x-axis limits to v2 fit range
  // fV2_p_fit->SetRange(PTMIN_V2_PROTON, PTMAX_V2_PROTON);
  // fV2_L_fit->SetRange(PTMIN_V2_LAMBDA, PTMAX_V2_LAMBDA);
  gV2_L->SetMarkerStyle(21);
  gV2_L->SetMarkerColor(ci[0]);
  gV2_L->SetLineColor(ci[0]);
  gV2_L->Draw("P SAME");
  fV2_p_fit->Draw("SAME");
  fV2_L_fit->Draw("SAME");
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg2->AddEntry(gV2_p, "Proton Data", "p");
  leg2->AddEntry(gV2_L, "Lambda Data", "p");
  leg2->AddEntry(fV2_p_fit, "Proton Fit", "l");
  leg2->AddEntry(fV2_L_fit, "Lambda Fit", "l");
  leg2->Draw();
  // Force x-axis to display 0 - 5 GeV

  c1->SaveAs(gSaveName.c_str());
}
