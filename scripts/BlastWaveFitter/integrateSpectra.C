#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResultPtr.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

// integrateSpectra.C
// -----------------------------------------------------------------------------
// Compute the integral yield of each Λ spectrum and plot it versus centrality.
// Usage examples:
//     root -l -q integrateSpectra.C
//     root -l -q 'integrateSpectra.C("mySpectra.root")'
// -----------------------------------------------------------------------------
void integrateSpectra(const char *filename = "spec_v2.root") {
  gStyle->SetOptStat(0);

  // --------------------------------------------------------------------------
  // 1. Open the ROOT file that contains the spectra
  // --------------------------------------------------------------------------
  TFile *f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    cerr << "integrateSpectra: Cannot open file " << filename << endl;
    return;
  }

  // --------------------------------------------------------------------------
  // 2. Define graph names and corresponding centrality intervals
  // --------------------------------------------------------------------------
  const char *gname[] = {"spec_lambda_0005", "spec_lambda_0510",
                         "spec_lambda_1020", "spec_lambda_2040",
                         "spec_lambda_4060", "spec_lambda_6080",
                         "spec_lambda_8090"};
  double cLow[] = {0, 5, 10, 20, 40, 60, 80};
  double cHigh[] = {5, 10, 20, 40, 60, 80, 90};
  const int nC = 7;

  // --------------------------------------------------------------------------
  // 3. Helper: trapezoidal integral of a TGraphAsymmErrors
  // --------------------------------------------------------------------------
  auto integrateGraph = [](const TGraphAsymmErrors *g) -> double {
    double sum = 0.;
    const int n = g->GetN();
    for (int i = 1; i < n; ++i) {
      double x1, y1, x2, y2;
      g->GetPoint(i - 1, x1, y1);
      g->GetPoint(i, x2, y2);
      sum += 0.5 * (y1 + y2) * (x2 - x1);
    }
    return sum;
  };

  // --------------------------------------------------------------------------
  // 4. Loop over spectra, integrate, and store results
  // --------------------------------------------------------------------------
  double cCenter[nC];
  double yield[nC];
  // Will later need integrals of 10 - 20 % and 20 - 40 % to build 20 - 30 %
  // spectrum
  double yield1020 = -1., yield2040 = -1.;
  TGraphAsymmErrors *g2030ptr =
      nullptr; // store the interpolated 20 - 30 % spectrum for later

  cout << "-------------------------------------------------------------"
       << endl;
  cout << std::left << std::setw(20) << "Graph name" << std::setw(14)
       << "Centrality"
       << "Integral yield" << endl;
  cout << "-------------------------------------------------------------"
       << endl;

  for (int i = 0; i < nC; ++i) {
    cCenter[i] = 0.5 * (cLow[i] + cHigh[i]);

    auto g = dynamic_cast<TGraphAsymmErrors *>(f->Get(gname[i]));
    if (!g) {
      cerr << "integrateSpectra: Cannot find graph " << gname[i] << endl;
      return;
    }

    yield[i] = integrateGraph(g);
    // Cache integrals needed for 20 - 30 % interpolation
    if (i == 2)
      yield1020 = yield[i]; // 10 - 20 %
    if (i == 3)
      yield2040 = yield[i]; // 20 - 40 %
    cout << std::left << std::setw(20) << gname[i] << std::right << std::setw(6)
         << std::fixed << std::setprecision(1) << cCenter[i] << "%     "
         << std::setprecision(6) << yield[i] << endl;
  }
  cout << "-------------------------------------------------------------"
       << endl;

  // --------------------------------------------------------------------------
  // 5. Load all spectra graphs and draw on a single canvas with legend
  // --------------------------------------------------------------------------
  // canvas with two pads: left = original spectra, right = 10%-wide bins
  TCanvas *cSpec = new TCanvas("cSpec", "Lambda Spectra", 1600, 600);
  cSpec->Divide(2, 1);

  // ---------------- left pad: original spectra ----------------
  cSpec->cd(1)->DrawFrame(0, 1e-6, 10, 100);
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TGraphAsymmErrors *graphs[nC] = {nullptr};

  int colors[nC] = {kRed,    kOrange + 7, kGreen + 2, kBlue,
                    kViolet, kMagenta,    kGray + 1};

  for (int i = 0; i < nC; ++i) {
    graphs[i] = dynamic_cast<TGraphAsymmErrors *>(f->Get(gname[i]));
    if (!graphs[i])
      continue;

    graphs[i]->SetMarkerColor(colors[i]);
    graphs[i]->SetLineColor(colors[i]);
    graphs[i]->SetMarkerStyle(20);
    graphs[i]->SetMarkerSize(1.0);

    if (i == 0) {
      graphs[i]->Draw("Same P");
      graphs[i]->GetYaxis()->SetRangeUser(1e-4, 10);
      graphs[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      graphs[i]->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T} dp_{T} dy)");
    } else {
      graphs[i]->Draw("P same");
    }

    std::string label = std::to_string(static_cast<int>(cLow[i])) + "-" +
                        std::to_string(static_cast<int>(cHigh[i])) + "%";
    leg->AddEntry(graphs[i], label.c_str(), "p");
  }

  leg->Draw();

  // --------------------------------------------------------------------------
  // 6. Right pad: build 10%-wide centrality slices purely by per‑pT pol3 fits
  //    For each pT bin, perform a pol3 fit of yield vs centrality (using the
  //    seven measured bins) and evaluate the fit at the desired centrality
  //    centres:  5, 15, 25, 35, 45, 55, 65 %.
  // --------------------------------------------------------------------------
  cSpec->cd(2)->DrawFrame(0, 1e-6, 10, 100);
  gPad->SetLogy();

  const int nPred = 7;
  double targetC[nPred] = {5., 15., 25., 35., 45., 55., 65.};
  const char *targetLabel[nPred] = {"0 - 10%",  "10 - 20%", "20 - 30%",
                                    "30 - 40%", "40 - 50%", "50 - 60%",
                                    "60 - 70%"};
  const char *targetName[nPred] = {"spec_lambda_0010", "spec_lambda_1020",
                                   "spec_lambda_2030", "spec_lambda_3040",
                                   "spec_lambda_4050", "spec_lambda_5060",
                                   "spec_lambda_6070"};

  // Prepare empty graphs for each predicted slice
  std::vector<TGraphAsymmErrors *> predGraphs;
  for (int j = 0; j < nPred; ++j) {
    predGraphs.push_back(new TGraphAsymmErrors());
    predGraphs.back()->SetNameTitle(targetName[j],
                                    targetName[j]); // use desired names
  }

  // Assume all original graphs share the same pT grid
  const int nPts = graphs[0]->GetN();
  double *xArr = graphs[0]->GetX();
  double *exLowArr = graphs[0]->GetEXlow();
  double *exHighArr = graphs[0]->GetEXhigh();

  // Temporary arrays for fit
  double yArr[7];

  for (int ipt = 0; ipt < nPts; ++ipt) {
    double x_pt;
    graphs[0]->GetPoint(ipt, x_pt, yArr[0]); // retrieve x coordinate

    // Gather yields at this pT for all 7 measured centralities
    for (int k = 0; k < 7; ++k) {
      double ytmp;
      graphs[k]->GetPoint(ipt, x_pt, ytmp);
      yArr[k] = ytmp;
    }

    // Build a small TGraph for this pT bin: yield vs centrality centre
    TGraph grCen(7, cCenter, yArr);
    grCen.SetName(Form("grTmp_%d", ipt));

    TFitResultPtr fitRes = grCen.Fit("pol3", "QES");
    TF1 *fC = grCen.GetFunction("pol3");
    if (!fC) { // fallback: linear if fit failed
      fC = new TF1(Form("lin%d", ipt), "pol1", 0, 90);
      grCen.Fit(fC, "QES");
    }
    // ---- Build error vs centrality and fit with the same polynomial scheme
    // ----
    double errArr[7];
    for (int k = 0; k < 7; ++k) {
      double eLow = graphs[k]->GetErrorYlow(ipt);
      double eHigh = graphs[k]->GetErrorYhigh(ipt);
      errArr[k] = 0.5 * (eLow + eHigh); // treat as symmetric
    }

    TGraph grErr(7, cCenter, errArr);
    grErr.SetName(Form("grErr_%d", ipt));

    // Try cubic first, then quadratic, then linear if needed
    TFitResultPtr fitErr = grErr.Fit("pol3", "QES");
    TF1 *fErr = grErr.GetFunction("pol3");
    if (!fitErr || static_cast<int>(fitErr) != 0) {
      fErr = new TF1(Form("quadErr%d", ipt), "pol2", 0, 90); // quadratic
      fitErr = grErr.Fit(fErr, "QES");
    }
    if (!fitErr || static_cast<int>(fitErr) != 0) {
      fErr = new TF1(Form("linErr%d", ipt), "pol1", 0, 90); // linear
      fitErr = grErr.Fit(fErr, "QES");
    }

    // Evaluate at each desired centrality centre
    for (int j = 0; j < nPred; ++j) {
      double yPred = fC->Eval(targetC[j]);
      predGraphs[j]->SetPoint(ipt, x_pt, yPred);
      // --- use the fitted error polynomial ---
      double errY =
          std::max(0.0, fErr->Eval(targetC[j])); // ensure non‑negative
      predGraphs[j]->SetPointError(ipt, exLowArr[ipt], exHighArr[ipt], errY,
                                   errY);
      cout << errY << "    " << endl;
    }
  }

  // ---- Style & draw the predicted graphs ----
  Color_t fitColors[] = {kRed,    kOrange + 7, kGreen + 2, kBlue,
                         kViolet, kMagenta,    kGray + 1};

  TLegend *legR = new TLegend(0.56, 0.58, 0.88, 0.88);
  legR->SetBorderSize(0);
  legR->SetFillStyle(0);

  for (int j = 0; j < nPred; ++j) {
    predGraphs[j]->SetMarkerStyle(24 + j % 10);
    predGraphs[j]->SetMarkerSize(1.0);
    predGraphs[j]->SetMarkerColor(fitColors[j % 7]);
    predGraphs[j]->SetLineColor(fitColors[j % 7]);

    if (j == 0) {
      predGraphs[j]->Draw("Same P");
      predGraphs[j]->GetYaxis()->SetRangeUser(1e-4, 10);
      predGraphs[j]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      predGraphs[j]->GetYaxis()->SetTitle("d^{2}N/(2#pi p_{T} dp_{T} dy)");
    } else {
      predGraphs[j]->Draw("P same");
    }
    legR->AddEntry(predGraphs[j], targetLabel[j], "p");
  }

  legR->Draw();
  cSpec->Update();
  cSpec->SaveAs("LambdaSpectra_TwoPads.pdf");

  // --------------------------------------------------------------------------
  // 7. Save the 10 %-wide predicted spectra to a separate ROOT file
  // --------------------------------------------------------------------------
  TFile *fpred = TFile::Open("PredictedSpectra.root", "RECREATE");
  if (fpred && !fpred->IsZombie()) {
    for (int j = 0; j < nPred; ++j) {
      if (predGraphs[j])
        predGraphs[j]->Write(); // write each TGraphAsymmErrors
    }
    fpred->Close();
    cout << "Predicted spectra written to PredictedSpectra.root" << endl;
  } else {
    cerr << "Warning: could not create PredictedSpectra.root" << endl;
  }
}
