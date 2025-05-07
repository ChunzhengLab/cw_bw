#include "AnalyzerCVE.h"

#include <TFile.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

AnalyzerCVE::AnalyzerCVE() {
  // Initialize all pointers to nullptr
  h_pt_lambda = nullptr;
  h_pt_proton = nullptr;
  h_eta_lambda = nullptr;
  h_eta_proton = nullptr;
  h_phi_lambda = nullptr;
  h_phi_proton = nullptr;

  p_v2_pt_lambda = nullptr;
  p_v2_pt_proton = nullptr;

  p_v2_lambda_sumPt = nullptr;
  p_v2_proton_sumPt = nullptr;
  p_v2_lambda_etaGap = nullptr;
  p_v2_proton_etaGap = nullptr;
}

AnalyzerCVE::~AnalyzerCVE() {
  // Delete single-particle histograms
  delete h_pt_lambda;
  delete h_pt_proton;
  delete h_eta_lambda;
  delete h_eta_proton;
  delete h_phi_lambda;
  delete h_phi_proton;

  // Delete v2 profiles
  delete p_v2_pt_lambda;
  delete p_v2_pt_proton;

  // Delete v2 vs sumPt and etaGap
  delete p_v2_lambda_sumPt;
  delete p_v2_proton_sumPt;
  delete p_v2_lambda_etaGap;
  delete p_v2_proton_etaGap;

  // Delete coarse centrality maps
  for (auto& kv : p_gamma) delete kv.second;
  for (auto& kv : p_delta) delete kv.second;

  for (auto& kv : p_gamma_LBCFriend) delete kv.second;
  for (auto& kv : p_delta_LBCFriend) delete kv.second;

  // Delete differential sumPt maps
  for (int b = 0; b < 3; ++b) {
    for (auto& kv : p_delta_sumPt[b]) delete kv.second;
    for (auto& kv : p_gamma_sumPt[b]) delete kv.second;
  }
  // Delete differential etaGap maps
  for (int b = 0; b < 4; ++b) {
    for (auto& kv : p_delta_etaGap[b]) delete kv.second;
    for (auto& kv : p_gamma_etaGap[b]) delete kv.second;
  }
  // Delete this-cent maps
  for (auto& kv : p_delta_vs_sumPt) delete kv.second;
  for (auto& kv : p_gamma_vs_sumPt) delete kv.second;
  for (auto& kv : p_delta_vs_etaGap) delete kv.second;
  for (auto& kv : p_gamma_vs_etaGap) delete kv.second;
}

void AnalyzerCVE::Init() {
  // Single-particle histograms
  h_pt_lambda = new TH1D("h_pt_lambda", "p_{T} #Lambda; p_{T} (GeV/c); Counts", 100, 0, 10);
  h_pt_proton = new TH1D("h_pt_proton", "p_{T} proton; p_{T} (GeV/c); Counts", 100, 0, 10);
  h_eta_lambda = new TH1D("h_eta_lambda", "#eta #Lambda; #eta; Counts", 50, -1.0, 1.0);
  h_eta_proton = new TH1D("h_eta_proton", "#eta proton; #eta; Counts", 50, -1.0, 1.0);
  h_phi_lambda = new TH1D("h_phi_lambda", "#phi #Lambda; #phi; Counts", 64, 0, 2 * M_PI);
  h_phi_proton = new TH1D("h_phi_proton", "#phi proton; #phi; Counts", 64, 0, 2 * M_PI);

  // v2 TProfiles vs pT
  p_v2_pt_lambda = new TProfile("p_v2_pt_lambda", "v_{2} vs p_{T} #Lambda; p_{T} (GeV/c); v_{2}", 100, 0, 10);
  p_v2_pt_proton = new TProfile("p_v2_pt_proton", "v_{2} vs p_{T} proton; p_{T} (GeV/c); v_{2}", 100, 0, 10);

  // v2 vs sumPt and etaGap
  p_v2_lambda_sumPt = new TProfile("p_v2_lambda_sumPt", "v_{2} #Lambda vs sum p_{T}; sum p_{T}; v_{2}", 100, 0, 10);
  p_v2_proton_sumPt = new TProfile("p_v2_proton_sumPt", "v_{2} proton vs sum p_{T}; sum p_{T}; v_{2}", 100, 0, 10);
  p_v2_lambda_etaGap = new TProfile("p_v2_lambda_etaGap", "v_{2} #Lambda vs #Delta#eta; #Delta#eta; v_{2}", 50, 0, 2);
  p_v2_proton_etaGap = new TProfile("p_v2_proton_etaGap", "v_{2} proton vs #Delta#eta; #Delta#eta; v_{2}", 50, 0, 2);

  // Initialize coarse species-pair profiles using pid-based map
  const std::vector<int> pids = {2212, 3122, -2212, -3122};
  // Helper to map PID to species name
  auto pidName = [](int pid) -> const char* {
    switch (pid) {
      case 2212:
        return "proton";
      case -2212:
        return "antiproton";
      case 3122:
        return "lambda";
      case -3122:
        return "antilambda";
      default:
        return Form("%d", pid);
    }
  };
  for (size_t ii = 0; ii < pids.size(); ++ii) {
    for (size_t jj = ii; jj < pids.size(); ++jj) {
      int pid1Sorted = pids[ii];
      int pid2Sorted = pids[jj];
      if (pid1Sorted > pid2Sorted) std::swap(pid1Sorted, pid2Sorted);
      auto key = std::make_pair(pid1Sorted, pid2Sorted);
      // gamma profile vs centrality
      {
        TString name = Form("p_gamma_%s_%s", pidName(pid1Sorted), pidName(pid2Sorted));
        p_gamma[key] = new TProfile(name, name, 10, 0, 100);
        p_gamma_LBCFriend[key] = new TProfile(name + "_LBCFriend", name+"_LBCFriend", 10, 0, 100);
      }
      // delta profile vs centrality
      {
        TString name = Form("p_delta_%s_%s", pidName(pid1Sorted), pidName(pid2Sorted));
        p_delta[key] = new TProfile(name, name, 10, 0, 100);
        p_delta_LBCFriend[key] = new TProfile(name + "_LBCFriend", name+"_LBCFriend", 10, 0, 100);
      }
    }
  }

  // Initialize differential sumPt maps
  for (size_t ii = 0; ii < pids.size(); ++ii) {
    for (size_t jj = ii; jj < pids.size(); ++jj) {
      int pid1Sorted = pids[ii];
      int pid2Sorted = pids[jj];
      if (pid1Sorted > pid2Sorted) std::swap(pid1Sorted, pid2Sorted);
      auto key = std::make_pair(pid1Sorted, pid2Sorted);
      // sumPt bins (3)
      for (int b = 0; b < 3; ++b) {
        TString name = Form("p_delta_%s_%s_sumPt_%d", pidName(pid1Sorted), pidName(pid2Sorted), b);
        p_delta_sumPt[b][key] = new TProfile(name, name, 10, 0, 100);
        name = Form("p_gamma_%s_%s_sumPt_%d", pidName(pid1Sorted), pidName(pid2Sorted), b);
        p_gamma_sumPt[b][key] = new TProfile(name, name, 10, 0, 100);
      }
      // etaGap bins (4)
      for (int b = 0; b < 4; ++b) {
        TString name = Form("p_delta_%s_%s_etaGap_%d", pidName(pid1Sorted), pidName(pid2Sorted), b);
        p_delta_etaGap[b][key] = new TProfile(name, name, 10, 0, 100);
        name = Form("p_gamma_%s_%s_etaGap_%d", pidName(pid1Sorted), pidName(pid2Sorted), b);
        p_gamma_etaGap[b][key] = new TProfile(name, name, 10, 0, 100);
      }
      // this-cent profiles
      {
        TString name = Form("p_delta_%s_%s_vs_sumPt", pidName(pid1Sorted), pidName(pid2Sorted));
        p_delta_vs_sumPt[key] = new TProfile(name, name, 50, 0, 20);
      }
      {
        TString name = Form("p_gamma_%s_%s_vs_sumPt", pidName(pid1Sorted), pidName(pid2Sorted));
        p_gamma_vs_sumPt[key] = new TProfile(name, name, 50, 0, 20);
      }
      {
        TString name = Form("p_delta_%s_%s_vs_etaGap", pidName(pid1Sorted), pidName(pid2Sorted));
        p_delta_vs_etaGap[key] = new TProfile(name, name, 50, 0, 1.6);
      }
      {
        TString name = Form("p_gamma_%s_%s_vs_etaGap", pidName(pid1Sorted), pidName(pid2Sorted));
        p_gamma_vs_etaGap[key] = new TProfile(name, name, 50, 0, 1.6);
      }
    }
  }
}

void AnalyzerCVE::Process(const Event& evt) {
  // Lazy init
  if (!h_pt_lambda) Init();

  const auto& particles = evt.GetParticles();
  float centrality = evt.Centrality();

  // Single-particle spectra and flow vs pT
  for (const auto& p : particles) {
    int pid = p.PID();
    if (pid != 2212 && pid != 3122 && pid != -2212 && pid != -3122) continue;
    float pt = p.Pt();
    float eta = p.Eta();
    float phi = p.Phi();
    if (phi < 0) phi += 2 * M_PI;
    if (pt < 0.2f || pt > 5.0f || std::fabs(eta) > 0.8f) continue;

    if (pid == 3122) {
      h_pt_lambda->Fill(pt);
      h_eta_lambda->Fill(eta);
      h_phi_lambda->Fill(phi);
      p_v2_pt_lambda->Fill(pt, std::cos(2 * phi));
    } else {
      h_pt_proton->Fill(pt);
      h_eta_proton->Fill(eta);
      h_phi_proton->Fill(phi);
      p_v2_pt_proton->Fill(pt, std::cos(2 * phi));
    }
  }

  // Two-particle correlations
  size_t n = particles.size();
  for (size_t i = 0; i < n; ++i) {
    const auto& p1 = particles[i];
    int pid1 = p1.PID();
    float pt1 = p1.Pt(), eta1 = p1.Eta(), phi1 = p1.Phi();
    if ((pid1 != 2212 && pid1 != 3122 && pid1 != -2212 && pid1 != -3122) || pt1 < 0.2f || pt1 > 5.0f ||
        std::fabs(eta1) > 0.8f)
      continue;

    for (size_t j = i + 1; j < n; ++j) {
      const auto& p2 = particles[j];
      int pid2 = p2.PID();
      float pt2 = p2.Pt(), eta2 = p2.Eta(), phi2 = p2.Phi();
      if ((pid2 != 2212 && pid2 != 3122 && pid2 != -2212 && pid2 != -3122) || pt2 < 0.2f || pt2 > 5.0f ||
          std::fabs(eta2) > 0.8f)
        continue;

      // compute obvs.
      float gamma = std::cos(phi1 + phi2);
      float delta = std::cos(phi1 - phi2);
      // compute differences
      float eta_gap = std::abs(eta1 - eta2);
      float sumPt = pt1 + pt2;

      // species-pair fill for all maps
      int pid1Sorted = pid1;
      int pid2Sorted = pid2;
      if (pid1Sorted > pid2Sorted) std::swap(pid1Sorted, pid2Sorted);
      auto key = std::make_pair(pid1Sorted, pid2Sorted);
      // coarse centrality
      p_gamma.at(key)->Fill(centrality, gamma);
      p_delta.at(key)->Fill(centrality, delta);
      if (p1.GetSerialNumberLBCFriend() == p2.GetSerialNumber() || p2.GetSerialNumberLBCFriend() == p1.GetSerialNumber()) {
          p_gamma_LBCFriend.at(key)->Fill(centrality,gamma);
          p_delta_LBCFriend.at(key)->Fill(centrality,delta);
      }
      // compute bins
      int etaGapBin{-1};
      if (eta_gap > 0.8f) etaGapBin = 3;
      else if (eta_gap > 0.6f) etaGapBin = 2;
      else if (eta_gap > 0.4f) etaGapBin = 1;
      else if (eta_gap > 0.2f) etaGapBin = 0;

      int sumPtBin{-1};
      if (sumPt > 5.0f && sumPt <= 8.0f) sumPtBin = 2;
      else if (sumPt > 3.0f && sumPt <= 5.0f) sumPtBin = 1;
      else if (sumPt > 1.0f && sumPt <= 3.0f) sumPtBin = 0;

      // differential fills
      if (sumPtBin > -1) {
        p_delta_sumPt[sumPtBin].at(key)->Fill(centrality, delta);
        p_gamma_sumPt[sumPtBin].at(key)->Fill(centrality, gamma);
      }
      if (etaGapBin > -1) {
        p_delta_etaGap[etaGapBin].at(key)->Fill(centrality, delta);
        p_gamma_etaGap[etaGapBin].at(key)->Fill(centrality, gamma);
      }
      // this-cent fills
      p_delta_vs_sumPt.at(key)->Fill(sumPt, delta);
      p_gamma_vs_sumPt.at(key)->Fill(sumPt, gamma);
      p_delta_vs_etaGap.at(key)->Fill(eta_gap, delta);
      p_gamma_vs_etaGap.at(key)->Fill(eta_gap, gamma);

      // Flow vs sumPt/etaGap
      if (pid1 == 3122) {
        p_v2_lambda_sumPt->Fill(sumPt, std::cos(2 * phi1));
        p_v2_lambda_etaGap->Fill(eta_gap, std::cos(2 * phi1));
      }
      if (pid1 == 2212) {
        p_v2_proton_sumPt->Fill(sumPt, std::cos(2 * phi1));
        p_v2_proton_etaGap->Fill(eta_gap, std::cos(2 * phi1));
      }
    }
  }
}

void AnalyzerCVE::Write(const std::string& filename) const {
  TFile f(filename.c_str(), "RECREATE");
  // Single-particle spectra
  h_pt_lambda->Write();
  h_pt_proton->Write();
  h_eta_lambda->Write();
  h_eta_proton->Write();
  h_phi_lambda->Write();
  h_phi_proton->Write();
  // Flow vs pT
  p_v2_pt_lambda->Write();
  p_v2_pt_proton->Write();
  // Flow vs sumPt and etaGap
  p_v2_lambda_sumPt->Write();
  p_v2_proton_sumPt->Write();
  p_v2_lambda_etaGap->Write();
  p_v2_proton_etaGap->Write();

  // Write coarse centrality profiles
  for (auto& kv : p_gamma) kv.second->Write();
  for (auto& kv : p_delta) kv.second->Write();

  for (auto& kv : p_gamma_LBCFriend) kv.second->Write();
  for (auto& kv : p_delta_LBCFriend) kv.second->Write();

  // Write differential sumPt profiles
  for (int b = 0; b < 3; ++b)
    for (auto& kv : p_delta_sumPt[b]) kv.second->Write();
  for (int b = 0; b < 3; ++b)
    for (auto& kv : p_gamma_sumPt[b]) kv.second->Write();
  // Write differential etaGap profiles
  for (int b = 0; b < 4; ++b)
    for (auto& kv : p_delta_etaGap[b]) kv.second->Write();
  for (int b = 0; b < 4; ++b)
    for (auto& kv : p_gamma_etaGap[b]) kv.second->Write();
  // Write this-cent profiles
  for (auto& kv : p_delta_vs_sumPt) kv.second->Write();
  for (auto& kv : p_gamma_vs_sumPt) kv.second->Write();
  for (auto& kv : p_delta_vs_etaGap) kv.second->Write();
  for (auto& kv : p_gamma_vs_etaGap) kv.second->Write();

  f.Close();
}
