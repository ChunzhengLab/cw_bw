#include "AnalyzerCVE.h"
#include <TFile.h>
#include <TString.h>
#include <algorithm>
#include <cmath>

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

  h_pair_tot_etaGap = nullptr;
  h_pair_tot_sumPt = nullptr;
  h_pair_lbc_etaGap = nullptr;
  h_pair_lbc_sumPt = nullptr;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      p_gamma_pair[i][j] = nullptr;
      p_delta_pair[i][j] = nullptr;
    }
  }

  // Initialize new differential arrays to nullptr
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        p_delta_pair_etaGap[i][j][k] = nullptr;
        p_gamma_pair_etaGap[i][j][k] = nullptr;
      }
      for (int k = 0; k < 3; ++k) {
        p_delta_pair_sumPt[i][j][k] = nullptr;
        p_gamma_pair_sumPt[i][j][k] = nullptr;
      }
      p_delta_pair_sumPt_thisCent[i][j] = nullptr;
      p_gamma_pair_sumPt_thisCent[i][j] = nullptr;
      p_delta_pair_etaGap_thisCent[i][j] = nullptr;
      p_gamma_pair_etaGap_thisCent[i][j] = nullptr;
    }
  }
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

  // Delete pair count histograms
  delete h_pair_tot_etaGap;
  delete h_pair_tot_sumPt;
  delete h_pair_lbc_etaGap;
  delete h_pair_lbc_sumPt;

  // Delete v2 vs sumPt and etaGap
  delete p_v2_lambda_sumPt;
  delete p_v2_proton_sumPt;
  delete p_v2_lambda_etaGap;
  delete p_v2_proton_etaGap;

  // Delete new differential arrays and this-cent arrays, and coarse matrices
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        delete p_delta_pair_etaGap[i][j][k];
        delete p_gamma_pair_etaGap[i][j][k];
      }
      for (int k = 0; k < 3; ++k) {
        delete p_delta_pair_sumPt[i][j][k];
        delete p_gamma_pair_sumPt[i][j][k];
      }
      delete p_delta_pair_sumPt_thisCent[i][j];
      delete p_gamma_pair_sumPt_thisCent[i][j];
      delete p_delta_pair_etaGap_thisCent[i][j];
      delete p_gamma_pair_etaGap_thisCent[i][j];
      // delete coarse correlation profiles
      delete p_gamma_pair[i][j];
      delete p_delta_pair[i][j];
    }
  }
}

void AnalyzerCVE::Init() {
  // Single-particle histograms
  h_pt_lambda = new TH1D("h_pt_lambda", "p_{T} #Lambda; p_{T} (GeV/c); Counts",
                         100, 0, 10);
  h_pt_proton = new TH1D("h_pt_proton", "p_{T} proton; p_{T} (GeV/c); Counts",
                         100, 0, 10);
  h_eta_lambda =
      new TH1D("h_eta_lambda", "#eta #Lambda; #eta; Counts", 50, -1.0, 1.0);
  h_eta_proton =
      new TH1D("h_eta_proton", "#eta proton; #eta; Counts", 50, -1.0, 1.0);
  h_phi_lambda =
      new TH1D("h_phi_lambda", "#phi #Lambda; #phi; Counts", 64, 0, 2 * M_PI);
  h_phi_proton =
      new TH1D("h_phi_proton", "#phi proton; #phi; Counts", 64, 0, 2 * M_PI);

  // v2 TProfiles vs pT
  p_v2_pt_lambda =
      new TProfile("p_v2_pt_lambda",
                   "v_{2} vs p_{T} #Lambda; p_{T} (GeV/c); v_{2}", 100, 0, 10);
  p_v2_pt_proton =
      new TProfile("p_v2_pt_proton",
                   "v_{2} vs p_{T} proton; p_{T} (GeV/c); v_{2}", 100, 0, 10);

  // v2 vs sumPt and etaGap
  p_v2_lambda_sumPt =
      new TProfile("p_v2_lambda_sumPt",
                   "v_{2} #Lambda vs sum p_{T}; sum p_{T}; v_{2}", 100, 0, 10);
  p_v2_proton_sumPt =
      new TProfile("p_v2_proton_sumPt",
                   "v_{2} proton vs sum p_{T}; sum p_{T}; v_{2}", 100, 0, 10);
  p_v2_lambda_etaGap =
      new TProfile("p_v2_lambda_etaGap",
                   "v_{2} #Lambda vs #Delta#eta; #Delta#eta; v_{2}", 50, 0, 2);
  p_v2_proton_etaGap =
      new TProfile("p_v2_proton_etaGap",
                   "v_{2} proton vs #Delta#eta; #Delta#eta; v_{2}", 50, 0, 2);

  // Particle species names for pairing
  static const char *spNames[4] = {"p", "L", "pBar", "LBar"};
  // Species-pair differentials: etaGap (4 bins) and sumPt (3 bins)
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      std::string tag = std::string(spNames[i]) + spNames[j];
      // 4 etaGap bins
      for (int k = 0; k < 4; ++k) {
        p_delta_pair_etaGap[i][j][k] = new TProfile(
            Form("p_delta_%s_etaGap_%d", tag.c_str(), k),
            Form("#delta %s vs #Delta#eta bin %d; #Delta#eta; #delta",
                 tag.c_str(), k),
            50, 0, 2);
        p_gamma_pair_etaGap[i][j][k] = new TProfile(
            Form("p_gamma_%s_etaGap_%d", tag.c_str(), k),
            Form("#gamma %s vs #Delta#eta bin %d; #Delta#eta; #gamma",
                 tag.c_str(), k),
            50, 0, 2);
      }
      // 3 sumPt bins
      for (int k = 0; k < 3; ++k) {
        p_delta_pair_sumPt[i][j][k] = new TProfile(
            Form("p_delta_%s_sumPt_%d", tag.c_str(), k),
            Form("#delta %s vs sum p_{T} bin %d; sum p_{T}; #delta",
                 tag.c_str(), k),
            50, 0, 20);
        p_gamma_pair_sumPt[i][j][k] = new TProfile(
            Form("p_gamma_%s_sumPt_%d", tag.c_str(), k),
            Form("#gamma %s vs sum p_{T} bin %d; sum p_{T}; #gamma",
                 tag.c_str(), k),
            50, 0, 20);
      }
      // this-cent differentials
      p_delta_pair_sumPt_thisCent[i][j] = new TProfile(
          Form("p_delta_%s_sumPt_thisCent", tag.c_str()),
          Form("#delta %s vs sum p_{T} (this cent); sum p_{T}; #delta",
               tag.c_str()),
          50, 0, 20);
      p_gamma_pair_sumPt_thisCent[i][j] = new TProfile(
          Form("p_gamma_%s_sumPt_thisCent", tag.c_str()),
          Form("#gamma %s vs sum p_{T} (this cent); sum p_{T}; #gamma",
               tag.c_str()),
          50, 0, 20);
      p_delta_pair_etaGap_thisCent[i][j] = new TProfile(
          Form("p_delta_%s_etaGap_thisCent", tag.c_str()),
          Form("#delta %s vs #Delta#eta (this cent); #Delta#eta; #delta",
               tag.c_str()),
          50, 0, 2);
      p_gamma_pair_etaGap_thisCent[i][j] = new TProfile(
          Form("p_gamma_%s_etaGap_thisCent", tag.c_str()),
          Form("#gamma %s vs #Delta#eta (this cent); #Delta#eta; #gamma",
               tag.c_str()),
          50, 0, 2);
    }
  }

  // Pair count histograms
  h_pair_tot_etaGap =
      new TH1D("h_pair_tot_etaGap",
               "Total pairs vs #Delta#eta; #Delta#eta; Counts", 50, 0, 2);
  h_pair_tot_sumPt =
      new TH1D("h_pair_tot_sumPt",
               "Total pairs vs sum p_{T}; sum p_{T}; Counts", 50, 0, 20);
  h_pair_lbc_etaGap =
      new TH1D("h_pair_lbc_etaGap",
               "LBC pairs vs #Delta#eta; #Delta#eta; Counts", 50, 0, 2);
  h_pair_lbc_sumPt =
      new TH1D("h_pair_lbc_sumPt", "LBC pairs vs sum p_{T}; sum p_{T}; Counts",
               50, 0, 20);

  // Species-pair correlation profiles
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      std::string tag = std::string(spNames[i]) + spNames[j];
      p_gamma_pair[i][j] = new TProfile(
          Form("p_gamma_%s", tag.c_str()),
          Form("#gamma %s; cent; #gamma", tag.c_str()), 20, 0, 100);
      p_delta_pair[i][j] = new TProfile(
          Form("p_delta_%s", tag.c_str()),
          Form("#delta %s; cent; #delta", tag.c_str()), 20, 0, 100);
    }
  }
}

void AnalyzerCVE::Process(const Event &evt) {
  // Lazy init
  if (!h_pt_lambda)
    Init();

  auto speciesIdx = [](int pid) {
    switch (pid) {
    case 2212:
      return 0;
    case 3122:
      return 1;
    case -2212:
      return 2;
    case -3122:
      return 3;
    default:
      return -1;
    }
  };

  const auto &particles = evt.GetParticles();
  float centrality = evt.Centrality();

  // Single-particle spectra and flow vs pT
  for (const auto &p : particles) {
    int pid = p.PID();
    if (pid != 2212 && pid != 3122)
      continue;
    float pt = p.Pt();
    float eta = p.Eta();
    float phi = TVector2::Phi_0_2pi(p.Phi());
    if (pt < 0.2f || pt > 5.0f || std::fabs(eta) > 0.8f)
      continue;

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
    const auto &p1 = particles[i];
    int pid1 = p1.PID();
    float pt1 = p1.Pt(), eta1 = p1.Eta(), phi1 = p1.Phi();
    if ((pid1 != 2212 && pid1 != 3122) || pt1 < 0.2f || pt1 > 5.0f ||
        std::fabs(eta1) > 0.8f)
      continue;

    for (size_t j = i + 1; j < n; ++j) {
      const auto &p2 = particles[j];
      int pid2 = p2.PID();
      float pt2 = p2.Pt(), eta2 = p2.Eta(), phi2 = p2.Phi();
      if ((pid2 != 2212 && pid2 != 3122) || pt2 < 0.2f || pt2 > 5.0f ||
          std::fabs(eta2) > 0.8f)
        continue;

      float gamma = std::cos(phi1 + phi2);
      float delta = std::cos(phi1 - phi2);
      // species-pair
      int idx1 = speciesIdx(pid1), idx2 = speciesIdx(pid2);
      if (idx1 >= 0 && idx2 >= 0) {
        int a = std::min(idx1, idx2), b = std::max(idx1, idx2);
        p_gamma_pair[a][b]->Fill(centrality, gamma);
        p_delta_pair[a][b]->Fill(centrality, delta);
      }

      float eta_gap = std::abs(eta1 - eta2);
      float sumPt = pt1 + pt2;

      // Count histograms
      h_pair_tot_etaGap->Fill(eta_gap);
      h_pair_tot_sumPt->Fill(sumPt);

      // LBC pairs
      if (p1.GetSerialNumberLBCFriend() >= 0 &&
          p2.GetSerialNumberLBCFriend() >= 0) {
        h_pair_lbc_etaGap->Fill(eta_gap);
        h_pair_lbc_sumPt->Fill(sumPt);
      }

      // Binning
      int etaBin = (eta_gap < 0.5f)   ? 0
                   : (eta_gap < 1.0f) ? 1
                   : (eta_gap < 1.5f) ? 2
                                      : 3;
      int sumBin = (sumPt < 5.0f) ? 0 : (sumPt < 10.0f) ? 1 : 2;

      // Differential fills per species-pair
      if (idx1 >= 0 && idx2 >= 0) {
        int a = std::min(idx1, idx2), b = std::max(idx1, idx2);
        // etaGap bins
        p_delta_pair_etaGap[a][b][etaBin]->Fill(centrality, delta);
        p_gamma_pair_etaGap[a][b][etaBin]->Fill(centrality, gamma);
        // sumPt bins
        p_delta_pair_sumPt[a][b][sumBin]->Fill(centrality, delta);
        p_gamma_pair_sumPt[a][b][sumBin]->Fill(centrality, gamma);
        // this-cent
        p_delta_pair_sumPt_thisCent[a][b]->Fill(sumPt, delta);
        p_gamma_pair_sumPt_thisCent[a][b]->Fill(sumPt, gamma);
        p_delta_pair_etaGap_thisCent[a][b]->Fill(eta_gap, delta);
        p_gamma_pair_etaGap_thisCent[a][b]->Fill(eta_gap, gamma);
      }

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

void AnalyzerCVE::Write(const std::string &filename) const {
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
  // Pair count histograms
  h_pair_tot_etaGap->Write();
  h_pair_tot_sumPt->Write();
  h_pair_lbc_etaGap->Write();
  h_pair_lbc_sumPt->Write();
  // Species-pair matrices
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      p_gamma_pair[i][j]->Write();
      p_delta_pair[i][j]->Write();
    }
  }
  // Write species-pair differential arrays
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        p_delta_pair_etaGap[i][j][k]->Write();
        p_gamma_pair_etaGap[i][j][k]->Write();
      }
      for (int k = 0; k < 3; ++k) {
        p_delta_pair_sumPt[i][j][k]->Write();
        p_gamma_pair_sumPt[i][j][k]->Write();
      }
      p_delta_pair_sumPt_thisCent[i][j]->Write();
      p_gamma_pair_sumPt_thisCent[i][j]->Write();
      p_delta_pair_etaGap_thisCent[i][j]->Write();
      p_gamma_pair_etaGap_thisCent[i][j]->Write();
    }
  }
  f.Close();
}
