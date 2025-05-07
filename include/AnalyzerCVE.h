#ifndef CW_BW_LBC_ANALYZERCVE_H
#define CW_BW_LBC_ANALYZERCVE_H

#include <TH1D.h>
#include <TProfile.h>
#include <TString.h>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>

#include "EventParticle.h"

// Custom hash for std::pair<int,int> on platforms where std::hash<pair> is not provided.
struct PairHash {
  size_t operator()(const std::pair<int, int>& p) const noexcept {
    uint64_t high = static_cast<uint32_t>(p.first);
    uint64_t low = static_cast<uint32_t>(p.second);
    return (high << 32) ^ low;
  }
};

/// AnalyzerCVE: fills physics observables into histograms and profiles
class AnalyzerCVE {
 public:
  AnalyzerCVE();
  ~AnalyzerCVE();

  /// Initialize all histograms and profiles
  void Init();

  /// Process a single event and fill observables
  void Process(const Event& evt);

  /// Write all histograms to a ROOT file
  void Write(const std::string& filename = "") const;

 private:
  // Single-particle spectra
  TH1D* h_pt_lambda;
  TH1D* h_pt_proton;
  TH1D* h_eta_lambda;
  TH1D* h_eta_proton;
  TH1D* h_phi_lambda;
  TH1D* h_phi_proton;

  // Flow (v2) vs pT
  TProfile* p_v2_pt_lambda;
  TProfile* p_v2_pt_proton;

  // Flow (v2) vs sumPt and EtaGap
  TProfile* p_v2_lambda_sumPt;
  TProfile* p_v2_proton_sumPt;
  TProfile* p_v2_lambda_etaGap;
  TProfile* p_v2_proton_etaGap;

  // Profile maps for pair correlations
  using ProfileMap = std::unordered_map<std::pair<int, int>, TProfile*, PairHash>;

  // Coarse centrality profiles
  ProfileMap p_gamma;
  ProfileMap p_delta;

  ProfileMap p_gamma_LBCFriend;
  ProfileMap p_delta_LBCFriend;

  // Differential by sumPt (3 bins)
  ProfileMap p_delta_sumPt[3];
  ProfileMap p_gamma_sumPt[3];

  // Differential by etaGap (4 bins)
  ProfileMap p_delta_etaGap[4];
  ProfileMap p_gamma_etaGap[4];

  // This-centrality profiles vs sumPt and vs etaGap
  ProfileMap p_delta_vs_sumPt;
  ProfileMap p_gamma_vs_sumPt;
  ProfileMap p_delta_vs_etaGap;
  ProfileMap p_gamma_vs_etaGap;
};

#endif  // CW_BW_LBC_ANALYZERCVE_H
