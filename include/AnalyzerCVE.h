#ifndef CW_BW_LBC_ANALYZERCVE_H
#define CW_BW_LBC_ANALYZERCVE_H

#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TProfile.h>
#include "EventParticle.h"

// AnalyzerCVE: fills physics observables into histograms and profiles
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

    // Species-pair differential by EtaGap and sumPt
    // Indexed [species1][species2][bin]
    TProfile* p_delta_pair_etaGap[4][4][4];
    TProfile* p_gamma_pair_etaGap[4][4][4];
    TProfile* p_delta_pair_sumPt[4][4][3];
    TProfile* p_gamma_pair_sumPt[4][4][3];

    // Species-pair differential per-centrality (this-cent)
    // Indexed [species1][species2]
    TProfile* p_delta_pair_sumPt_thisCent[4][4];
    TProfile* p_gamma_pair_sumPt_thisCent[4][4];
    TProfile* p_delta_pair_etaGap_thisCent[4][4];
    TProfile* p_gamma_pair_etaGap_thisCent[4][4];

    // Flow (v2) vs sumPt and EtaGap
    TProfile* p_v2_lambda_sumPt;
    TProfile* p_v2_proton_sumPt;
    TProfile* p_v2_lambda_etaGap;
    TProfile* p_v2_proton_etaGap;

    // Pair count histograms
    TH1D* h_pair_tot_etaGap;
    TH1D* h_pair_tot_sumPt;
    TH1D* h_pair_lbc_etaGap;
    TH1D* h_pair_lbc_sumPt;

    // Species-pair correlation matrix (p, Lambda, anti-p, anti-Lambda)
    TProfile* p_gamma_pair[4][4];
    TProfile* p_delta_pair[4][4];
};

#endif // CW_BW_LBC_ANALYZERCVE_H
