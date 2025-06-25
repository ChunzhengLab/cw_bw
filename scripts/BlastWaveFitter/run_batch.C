// 工作流程：
// 1) 按给定中心度码(如"2030")拼出图名：
//       spec_proton_2030, spec_lambda_2030,
//       v2_proton_2030,   v2_lambda_2030
// 2) 先到 spec_v2.root 取图；若缺失则到 spec_lambda_infered.root 补。
// 3) 全部找到后  SetDataGraphs(...) → BlastWaveFitter().
// -----------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "BlastWaveFitter.C"
#include "TFile.h"
#include "TGraphAsymmErrors.h"

extern void SetSaveName(const char*);
extern void SetCentralityLabel(const char*);

void run(const char* centCode = "3040") {
  // 编译并载入拟合器（只做一次）
  // gROOT->ProcessLine(".L BlastWaveFitter.C+");

  // 取图
  auto fetch = [](TFile* f, const std::string& name) -> TGraphAsymmErrors* { return f ? dynamic_cast<TGraphAsymmErrors*>(f->Get(name.c_str())) : nullptr; };

  std::string code(centCode);
  std::string gSpec_proton_name = "spec_proton_" + code;
  std::string gSpec_kaon_name = "spec_kaon_" + code;
  std::string gSpec_pion_name = "spec_pion_" + code;

  std::string gV2_proton_name = "v2_proton_" + code;
  std::string gV2_kaon_name = "v2_kaon_" + code;
  std::string gV2_pion_name = "v2_pion_" + code;

  // -------- 先试主文件 --------
  TFile* fMain = TFile::Open("spec_v2.root");
  std::cout << "[run] Trying to load graphs from spec_v2.root ..." << std::endl;
  if (!fMain || fMain->IsZombie()) std::cerr << "Warning: cannot open spec_v2.root\n";

  TGraphAsymmErrors* gSpec_proton = fetch(fMain, gSpec_proton_name);
  TGraphAsymmErrors* gSpec_pion = fetch(fMain, gSpec_pion_name);
  TGraphAsymmErrors* gSpec_kaon = fetch(fMain, gSpec_kaon_name);
  TGraphAsymmErrors* gV2_proton = fetch(fMain, gV2_proton_name);
  TGraphAsymmErrors* gV2_pion = fetch(fMain, gV2_pion_name);
  TGraphAsymmErrors* gV2_kaon = fetch(fMain, gV2_kaon_name);

  // -------- 检查 --------
  if (!gSpec_proton || !gSpec_pion || !gSpec_kaon || !gV2_proton || !gV2_pion || !gV2_kaon) {
    std::cerr << "ERROR: graphs for " << code << " not found in either file.\n";
    return;
  }

  // Set output filename: BlastWaveFit_<cent>.pdf
  std::string outname = "BlastWaveFit_" + code + ".pdf";
  SetSaveName(outname.c_str());
  std::string centLabel = code.substr(0, 2) + "-" + code.substr(2, 2) + "%";
  SetCentralityLabel(centLabel.c_str());
  SetDataGraphs(gSpec_proton, gSpec_pion, gSpec_kaon, gV2_proton, gV2_pion, gV2_kaon);
  BlastWaveFitter();
}

// -----------------------------------------------------------
// Batch driver: fit several centrality windows in one go
void run_batch() {
  // const char* codes[] = {"1020", "2030", "3040", "4050", "5060"};
  const char* codes[] = {"3040"};

  const int Ncodes = sizeof(codes) / sizeof(codes[0]);

  std::vector<std::string> centName;
  std::vector<BWFitResult> res;

  for (int i = 0; i < Ncodes; ++i) {
    std::cout << "\n===========================\n";
    std::cout << "[run_batch] Fitting centrality " << codes[i] << " ...\n";
    run(codes[i]);                      // reuse the single‑centrality driver
    res.push_back(GetLastFitResult());  // harvest summary
    centName.emplace_back(codes[i]);
  }

  std::ofstream fout("BlastWaveFitResults.csv");
  fout << "cent,fit_ok,chi2,ndf,"
          "betaT,Tkin,n_flow,rho2_p,rho2_L,A_p,A_L,Rx\n";
  for (size_t i = 0; i < res.size(); ++i) {
    const auto& r = res[i];
    fout << centName[i] << ',' << (r.valid ? 1 : 0) << ',' << r.chi2 << ',' << r.ndf;
    for (int p = 0; p < 8; ++p) fout << ',' << r.params[p];
    fout << '\n';
  }
  std::cout << "[run_batch] Done. Results written to BlastWaveFitResults.csv\n";
}
