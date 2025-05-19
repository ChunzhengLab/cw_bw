// AnalyzerV0.cxx
#include <TFile.h>
#include <TMath.h>

#include <algorithm>

#include "AnalyzerV0.h"

// 构造函数
AnalyzerV0::AnalyzerV0()
    : fEtaGap(0.4),
      p_meanPt_Pos(nullptr),
      p_meanPt_Neg(nullptr),
      p_meanPt_Neg_mul_Pos(nullptr),
      h_mult(nullptr),
      h_mult_pass(nullptr),
      p_meanPt_Pos_thisEvt(nullptr),
      p_meanPt_Neg_thisEvt(nullptr) {
  // 全局 TProfile 和 TH1D 指针留空，延迟到 Init()
  std::fill_n(f_pt_Neg, kNSpecies, nullptr);
  std::fill_n(f_pt_Neg_mul_Pos, kNSpecies, nullptr);
  std::fill_n(p_v2_pt, kNSpecies, nullptr);
  std::fill_n(f_pt_Neg_thisEvt, kNSpecies, nullptr);
}

// 析构函数
AnalyzerV0::~AnalyzerV0() {
  delete p_meanPt_Pos;
  delete p_meanPt_Neg;
  delete p_meanPt_Neg_mul_Pos;
  for (int i = 0; i < kNSpecies; i++) {
    delete f_pt_Neg[i];
    delete f_pt_Neg_mul_Pos[i];
    delete p_v2_pt[i];
    delete f_pt_Neg_thisEvt[i];
  }
  delete h_mult;
  delete h_mult_pass;
  delete p_meanPt_Pos_thisEvt;
  delete p_meanPt_Neg_thisEvt;
}

// 延迟初始化所有 histogram/profile
void AnalyzerV0::Init() {
  int fNPtBins = 21;
  std::vector<double> fPtEdges = {0.,   0.25, 0.5,  0.75, 1.0,  1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
                                  2.75, 3.0,  3.25, 3.5,  3.75, 4.0,  4.5, 5.0,  6.0, 7.0,  10};
  // 全局 meanPt profiles
  p_meanPt_Pos = new TProfile("p_meanPt_Pos", "meanPt (#eta > fEtaGap)", 1, 0, 1);
  p_meanPt_Neg = new TProfile("p_meanPt_Neg", "meanPt (#eta < 0)", 1, 0, 1);
  p_meanPt_Neg_mul_Pos = new TProfile("p_meanPt_Neg_mul_Pos", "meanPtNeg*meanPtPos", 1, 0, 1);

  // QA 多重直方图
  h_mult = new TH1D("h_mult", "Event Multiplicity;mult;Counts", 500, 0, 5000);
  h_mult_pass = new TH1D("h_mult_pass", "Passed Multiplicity;mult;Counts", 500, 0, 5000);

  // PID 与 all-charged 共用的 f(pt) 与 v2(pt)
  for (int i = 0; i < kNSpecies; i++) {
    TString fneg_name = TString::Format("f_pt_Neg_%s", fName[i]);
    TString fmnp_name = TString::Format("f_pt_Neg_mul_Pos_%s", fName[i]);
    f_pt_Neg[i] = new TProfile(fneg_name, fneg_name, fNPtBins, &fPtEdges[0]);
    f_pt_Neg_mul_Pos[i] = new TProfile(fmnp_name, fmnp_name, fNPtBins, &fPtEdges[0]);
    TString v2_name = TString::Format("p_v2_pt_%s", fName[i]);
    p_v2_pt[i] = new TProfile(v2_name, v2_name, fNPtBins, &fPtEdges[0]);
    // 单事件临时容器
    TString tmp_name = TString::Format("f_pt_Neg_thisEvt_%s", fName[i]);
    f_pt_Neg_thisEvt[i] = new TH1D(tmp_name, tmp_name, fNPtBins, &fPtEdges[0]);
  }

  // 单事件 meanPt profiles
  p_meanPt_Pos_thisEvt = new TProfile("p_meanPt_Pos_thisEvt", "", 1, 0, 1);
  p_meanPt_Neg_thisEvt = new TProfile("p_meanPt_Neg_thisEvt", "", 1, 0, 1);
}

// 每个事件调用
void AnalyzerV0::Process(const Event& evt) {
  int fNPtBins = 21;
  std::vector<double> fPtEdges = {0.,   0.25, 0.5,  0.75, 1.0,  1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
                                  2.75, 3.0,  3.25, 3.5,  3.75, 4.0,  4.5, 5.0,  6.0, 7.0,  10};

  if (!p_meanPt_Pos) Init();  // 延迟初始化

  // 重置单事件容器
  p_meanPt_Pos_thisEvt->Reset();
  p_meanPt_Neg_thisEvt->Reset();
  for (int i = 0; i < kNSpecies; i++) { f_pt_Neg_thisEvt[i]->Reset(); }

  // 统计多重
  int mult_total = 0, mult_pass = 0;
  // 遍历所有粒子
  for (const auto& p : evt.GetParticles()) {
    double pt = p.Pt();
    double eta = p.Eta();
    mult_total++;
    if (pt < 0.2 || TMath::Abs(eta) > 0.8) continue;
    mult_pass++;
    if (eta > fEtaGap) {
      p_meanPt_Pos_thisEvt->Fill(0.5, pt);
    } else if (eta < 0.0) {
      p_meanPt_Neg_thisEvt->Fill(0.5, pt);
    }
    // 根据 η 和 PDG 判断是哪类 species
    int idx = 0;  // all charged
    if (eta < 0.0) { f_pt_Neg_thisEvt[idx]->Fill(pt); }
    p_v2_pt[idx]->Fill(pt, TMath::Cos(2 * p.Phi()));
    ////PID
    if (TMath::Abs(p.PID()) == 211)
      idx = 1;  // pion
    else if (TMath::Abs(p.PID()) == 321)
      idx = 2;  // kaon
    else if (TMath::Abs(p.PID()) == 2212)
      idx = 3;  // proton
    // 负 η 区填 f_pt_thisEvt
    if (eta < 0.0) { f_pt_Neg_thisEvt[idx]->Fill(pt); }
    // v2(pt)
    p_v2_pt[idx]->Fill(pt, TMath::Cos(2 * p.Phi()));
  }

  // 填 QA
  h_mult->Fill(mult_total);
  h_mult_pass->Fill(mult_pass);
  if (mult_pass < 1) return;

  // 单事件均值
  double meanPos = p_meanPt_Pos_thisEvt->GetBinContent(1);
  double meanNeg = p_meanPt_Neg_thisEvt->GetBinContent(1);

  // 全局累积
  p_meanPt_Pos->Fill(0.5, meanPos);
  p_meanPt_Neg->Fill(0.5, meanNeg);
  p_meanPt_Neg_mul_Pos->Fill(0.5, meanNeg * meanPos);

  for (int i = 0; i < kNSpecies; i++) {
    // 归一化负侧分布
    double nneg = f_pt_Neg_thisEvt[i]->Integral(1, fNPtBins);
    if (nneg < 1) continue;
    f_pt_Neg_thisEvt[i]->Scale(1.0 / nneg);
    for (int ib = 1; ib <= fNPtBins; ib++) {
      double f_pt_thisbin = f_pt_Neg_thisEvt[i]->GetBinContent(ib);
      double binCenter = f_pt_Neg_thisEvt[i]->GetBinCenter(ib);
      f_pt_Neg[i]->Fill(binCenter, f_pt_thisbin);
      f_pt_Neg_mul_Pos[i]->Fill(binCenter, f_pt_thisbin * meanPos);
    }
  }
}

// 写到文件
void AnalyzerV0::Write(const std::string& outname) const {
  TFile f(outname.c_str(), "RECREATE");
  p_meanPt_Pos->Write();
  p_meanPt_Neg->Write();
  p_meanPt_Neg_mul_Pos->Write();
  h_mult->Write();
  h_mult_pass->Write();
  for (int i = 0; i < kNSpecies; i++) {
    f_pt_Neg[i]->Write();
    f_pt_Neg_mul_Pos[i]->Write();
    p_v2_pt[i]->Write();
  }
  f.Close();
}
