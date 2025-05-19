#ifndef ANALYZER_V0_H
#define ANALYZER_V0_H

#include <TH1D.h>
#include <TProfile.h>

#include <string>
#include <vector>

#include "EventParticle.h"

class AnalyzerV0 {
 public:
  // 构造：传入 pt bin 数和边界
  AnalyzerV0();
  ~AnalyzerV0();

  // 每个事件调用
  void Process(const Event& evt);
  // 最后写文件
  void Write(const std::string& outname) const;
  // 延迟初始化
  void Init();

 private:
  double fEtaGap;
  // 4 种“粒子”：0=all charged,1=pi,2=K,3=p
  static constexpr int kNSpecies = 4;
  const char* fName[kNSpecies] = {"charge", "pion", "kaon", "proton"};

  // ---- 全局 & PID 共用 ----
  // meanPt(+) 与 meanPt(-)
  TProfile* p_meanPt_Pos;
  TProfile* p_meanPt_Neg;
  TProfile* p_meanPt_Neg_mul_Pos;

  // f(pt)(-) 与 f(pt)*meanPt(+)
  TProfile* f_pt_Neg[kNSpecies];
  TProfile* f_pt_Neg_mul_Pos[kNSpecies];

  TH1D* h_mult;
  TH1D* h_mult_pass;

  // v2(pt) = <cos2φ>
  TProfile* p_v2_pt[kNSpecies];

  // 单事件容器
  TH1D* f_pt_Neg_thisEvt[kNSpecies];
  TProfile* p_meanPt_Pos_thisEvt;
  TProfile* p_meanPt_Neg_thisEvt;
};

#endif  // ANALYZER_V0_H
