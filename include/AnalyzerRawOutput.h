#ifndef CW_BW_LBC_ANALYZERRAWOUTPUT_H
#define CW_BW_LBC_ANALYZERRAWOUTPUT_H

#include <TTree.h>

#include <vector>

#include "EventParticle.h"

/// AnalyzerRawOutput: 将事件的原始粒子信息写入 ROOT TTree
/// 不管理 TFile，由调用方负责打开/关闭文件。
/// TTree 创建时自动注册到当前 TFile。
class AnalyzerRawOutput {
 public:
  AnalyzerRawOutput() = default;
  ~AnalyzerRawOutput() = default;

  /// 创建 TTree 及分支（需在 TFile 已打开后调用）
  void Init();

  /// 将一个事件填入 TTree
  void Process(const Event& evt);

 private:
  TTree* tree_ = nullptr;  // owned by current TFile

  // Event-level branches
  float centrality_{};
  int centBin_{};
  int nParticles_{};

  // Particle-level branches (per event, variable length)
  std::vector<int> pid_;
  std::vector<float> pt_;
  std::vector<float> eta_;
  std::vector<float> phi_;
  std::vector<float> px_;
  std::vector<float> py_;
  std::vector<float> pz_;
  std::vector<float> energy_;
  std::vector<int> charge_;
  std::vector<int> serialNumber_;
  std::vector<int> serialNumberLBCFriend_;
};

#endif  // CW_BW_LBC_ANALYZERRAWOUTPUT_H
