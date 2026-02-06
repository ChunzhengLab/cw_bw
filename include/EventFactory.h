#ifndef CW_BW_LBC_EVENTFACTORY_H
#define CW_BW_LBC_EVENTFACTORY_H

#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

#include <array>
#include <memory>
#include <random>
#include <vector>

#include "Config.h"
#include "EventParticle.h"

/// 生成模式：pikp = π/K/p, plambda = p/Λ (原有行为)
enum class GenerationMode { kPiKP, kPLambda };

/// EventFactory: 根据配置和中心度生成完整的 Event 对象
class EventFactory {
 public:
  /// 构造：传入全局配置和生成模式
  explicit EventFactory(const Config& cfg, GenerationMode mode = GenerationMode::kPLambda);

  /// 生成并返回一个事件，内部完成分箱和粒子构建
  Event GetEvent(float centrality);

 private:
  // Configuration and random number generator
  const Config& cfg_;
  GenerationMode mode_;
  mutable std::mt19937_64 rng_{std::random_device{}()};
  // Uniform distributions for random sampling
  mutable std::uniform_real_distribution<double> dist01_{0.0, 1.0};
  mutable std::uniform_real_distribution<double> dist2pi_{0.0, 2.0 * M_PI};

  // Prebuilt energy samplers (Maxwell–Jüttner) per bin/species
  std::vector<std::unique_ptr<TF1>> fE_pion_;
  std::vector<std::unique_ptr<TF1>> fE_kaon_;
  std::vector<std::unique_ptr<TF1>> fE_proton_;
  std::vector<std::unique_ptr<TF1>> fE_lambda_;

  void DetermineCentBin(Event& evt) const;
  void BuildParticles(Event& evt);
  inline void SeedEmissionPoint(float& x, float& y, int bin) const;
  int GivePidBasedOnRatio(float ratio) const;
  static double GetMassFromPid(int pid);
  TVector3 GetLocalBoostVector(float x, float y, int bin, int pid) const;
  double SampleEnergy(int bin, int pid) const;
  inline std::array<double, 3> SampleDirection() const;
};

#endif  // CW_BW_LBC_EVENTFACTORY_H
