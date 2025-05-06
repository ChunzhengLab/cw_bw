#include "EventFactory.h"

#include <TString.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

#include "Config.h"
#include "EventParticle.h"

// Constants for maximum momentum and default Ry
static constexpr double kMaxMomentum = 10.0;  // GeV/c, upper p* cutoff
static constexpr double kRyDefault = 10.0;    // fm, fixed Ry

// Maxwell–Jüttner kernel for TF1: parameters[0]=mass, [1]=temp
static double MaxwellJuttnerKernel(double* x, double* par) {
  double E = x[0];
  double mass = par[0];
  double temp = par[1];
  double p = std::sqrt(E * E - mass * mass);
  return p * E * std::exp(-E / temp);
}

// Compute the local flow angle φ_b given emission angle φ_s and ellipse axes
// Rx,Ry
inline double ComputePhiBFromPhiS(double phi_s, double Rx, double Ry) {
  // tan φ_b = (Ry²/Rx²) tan φ_s
  double y = std::sin(phi_s) / (Ry * Ry);
  double x = std::cos(phi_s) / (Rx * Rx);
  double phi_b = std::atan2(y, x);
  return phi_b;
}

EventFactory::EventFactory(const Config& cfg) : cfg_(cfg) {
  int nBins = cfg_.Tkin.size();
  fE_proton_.reserve(nBins);
  fE_lambda_.reserve(nBins);
  const double pMax = kMaxMomentum;
  for (int i = 0; i < nBins; ++i) {
    double mass_p = Particle::kMassProton;
    double mass_L = Particle::kMassLambda;
    double temp = cfg_.Tkin[i];
    double Emax_p = std::sqrt(mass_p * mass_p + pMax * pMax);
    auto fp = std::make_unique<TF1>(Form("fE_p_%d", i), MaxwellJuttnerKernel, mass_p, Emax_p, 2);
    fp->SetParameter(0, mass_p);
    fp->SetParameter(1, temp);
    fE_proton_.push_back(std::move(fp));
    double Emax_L = std::sqrt(mass_L * mass_L + pMax * pMax);
    auto fL = std::make_unique<TF1>(Form("fE_L_%d", i), MaxwellJuttnerKernel, mass_L, Emax_L, 2);
    fL->SetParameter(0, mass_L);
    fL->SetParameter(1, temp);
    fE_lambda_.push_back(std::move(fL));
  }
}

// Public: 生成事件
Event EventFactory::GetEvent(float centrality) {
  Event evt;
  evt.centrality = centrality;
  DetermineCentBin(evt);
  BuildParticles(evt);
  return evt;
}

// 根据中心度设定分箱
void EventFactory::DetermineCentBin(Event& evt) const {
  float cent = evt.centrality;
  if (cent < 10.0f || cent >= 60.0f) { throw std::runtime_error("Unsupported centrality range: must be in [10,60)%"); }

  if (cent >= 10.0f && cent < 20.0f)
    evt.centBin = 0;
  else if (cent >= 20.0f && cent < 30.0f)
    evt.centBin = 1;
  else if (cent >= 30.0f && cent < 40.0f)
    evt.centBin = 2;
  else if (cent >= 40.0f && cent < 50.0f)
    evt.centBin = 3;
  else if (cent >= 50.0f && cent < 60.0f)
    evt.centBin = 4;
}

// 构建该事件的所有粒子
void EventFactory::BuildParticles(Event& evt) {
  int bin = evt.centBin;
  // 使用 STL 负二项分布和 uniform 分布
  // 考虑 pT/η 裁剪后的保留比例
  double mu_pre = cfg_.mu5TeV[bin];
  double mu = mu_pre / cfg_.KinCutRatio[bin];
  double sigma = cfg_.NBDSigma[bin];
  double variance = sigma * sigma;
  double p = mu / variance;
  int r_int = static_cast<int>(std::round(mu * mu / (variance - mu)));
  std::negative_binomial_distribution<int> nbdist(r_int, p);

  unsigned int multiplicity;
  do {
    multiplicity = static_cast<unsigned int>(nbdist(rng_));
  } while (multiplicity < static_cast<unsigned int>(cfg_.NBDLow[bin]) ||
           multiplicity > static_cast<unsigned int>(cfg_.NBDHigh[bin]));

  // 只生成 Lambda 和 Proton，修正多重数
  multiplicity = static_cast<unsigned int>(multiplicity * cfg_.ratioProtonInclusive * (1.0 + cfg_.ratioProtonLambda) /
                                           cfg_.ratioProtonLambda);

  evt.particles.clear();
  evt.particles.reserve(2 * multiplicity);
  int nextSN = 0;

  for (unsigned int i = 0; i < multiplicity; ++i) {
    int sn = ++nextSN;

    // 1) 发射点抽样
    float x, y;
    SeedEmissionPoint(x, y, bin);

    // 2) PID 抽样
    int pid = GivePidBasedOnRatio(cfg_.ratioProtonLambda);

    // 3) 计算流场 boost (species-dependent ρ2)
    TVector3 boost = GetLocalBoostVector(x, y, bin, pid);

    // 4) 质量/温度
    double mass = (std::abs(pid) == 2212 ? Particle::kMassProton : Particle::kMassLambda);
    double temp = cfg_.Tkin[bin];

    // 5) 动量抽样
    double E = SampleEnergy(bin, pid);
    // compute momentum magnitude p* from energy and mass
    double p_star = std::sqrt(E * E - mass * mass);
    auto dir = SampleDirection();
    TLorentzVector momentum(dir[0] * p_star, dir[1] * p_star, dir[2] * p_star, E);
    momentum.Boost(boost);

    Particle p(pid, momentum);
    p.SetSerialNumber(sn);

    // 6) LBC 配对 (inline probability)
    if (dist01_(rng_) < cfg_.fracLBC[bin]) {
      int sn_friend = ++nextSN;
      int pid_friend;
      if (pid == 2212)       pid_friend = -3122;
      else if (pid == -2212) pid_friend =  3122;
      else if (pid == 3122)  pid_friend = -2212;
      else if (pid == -3122) pid_friend =  2212;
      else pid_friend = -pid;

      double E_friend = SampleEnergy(bin, pid_friend);
      double p_star_friend = std::sqrt(E_friend * E_friend - mass * mass);
      auto dir_friend = SampleDirection();
      TLorentzVector mom_friend(dir_friend[0] * p_star_friend, dir_friend[1] * p_star_friend,
                                dir_friend[2] * p_star_friend, E_friend);
      mom_friend.Boost(boost);

      Particle p_friend(pid_friend, mom_friend);
      p_friend.SetSerialNumber(sn_friend);
      p_friend.SetSerialNumberLBCFriend(p.GetSerialNumber());
      p.SetSerialNumberLBCFriend(p_friend.GetSerialNumber());
      evt.particles.push_back(std::move(p_friend));
    } else {
      p.SetSerialNumberLBCFriend(-1);
    }

    // 7) 添加主粒子
    evt.particles.push_back(std::move(p));
  }
  // 随机打乱粒子顺序后截断至 multiplicity，不再区分主粒子和 LBC 伴侣
  std::shuffle(evt.particles.begin(), evt.particles.end(), rng_);
  if (evt.particles.size() > multiplicity) {
    evt.particles.erase(evt.particles.begin() + multiplicity, evt.particles.end());
    evt.particles.shrink_to_fit();
  }
}

// 抽样发射点 (x,y)
inline void EventFactory::SeedEmissionPoint(float& x, float& y, int bin) const {
  float u = dist01_(rng_);
  float phi_s = dist2pi_(rng_);
  // Position inside ellipse with semi-axes Rx and Ry
  x = std::sqrt(u) * static_cast<float>(cfg_.Rx[bin]) * std::cos(phi_s);
  y = std::sqrt(u) * static_cast<float>(kRyDefault) * std::sin(phi_s);
}

// 计算流场 boost 向量
TVector3 EventFactory::GetLocalBoostVector(float x, float y, int bin, int pid) const {
  // Elliptical coordinates & normalized radius
  double Rx = cfg_.Rx[bin];
  double Ry = kRyDefault;
  // normalized radius in ellipse
  double r_norm = std::hypot(x / Rx, y / Ry);
  // emission angle in ellipse coordinate
  double phi_s = std::atan2(y, x); // SUGGESTED by WAN JIE Thanks!
  // double phi_s = std::atan2(y / Ry, x / Rx);

  // flow rapidity parameters
  double rho0 = std::atanh(cfg_.betaT[bin]);
  // choose anisotropy based on particle type
  double rho2 = (std::abs(pid) == 2212 ? cfg_.rho2_p[bin] : cfg_.rho2_L[bin]);
  double phi_b = ComputePhiBFromPhiS(phi_s, Rx, Ry);
  double rhob = std::pow(r_norm, cfg_.n[bin]) * (rho0 + rho2 * std::cos(2 * phi_b));

  // longitudinal pseudorapidity sampling
  double etas = std::atanh(dist01_(rng_) * 2.0 - 1.0);

  // construct fluid four-velocity and return its boost vector
  TLorentzVector u;
  u.SetXYZT(std::sinh(rhob) * std::cos(phi_b), std::sinh(rhob) * std::sin(phi_b), std::cosh(rhob) * std::sinh(etas),
            std::cosh(rhob) * std::cosh(etas));
  return u.BoostVector();
}

// 根据比例抽 PID
inline int EventFactory::GivePidBasedOnRatio(float ratio) const {
  double p = dist01_(rng_);
  int pid = (p < ratio / (ratio + 1.0f) ? 2212 : 3122);
  if (dist01_(rng_) < 0.5) pid = -pid;
  return pid;
}

double EventFactory::SampleEnergy(int bin, int pid) const {
  // choose sampler based on species
  auto& sampler = (std::abs(pid) == 2212 ? fE_proton_[bin] : fE_lambda_[bin]);
  return sampler->GetRandom();
}

// 均匀方向抽样
inline std::array<double, 3> EventFactory::SampleDirection() const {
  double cosTheta = dist01_(rng_) * 2.0 - 1.0;
  double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
  double phi = dist2pi_(rng_);
  return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
}
