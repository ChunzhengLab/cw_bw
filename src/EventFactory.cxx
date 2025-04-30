#include <memory>
#include <TF1.h>
#include "EventFactory.h"
#include "Config.h"
#include <cmath>
#include <TF1.h>
#include <TString.h>
#include "EventParticle.h"

// Constants for maximum momentum and default Ry
static constexpr double kMaxMomentum = 10.0;  // GeV/c, upper p* cutoff
static constexpr double kRyDefault    = 10.0;  // fm, fixed Ry

// Maxwell–Jüttner kernel for TF1: parameters[0]=mass, [1]=temp
static double MaxwellJuttnerKernel(double *x, double *par) {
    double E    = x[0];
    double mass = par[0];
    double temp = par[1];
    double p    = std::sqrt(E*E - mass*mass);
    return p * E * std::exp(-E / temp);
}

// Compute the local flow angle φ_b given emission angle φ_s and ellipse axes Rx,Ry
inline double ComputePhiBFromPhiS(double phi_s, double Rx, double Ry) {
    // tan φ_b = (Ry²/Rx²) tan φ_s
    double tanb = (Ry*Ry)/(Rx*Rx) * std::tan(phi_s);
    double phi_b = std::atan(tanb);
    // adjust into correct quadrant
    if (phi_s > M_PI/2 && phi_s < 3*M_PI/2) phi_b += M_PI;
    return phi_b;
}

EventFactory::EventFactory(const Config& cfg)
  : cfg_(cfg), rnd_(std::random_device{}()) 
{
    int nBins = cfg_.Tkin.size();
    fE_proton_.reserve(nBins);
    fE_lambda_.reserve(nBins);
    const double pMax = kMaxMomentum;
    for (int i = 0; i < nBins; ++i) {
        double mass_p = Particle::kMassProton;
        double mass_L = Particle::kMassLambda;
        double temp   = cfg_.Tkin[i];
        double Emax_p = std::sqrt(mass_p*mass_p + pMax*pMax);
        auto fp = std::make_unique<TF1>(Form("fE_p_%d", i),
                                        MaxwellJuttnerKernel,
                                        mass_p, Emax_p, 2);
        fp->SetParameter(0, mass_p);
        fp->SetParameter(1, temp);
        fE_proton_.push_back(std::move(fp));
        double Emax_L = std::sqrt(mass_L*mass_L + pMax*pMax);
        auto fL = std::make_unique<TF1>(Form("fE_L_%d", i),
                                        MaxwellJuttnerKernel,
                                        mass_L, Emax_L, 2);
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
    evt.centBin = 0;
    for (int i = 0; i < static_cast<int>(cfg_.NBDLow.size()); ++i) {
        if (evt.centrality >= cfg_.NBDLow[i] && evt.centrality < cfg_.NBDHigh[i]) {
            evt.centBin = i;
            break;
        }
    }
}

// 构建该事件的所有粒子
void EventFactory::BuildParticles(Event& evt) {
    int bin = evt.centBin;
    int multiplicity = static_cast<int>(cfg_.mu5TeV[bin]);
    evt.particles.clear();
    evt.particles.reserve(multiplicity);
    int nextSN = 0;

    for (int i = 0; i < multiplicity; ++i) {
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
        double E  = SampleEnergy(bin, pid);
        auto   dir = SampleDirection();
        TLorentzVector momentum(dir[0]*E, dir[1]*E, dir[2]*E, E);
        momentum.Boost(boost);

        Particle p(pid, momentum);
        p.SetSerialNumber(sn);

        // 6) LBC 配对 (inline probability)
        if (rnd_.Rndm() < cfg_.fracLBC[bin]) {
            int sn_friend = ++nextSN;
            int pid_friend = -pid;
            double E_friend    = SampleEnergy(bin, pid_friend);
            auto   dir_friend  = SampleDirection();
            TLorentzVector mom_friend(dir_friend[0]*E_friend,
                                      dir_friend[1]*E_friend,
                                      dir_friend[2]*E_friend,
                                      E_friend);
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
}

// 抽样发射点 (x,y)
void EventFactory::SeedEmissionPoint(float& x, float& y, int bin) const {
    // Elliptical uniform area sampling: r_norm = sqrt(u), phi_s uniform
    float u      = rnd_.Rndm();
    float r_norm = std::sqrt(u);
    float phi_s  = 2 * static_cast<float>(M_PI) * rnd_.Rndm();
    // Position inside ellipse with semi-axes Rx and Ry
    x = r_norm * static_cast<float>(cfg_.Rx[bin]) * std::cos(phi_s);
    y = r_norm * static_cast<float>(kRyDefault)   * std::sin(phi_s);
}

// 计算流场 boost 向量
TVector3 EventFactory::GetLocalBoostVector(float x, float y, int bin, int pid) const {
    // Elliptical coordinates & normalized radius
    double Rx = cfg_.Rx[bin];
    double Ry = kRyDefault;
    // normalized radius in ellipse
    double r_norm = std::hypot(x/Rx, y/Ry);
    // emission angle in ellipse coordinate
    double phi_s = std::atan2(y/Ry, x/Rx);

    // flow rapidity parameters
    double rho0 = std::atanh(cfg_.betaT[bin]);
    // choose anisotropy based on particle type
    double rho2 = (std::abs(pid) == 2212 ? cfg_.rho2_p[bin] : cfg_.rho2_L[bin]);
    double phi_b = ComputePhiBFromPhiS(phi_s, Rx, Ry);
    double rhob = std::pow(r_norm, cfg_.n[bin]) * (rho0 + rho2 * std::cos(2 * phi_b));

    // longitudinal pseudorapidity sampling
    double cstheta = 2.0*(rnd_.Rndm() - 0.5);
    double thetas = std::acos(cstheta);
    double etas   = -std::log(std::tan(thetas/2.0));

    // construct fluid four-velocity and return its boost vector
    TLorentzVector u;
    u.SetXYZT(std::sinh(rhob)*std::cos(phi_b),
              std::sinh(rhob)*std::sin(phi_b),
              std::cosh(rhob)*std::sinh(etas),
              std::cosh(rhob)*std::cosh(etas));
    return u.BoostVector();
}

// 根据比例抽 PID
int EventFactory::GivePidBasedOnRatio(float ratio) const {
    float p = rnd_.Rndm();
    int pid = (p < ratio/(ratio+1.0f) ? 2212 : 3122);
    if (rnd_.Rndm() < 0.5f) pid = -pid;
    return pid;
}

double EventFactory::SampleEnergy(int bin, int pid) const {
    // choose sampler based on species
    auto &sampler = (std::abs(pid) == 2212 ? fE_proton_[bin] : fE_lambda_[bin]);
    return sampler->GetRandom();
}

// 均匀方向抽样
std::array<double,3> EventFactory::SampleDirection() const {
    double cosTheta = 2.0 * rnd_.Rndm() - 1.0;
    double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);
    double phi      = 2.0 * M_PI * rnd_.Rndm();
    return { sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta };
}
