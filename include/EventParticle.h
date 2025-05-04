#ifndef CW_BW_LBC_EVENTPARTICLE_H
#define CW_BW_LBC_EVENTPARTICLE_H

#include <TLorentzVector.h>

#include <vector>

class Particle {
 public:
  static constexpr double kMassProton = 0.938272;
  static constexpr double kMassLambda = 1.115683;
  Particle(int pid, const TLorentzVector& momentum) : pid_(pid), momentum_(momentum) {
  }

  inline int PID() const {
    return pid_;
  }
  inline const TLorentzVector& Momentum() const {
    return momentum_;
  }
  inline double Pt() const {
    return momentum_.Pt();
  }
  inline double Eta() const {
    return momentum_.PseudoRapidity();
  }
  inline double Phi() const {
    return momentum_.Phi();
  }
  inline int Charge() const {
    return (pid_ > 0 ? 1 : -1);
  }

  inline void SetSerialNumber(int sn) {
    serialNumber_ = sn;
  }
  inline int GetSerialNumber() const {
    return serialNumber_;
  }

  inline void SetSerialNumberLBCFriend(int sn) {
    serialNumberLBCFriend_ = sn;
  }
  inline int GetSerialNumberLBCFriend() const {
    return serialNumberLBCFriend_;
  }

 private:
  int pid_;
  TLorentzVector momentum_;
  int serialNumber_{-1};
  int serialNumberLBCFriend_{-1};
};

struct Event {
  float centrality;                 // Event centrality (0–100)
  int centBin;                      // Centrality bin index
  std::vector<Particle> particles;  // List of particles in this event

  inline size_t Multiplicity() const {
    return particles.size();
  }
  inline float Centrality() const {
    return centrality;
  }
  inline const std::vector<Particle>& GetParticles() const {
    return particles;
  }
};

#endif  // CW_BW_LBC_EVENTPARTICLE_H
