#include "AnalyzerRawOutput.h"

void AnalyzerRawOutput::Init() {
  tree_ = new TTree("Events", "Blast-wave generated events");

  // Event-level
  tree_->Branch("centrality", &centrality_);
  tree_->Branch("centBin", &centBin_);
  tree_->Branch("nParticles", &nParticles_);

  // Particle-level vectors
  tree_->Branch("pid", &pid_);
  tree_->Branch("pt", &pt_);
  tree_->Branch("eta", &eta_);
  tree_->Branch("phi", &phi_);
  tree_->Branch("px", &px_);
  tree_->Branch("py", &py_);
  tree_->Branch("pz", &pz_);
  tree_->Branch("energy", &energy_);
  tree_->Branch("charge", &charge_);
  tree_->Branch("serialNumber", &serialNumber_);
  tree_->Branch("serialNumberLBCFriend", &serialNumberLBCFriend_);
}

void AnalyzerRawOutput::Process(const Event& evt) {
  centrality_ = evt.Centrality();
  centBin_ = evt.centBin;

  const auto& particles = evt.GetParticles();
  nParticles_ = static_cast<int>(particles.size());

  pid_.clear();
  pt_.clear();
  eta_.clear();
  phi_.clear();
  px_.clear();
  py_.clear();
  pz_.clear();
  energy_.clear();
  charge_.clear();
  serialNumber_.clear();
  serialNumberLBCFriend_.clear();

  pid_.reserve(nParticles_);
  pt_.reserve(nParticles_);
  eta_.reserve(nParticles_);
  phi_.reserve(nParticles_);
  px_.reserve(nParticles_);
  py_.reserve(nParticles_);
  pz_.reserve(nParticles_);
  energy_.reserve(nParticles_);
  charge_.reserve(nParticles_);
  serialNumber_.reserve(nParticles_);
  serialNumberLBCFriend_.reserve(nParticles_);

  for (const auto& p : particles) {
    pid_.push_back(p.PID());
    pt_.push_back(static_cast<float>(p.Pt()));
    eta_.push_back(static_cast<float>(p.Eta()));
    phi_.push_back(static_cast<float>(p.Phi()));
    px_.push_back(static_cast<float>(p.Px()));
    py_.push_back(static_cast<float>(p.Py()));
    pz_.push_back(static_cast<float>(p.Momentum().Pz()));
    energy_.push_back(static_cast<float>(p.Momentum().E()));
    charge_.push_back(p.Charge());
    serialNumber_.push_back(p.GetSerialNumber());
    serialNumberLBCFriend_.push_back(p.GetSerialNumberLBCFriend());
  }

  tree_->Fill();
}
