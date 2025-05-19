#include "Config.h"

#include <yaml-cpp/yaml.h>

#include <stdexcept>

Config Config::Load(const std::string& filename) {
  YAML::Node config = YAML::LoadFile(filename);

  Config cfg;

  if (!config["nEvents"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: nEvents");
  cfg.nEvents = config["nEvents"].as<int>();

  if (!config["mu5TeV"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: mu5TeV");
  for (const auto& val : config["mu5TeV"]) { cfg.mu5TeV.push_back(val.as<double>()); }

  if (!config["KinCutRatio"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: KinCutRatio");
  for (const auto& val : config["KinCutRatio"]) { cfg.KinCutRatio.push_back(val.as<double>()); }

  if (!config["NBDLow"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: NBDLow");
  for (const auto& val : config["NBDLow"]) { cfg.NBDLow.push_back(val.as<double>()); }

  if (!config["NBDHigh"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: NBDHigh");
  for (const auto& val : config["NBDHigh"]) { cfg.NBDHigh.push_back(val.as<double>()); }

  if (!config["NBDSigma"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: NBDSigma");
  for (const auto& val : config["NBDSigma"]) { cfg.NBDSigma.push_back(val.as<double>()); }

  if (!config["Tkin"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: Tkin");
  for (const auto& val : config["Tkin"]) cfg.Tkin.push_back(val.as<double>());

  if (!config["TkinSigma"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: TkinSigma");
  for (const auto& val : config["TkinSigma"]) cfg.TkinSigma.push_back(val.as<double>());

  if (!config["betaT"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: betaT");
  for (const auto& val : config["betaT"]) cfg.betaT.push_back(val.as<double>());

  if (!config["betaTSigma"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: betaTSigma");
  for (const auto& val : config["betaTSigma"]) cfg.betaTSigma.push_back(val.as<double>());

  if (!config["n"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: n");
  for (const auto& val : config["n"]) cfg.n.push_back(val.as<double>());

  if (!config["rho2_p"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: rho2_p");
  for (const auto& val : config["rho2_p"]) cfg.rho2_p.push_back(val.as<double>());

  if (!config["rho2_pion"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: rho2_pion");
  for (const auto& val : config["rho2_pion"]) cfg.rho2_pion.push_back(val.as<double>());

  if (!config["rho2_kaon"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: rho2_kaon");
  for (const auto& val : config["rho2_kaon"]) cfg.rho2_kaon.push_back(val.as<double>());

  if (!config["Rx"]) throw std::runtime_error("Config::Load(" + filename + "): Missing required node: Rx");
  for (const auto& val : config["Rx"]) cfg.Rx.push_back(val.as<double>());

  // Particle ratio and debug
  if (!config["ratioProtonPion"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: ratioProtonPion");
  cfg.ratioProtonPion = config["ratioProtonPion"].as<double>();

  if (!config["ratioKaonPion"])
    throw std::runtime_error("Config::Load(" + filename + "): Missing required node: ratioKaonPion");
  cfg.ratioKaonPion = config["ratioKaonPion"].as<double>();

  // Verify that all centrality-binned vectors have the same length
  {
    size_t nBins = cfg.NBDLow.size();
    auto check = [&](const std::vector<double>& v, const std::string& name) {
      if (v.size() != nBins)
        throw std::runtime_error("Config::Load(" + filename + "): vector size mismatch for '" + name + "'");
    };
    check(cfg.mu5TeV, "mu5TeV");
    check(cfg.KinCutRatio, "KinCutRatio");
    check(cfg.NBDHigh, "NBDHigh");
    check(cfg.NBDSigma, "NBDSigma");
    check(cfg.Tkin, "Tkin");
    check(cfg.TkinSigma, "TkinSigma");
    check(cfg.betaT, "betaT");
    check(cfg.betaTSigma, "betaTSigma");
    check(cfg.n, "n");
    check(cfg.rho2_p, "rho2_p");
    check(cfg.rho2_pion, "rho2_pion");
    check(cfg.rho2_kaon, "rho2_kaon");
    check(cfg.Rx, "Rx");
  }

  return cfg;
}
