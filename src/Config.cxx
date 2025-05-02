#include "Config.h"
#include <stdexcept>
#include <yaml-cpp/yaml.h>

Config Config::Load(const std::string &filename) {
  YAML::Node config = YAML::LoadFile(filename);

  Config cfg;

  if (!config["nEvents"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: nEvents");
  cfg.nEvents = config["nEvents"].as<int>();

  if (!config["mu5TeV"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: mu5TeV");
  for (const auto &val : config["mu5TeV"]) {
    cfg.mu5TeV.push_back(val.as<double>());
  }

  if (!config["k5TeV"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: k5TeV");
  for (const auto &val : config["k5TeV"]) {
    cfg.k5TeV.push_back(val.as<double>());
  }

  if (!config["NBDLow"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: NBDLow");
  for (const auto &val : config["NBDLow"]) {
    cfg.NBDLow.push_back(val.as<double>());
  }

  if (!config["NBDHigh"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: NBDHigh");
  for (const auto &val : config["NBDHigh"]) {
    cfg.NBDHigh.push_back(val.as<double>());
  }

  if (!config["fracLBC"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: fracLBC");
  for (const auto &val : config["fracLBC"]) {
    cfg.fracLBC.push_back(val.as<double>());
  }

  if (!config["Tkin"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: Tkin");
  for (const auto &val : config["Tkin"])
    cfg.Tkin.push_back(val.as<double>());

  if (!config["betaT"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: betaT");
  for (const auto &val : config["betaT"])
    cfg.betaT.push_back(val.as<double>());

  if (!config["n"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: n");
  for (const auto &val : config["n"])
    cfg.n.push_back(val.as<double>());

  if (!config["rho2_p"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: rho2_p");
  for (const auto &val : config["rho2_p"])
    cfg.rho2_p.push_back(val.as<double>());

  if (!config["rho2_L"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: rho2_L");
  for (const auto &val : config["rho2_L"])
    cfg.rho2_L.push_back(val.as<double>());

  if (!config["Rx"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: Rx");
  for (const auto &val : config["Rx"])
    cfg.Rx.push_back(val.as<double>());

  // Particle ratio and debug
  if (!config["ratioProtonLambda"])
    throw std::runtime_error("Config::Load(" + filename +
                             "): Missing required node: ratioProtonLambda");
  cfg.ratioProtonLambda = config["ratioProtonLambda"].as<double>();

  // Optional: ratioProtonInclusive
  if (config["ratioProtonInclusive"])
    cfg.ratioProtonInclusive = config["ratioProtonInclusive"].as<double>();
  else
    cfg.ratioProtonInclusive = 5.0 / 120.0;

  // Optional: isDebug
  cfg.isDebug = false;
  if (config["isDebug"])
    cfg.isDebug = config["isDebug"].as<bool>();

  // Verify that all centrality-binned vectors have the same length
  {
    size_t nBins = cfg.NBDLow.size();
    auto check = [&](const std::vector<double> &v, const std::string &name) {
      if (v.size() != nBins)
        throw std::runtime_error("Config::Load(" + filename +
                                 "): vector size mismatch for '" + name + "'");
    };
    check(cfg.mu5TeV, "mu5TeV");
    check(cfg.k5TeV, "k5TeV");
    check(cfg.NBDHigh, "NBDHigh");
    check(cfg.Tkin, "Tkin");
    check(cfg.betaT, "betaT");
    check(cfg.n, "n");
    check(cfg.rho2_p, "rho2_p");
    check(cfg.rho2_L, "rho2_L");
    check(cfg.Rx, "Rx");
  }

  return cfg;
}
