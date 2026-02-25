#include "Config.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <map>

// Trim leading and trailing whitespace
static std::string trim(const std::string& s) {
  size_t start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  size_t end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

// Strip inline comment: everything after '#' that is not inside '[...]'
static std::string stripComment(const std::string& s) {
  bool inBracket = false;
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '[') inBracket = true;
    else if (s[i] == ']') inBracket = false;
    else if (s[i] == '#' && !inBracket) return s.substr(0, i);
  }
  return s;
}

// Parse a bracketed list "[v1, v2, ...]" into vector<double>
static std::vector<double> parseVector(const std::string& s) {
  std::vector<double> vec;
  // Find content between '[' and ']'
  size_t open = s.find('[');
  size_t close = s.find(']');
  if (open == std::string::npos || close == std::string::npos || close <= open)
    throw std::runtime_error("parseVector: malformed array: " + s);
  std::string inner = s.substr(open + 1, close - open - 1);
  std::istringstream iss(inner);
  std::string token;
  while (std::getline(iss, token, ',')) {
    std::string t = trim(token);
    if (!t.empty()) vec.push_back(std::stod(t));
  }
  return vec;
}

Config Config::Load(const std::string& filename) {
  std::ifstream fin(filename.c_str());
  if (!fin.is_open())
    throw std::runtime_error("Config::Load: cannot open file: " + filename);

  // Store all key-value pairs: key -> raw value string
  std::map<std::string, std::string> kvMap;
  std::string line;
  while (std::getline(fin, line)) {
    std::string stripped = trim(stripComment(line));
    if (stripped.empty()) continue;

    size_t colonPos = stripped.find(':');
    if (colonPos == std::string::npos) continue;

    std::string key = trim(stripped.substr(0, colonPos));
    std::string val = trim(stripped.substr(colonPos + 1));
    if (key.empty()) continue;
    kvMap[key] = val;
  }

  // Helper lambdas
  auto require = [&](const std::string& key) -> const std::string& {
    auto it = kvMap.find(key);
    if (it == kvMap.end())
      throw std::runtime_error("Config::Load(" + filename + "): Missing required key: " + key);
    return it->second;
  };

  auto getInt = [&](const std::string& key) -> int {
    return std::stoi(require(key));
  };

  auto getDouble = [&](const std::string& key) -> double {
    return std::stod(require(key));
  };

  auto getVec = [&](const std::string& key) -> std::vector<double> {
    return parseVector(require(key));
  };

  Config cfg;

  cfg.nEvents = getInt("nEvents");

  cfg.mu5TeV     = getVec("mu5TeV");
  cfg.KinCutRatio = getVec("KinCutRatio");
  cfg.NBDLow     = getVec("NBDLow");
  cfg.NBDHigh    = getVec("NBDHigh");
  cfg.NBDSigma   = getVec("NBDSigma");
  cfg.Tkin       = getVec("Tkin");
  cfg.betaT      = getVec("betaT");
  cfg.n          = getVec("n");
  cfg.rho2_p     = getVec("rho2_p");
  cfg.rho2_L     = getVec("rho2_L");
  cfg.Rx         = getVec("Rx");
  cfg.fracLBC    = getVec("fracLBC");

  cfg.ratioProtonLambda = getDouble("ratioProtonLambda");

  if (kvMap.count("ratioProtonInclusive"))
    cfg.ratioProtonInclusive = std::stod(kvMap["ratioProtonInclusive"]);
  else
    cfg.ratioProtonInclusive = 5.0 / 120.0;

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
    check(cfg.betaT, "betaT");
    check(cfg.n, "n");
    check(cfg.rho2_p, "rho2_p");
    check(cfg.rho2_L, "rho2_L");
    check(cfg.Rx, "Rx");
    check(cfg.fracLBC, "fracLBC");
  }

  return cfg;
}
